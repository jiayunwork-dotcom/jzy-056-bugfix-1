package geom

import (
	"math"
	"sort"
)

// IsCollinear reports whether every point of the set lies on one straight
// line (within the numeric tolerance band). Sets of 0-2 points are collinear.
func IsCollinear(pts []Point) bool {
	if len(pts) < 3 {
		return true
	}
	// Anchor on the two most separated points along the principal axis; a
	// wide baseline makes the degeneracy test numerically strongest.
	i0, i1 := widestPair(pts)
	a, b := pts[i0], pts[i1]
	if a == b {
		return true
	}
	for i, p := range pts {
		if i == i0 || i == i1 {
			continue
		}
		if rawOrient(a, b, p) != On {
			return false
		}
	}
	return true
}

// widestPair returns indices of the pair with the largest axis-aligned spread
// (using whichever of x or y dominates), which is a cheap stand-in for the
// diameter pair.
func widestPair(pts []Point) (int, int) {
	var minX, maxX, minY, maxY int
	for i := 1; i < len(pts); i++ {
		if pts[i].X < pts[minX].X {
			minX = i
		}
		if pts[i].X > pts[maxX].X {
			maxX = i
		}
		if pts[i].Y < pts[minY].Y {
			minY = i
		}
		if pts[i].Y > pts[maxY].Y {
			maxY = i
		}
	}
	if pts[maxX].X-pts[minX].X >= pts[maxY].Y-pts[minY].Y {
		return minX, maxX
	}
	return minY, maxY
}

// ConvexHull returns the vertices of the convex hull in counter-clockwise
// order using Andrew's monotone chain. Points on a hull edge but not at its
// corners are not included (the result are the true hull corners).
func ConvexHull(pts []Point) []Point {
	sorted := make([]Point, len(pts))
	copy(sorted, pts)
	sort.Slice(sorted, func(i, j int) bool {
		if sorted[i].X != sorted[j].X {
			return sorted[i].X < sorted[j].X
		}
		return sorted[i].Y < sorted[j].Y
	})
	// Exact duplicates must already have been rejected; defensively collapse.
	uniq := sorted[:0]
	for i, p := range sorted {
		if i > 0 && p == sorted[i-1] {
			continue
		}
		uniq = append(uniq, p)
	}
	if len(uniq) <= 1 {
		return append([]Point(nil), uniq...)
	}

	// crossOn: pop while the turn is not a strict left turn; the ON band is
	// treated as collinear so points sitting on a hull edge get dropped.
	cross := func(o, a, b Point) Sign { return rawOrient(o, a, b) }

	lower := make([]Point, 0, len(uniq))
	for _, p := range uniq {
		for len(lower) >= 2 && cross(lower[len(lower)-2], lower[len(lower)-1], p) != Positive {
			lower = lower[:len(lower)-1]
		}
		lower = append(lower, p)
	}
	upper := make([]Point, 0, len(uniq))
	for i := len(uniq) - 1; i >= 0; i-- {
		p := uniq[i]
		for len(upper) >= 2 && cross(upper[len(upper)-2], upper[len(upper)-1], p) != Positive {
			upper = upper[:len(upper)-1]
		}
		upper = append(upper, p)
	}
	// First/last points are duplicated between chains.
	return append(lower[:len(lower)-1], upper[:len(upper)-1]...)
}

// PolygonArea returns the area of a simple polygon given in order (any
// winding; the absolute value is returned).
func PolygonArea(ring []Point) float64 {
	var s float64
	n := len(ring)
	for i := 0; i < n; i++ {
		a := ring[i]
		b := ring[(i+1)%n]
		s += a.X*b.Y - b.X*a.Y
	}
	return math.Abs(s) * 0.5
}

// SignedArea2 returns the doubled signed area of triangle (a,b,c).
func SignedArea2(a, b, c Point) float64 {
	return (b.X-a.X)*(c.Y-a.Y) - (b.Y-a.Y)*(c.X-a.X)
}

// HullClearance returns the smallest distance from a point strictly INSIDE
// the convex hull to any hull edge SEGMENT. Points on the hull itself are
// excluded: a very acute hull corner projects just beyond an adjacent edge's
// endpoint and would otherwise report a meaningless near-zero distance to the
// infinite supporting line. The segment (not infinite-line) distance is used.
func HullClearance(pts []Point, hull []Point) float64 {
	best := math.Inf(1)
	isHull := make(map[Point]struct{}, len(hull))
	for _, h := range hull {
		isHull[h] = struct{}{}
	}
	for _, p := range pts {
		if _, on := isHull[p]; on {
			continue
		}
		for i := 0; i < len(hull); i++ {
			a := hull[i]
			b := hull[(i+1)%len(hull)]
			if d := pointSegmentDistance(p, a, b); d > 0 && d < best {
				best = d
			}
		}
	}
	if !math.IsInf(best, 0) {
		return best
	}
	_, _, d := BoundsOf(hull).Span()
	if d == 0 {
		return 1
	}
	return d
}

// pointSegmentDistance is the Euclidean distance from p to segment [a,b].
func pointSegmentDistance(p, a, b Point) float64 {
	abx := b.X - a.X
	aby := b.Y - a.Y
	len2 := abx*abx + aby*aby
	if len2 == 0 {
		return math.Hypot(p.X-a.X, p.Y-a.Y)
	}
	t := ((p.X-a.X)*abx + (p.Y-a.Y)*aby) / len2
	if t < 0 {
		t = 0
	} else if t > 1 {
		t = 1
	}
	qx := a.X + t*abx
	qy := a.Y + t*aby
	return math.Hypot(p.X-qx, p.Y-qy)
}
