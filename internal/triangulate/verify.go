package triangulate

import (
	"math"

	"github.com/example/delaunay/internal/geom"
)

// MeshReport summarises geometric properties of a finished triangulation.
type MeshReport struct {
	TriangleCount int
	// EmptyCircle is true when no input point lies strictly inside any
	// triangle's circumcircle (points within the tolerance band count as
	// cocircular and pass).
	EmptyCircle bool
	// AllCCW is true when every triangle is strictly counter-clockwise.
	AllCCW bool
	// AreaSum is the sum of unsigned triangle areas.
	AreaSum float64
	// NonOverlapping is true when no two triangle interiors intersect
	// (shared vertices/edges allowed).
	NonOverlapping bool
}

// Verify recomputes, from scratch and with the same tolerant predicates,
// whether tris over pts satisfies the Delaunay empty-circumcircle property
// and basic mesh integrity. The predicate context mirrors the one used while
// building: recentred on pts[0] with the translation-noise radius.
func Verify(pts []geom.Point, tris []Tri) Mesh {
	anchor := pts[0]
	radius := geom.TranslationNoiseRadius(anchor)
	pc := geom.PredCtx{Radius: radius}
	local := make([]geom.Point, len(pts))
	for i := range pts {
		local[i] = pts[i].Sub(anchor)
	}

	r := Mesh{
		TriangleCount:  len(tris),
		EmptyCircle:    true,
		AllCCW:         true,
		NonOverlapping: true,
	}
	for _, t := range tris {
		a, b, c := local[t.A], local[t.B], local[t.C]
		if pc.Orient2d(a, b, c) != geom.Positive {
			r.AllCCW = false
		}
		ar := math.Abs(geom.SignedArea2(a, b, c)) * 0.5
		r.AreaSum += ar
		for p := range local {
			if p == t.A || p == t.B || p == t.C {
				continue
			}
			if pc.InCircle(a, b, c, local[p]) == geom.Positive {
				r.EmptyCircle = false
			}
		}
	}

	// Pairwise interior-intersection check (quadratic; this is a verifier).
	for i := 0; i < len(tris); i++ {
		for j := i + 1; j < len(tris); j++ {
			if trianglesIntersect(local, pc, tris[i], tris[j]) {
				r.NonOverlapping = false
			}
		}
	}
	return r
}

// Mesh is the exported report type (alias keeps the construction concise).
type Mesh = MeshReport

// trianglesIntersect reports whether two triangles' interiors overlap, given
// that a valid mesh only shares full edges or vertices. We check for any
// pair of edges crossing in their interiors.
func trianglesIntersect(pts []geom.Point, pc geom.PredCtx, t1, t2 Tri) bool {
	e1 := [3][2]int{{t1.A, t1.B}, {t1.B, t1.C}, {t1.C, t1.A}}
	e2 := [3][2]int{{t2.A, t2.B}, {t2.B, t2.C}, {t2.C, t2.A}}
	for _, a := range e1 {
		for _, b := range e2 {
			// Shared endpoint -> edges meet legitimately.
			if a[0] == b[0] || a[0] == b[1] || a[1] == b[0] || a[1] == b[1] {
				continue
			}
			if segmentsCross(pc, pts[a[0]], pts[a[1]], pts[b[0]], pts[b[1]]) {
				return true
			}
		}
	}
	return false
}

// segmentsCross: proper intersection of two closed segments with no shared
// endpoint (straddling on both sides).
func segmentsCross(pc geom.PredCtx, p1, p2, p3, p4 geom.Point) bool {
	d1 := pc.Orient2d(p3, p4, p1)
	d2 := pc.Orient2d(p3, p4, p2)
	d3 := pc.Orient2d(p1, p2, p3)
	d4 := pc.Orient2d(p1, p2, p4)
	return ((d1 == geom.Positive && d2 == geom.Negative) || (d1 == geom.Negative && d2 == geom.Positive)) &&
		((d3 == geom.Positive && d4 == geom.Negative) || (d3 == geom.Negative && d4 == geom.Positive))
}
