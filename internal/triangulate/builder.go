package triangulate

import (
	"github.com/example/delaunay/internal/geom"
)

// edge is an undirected pair of vertex indices.
type edge struct{ u, v int }

// Builder incrementally builds a Delaunay mesh. Construct it with NewBuilder
// on a validated point set; InsertPoint may then be called any number of
// times (including for points beyond the original set). Finish returns the
// ghost-free triangulation.
type Builder struct {
	// local holds all coordinates in the recentred frame: local[i] is
	// pts[i]-anchor. Indices never move, so topology is translation
	// independent.
	local []geom.Point
	// anchor is the point that was subtracted; circumcentres are converted
	// back to the caller's frame by adding it.
	anchor geom.Point
	pred   geom.PredCtx
	tris   []Tri
	nReal  int // number of non-super vertices at construction time
}

// NewBuilder prepares the triangulator: recentres the set, allocates the
// super triangle and inserts every given point in order.
func NewBuilder(pts []geom.Point) *Builder {
	// Recentre on the first point. Subtraction of same-sign nearby floats is
	// exact (Sterbenz) when the input is far translated, and in general the
	// residual error of p-q in float64 is bounded by 2*eps*|p|.
	anchor := pts[0]
	radius := geom.TranslationNoiseRadius(anchor)
	local := make([]geom.Point, len(pts))
	for i, p := range pts {
		local[i] = p.Sub(anchor)
	}

	b := &Builder{
		local:  local,
		anchor: anchor,
		pred:   geom.PredCtx{Radius: radius},
	}

	// Place the super triangle far enough that ghost circumcircles cannot
	// bulge inward past a hull edge; the bound uses the closest approach of
	// any input point to the convex hull.
	bounds := geom.BoundsOf(local)
	hull := geom.ConvexHull(local)
	clearance := geom.HullClearance(local, hull)
	sp, seed, nReal := superTriangle(local, bounds, clearance)
	b.local = append(b.local, sp...)
	b.nReal = nReal
	b.tris = []Tri{seed}

	for i := 0; i < nReal; i++ {
		b.insert(i)
	}
	return b
}

// InsertPoint adds a brand-new user point (in the original coordinate frame)
// into an existing triangulation and returns its vertex index. Before
// insertion the super vertices occupy the final three slots nReal..nReal+2;
// afterwards they occupy nReal+1..nReal+3, so only references strictly at or
// above the old super-base are shifted — real vertex indices are untouched.
func (b *Builder) InsertPoint(p geom.Point) int {
	idx := b.nReal
	oldGhostBase := b.nReal
	for i := range b.tris {
		t := &b.tris[i]
		if t.A >= oldGhostBase {
			t.A++
		}
		if t.B >= oldGhostBase {
			t.B++
		}
		if t.C >= oldGhostBase {
			t.C++
		}
	}
	b.local = insertBeforeSuper(b.local, p.Sub(b.anchor))
	b.nReal++
	b.insert(idx)
	return idx
}

// insertBeforeSuper inserts q in the slot just before the final three super
// vertices and shifts the super vertices one slot to the right.
func insertBeforeSuper(pts []geom.Point, q geom.Point) []geom.Point {
	n := len(pts) - 3
	out := make([]geom.Point, 0, len(pts)+1)
	out = append(out, pts[:n]...)
	out = append(out, q)
	out = append(out, pts[n:]...)
	return out
}

// pointInTriangleCCW reports whether p lies in the closed CCW triangle
// (a,b,c): p must be on the left side of (or on) every directed edge.
func pointInTriangleCCW(a, b, c, p geom.Point, pc geom.PredCtx) bool {
	return pc.Orient2d(a, b, p) != geom.Negative &&
		pc.Orient2d(b, c, p) != geom.Negative &&
		pc.Orient2d(c, a, p) != geom.Negative
}

// insert runs one Bowyer-Watson step for vertex pi.
func (b *Builder) insert(pi int) {
	p := b.local[pi]

	// 1. Find the bad triangles.
	//    - strictly inside the circumdisk => bad.
	//    - For REAL triangles the tolerance ON band (cocircular) is bad only
	//      when p is in the triangle's closed interior: that tie-break keeps
	//      the pre-existing diagonal for exactly cocircular quadruples and
	//      prevents jitter flip-flopping.
	//    - For GHOST triangles (touching a super vertex) the STRICT predicate
	//      (round-off interval only, no geometric fuzzy band) is used. The
	//      ghost vertices are extremely far away, so a relative fuzzy band
	//      would swallow every sign and corrupt the cavity near hull edges.
	ghostBase := len(b.local) - 3
	isGhost := func(t Tri) bool {
		return t.A >= ghostBase || t.B >= ghostBase || t.C >= ghostBase
	}
	bad := make(map[int]struct{}, 8)
	for ti, t := range b.tris {
		a, c, d := t.verts(b.local)
		s := b.pred.InCircle(a, c, d, p)
		if isGhost(t) {
			s = b.pred.InCircleStrict(a, c, d, p)
		}
		switch s {
		case geom.Positive:
			bad[ti] = struct{}{}
		case geom.On:
			if pointInTriangleCCW(a, c, d, p, b.pred) {
				bad[ti] = struct{}{}
			}
		}
	}

	// 2. Collect each cavity edge and its bad incident triangles.
	type edgeInfo struct{ badTris []int }
	infos := make(map[edge]*edgeInfo)
	order := make([]edge, 0, 3*len(bad))
	for ti := range bad {
		t := b.tris[ti]
		for _, k0 := range [3]edge{{t.A, t.B}, {t.B, t.C}, {t.C, t.A}} {
			k := k0
			if k.u > k.v {
				k.u, k.v = k.v, k.u
			}
			info, ok := infos[k]
			if !ok {
				info = &edgeInfo{}
				infos[k] = info
				order = append(order, k)
			}
			info.badTris = append(info.badTris, ti)
		}
	}

	// 3. Re-stitch: fan one CCW triangle (pi,u,v) onto each cavity boundary
	//    edge (an edge with exactly one bad owner). The edge is oriented
	//    explicitly so pi lies on its left; this is independent of which bad
	//    triangle carried it and of map iteration order.
	newTris := make([]Tri, 0, len(order))
	for _, k := range order {
		if len(infos[k].badTris) != 1 {
			continue // interior edge shared by two bad triangles
		}
		u, v := k.u, k.v
		switch b.pred.Orient2d(b.local[u], b.local[v], p) {
		case geom.Negative:
			u, v = v, u
		case geom.On:
			continue // pi on the edge: zero-area fan triangle
		}
		newTris = append(newTris, Tri{pi, u, v})
	}

	// 4. Keep surviving triangles in their old slots (skip bad), append new.
	kept := b.tris[:0]
	for ti, t := range b.tris {
		if _, isBad := bad[ti]; !isBad {
			kept = append(kept, t)
		}
	}
	b.tris = append(kept, newTris...)
}

// Finish strips every triangle touching a super vertex and returns the
// ghost-free triangles together with the local-coordinate table and anchor
// (so Voronoi centres can be converted back). The result order is canonical
// (lexicographic on the index triple) so responses are deterministic.
func (b *Builder) Finish() ([]Tri, []geom.Point, geom.Point) {
	n0 := b.nReal // super vertices occupy nReal..nReal+2
	s0, s1, s2 := n0, n0+1, n0+2
	out := make([]Tri, 0, len(b.tris))
	for _, t := range b.tris {
		if t.A == s0 || t.A == s1 || t.A == s2 ||
			t.B == s0 || t.B == s1 || t.B == s2 ||
			t.C == s0 || t.C == s1 || t.C == s2 {
			continue
		}
		out = append(out, t)
	}
	// Return only real-point local coordinates.
	realLocal := append([]geom.Point(nil), b.local[:n0]...)
	sortTris(out)
	return out, realLocal, b.anchor
}

// sortTris orders triangles lexicographically.
func sortTris(ts []Tri) {
	for i := 1; i < len(ts); i++ {
		for j := i; j > 0 && lessTri(ts[j], ts[j-1]); j-- {
			ts[j-1], ts[j] = ts[j], ts[j-1]
		}
	}
}

func lessTri(a, b Tri) bool {
	if a.A != b.A {
		return a.A < b.A
	}
	if a.B != b.B {
		return a.B < b.B
	}
	return a.C < b.C
}

// Triangulate is the convenience one-shot entry point. It returns the
// ghost-free triangles and the validated points (unchanged order).
func Triangulate(pts []geom.Point) ([]Tri, error) {
	b := NewBuilder(pts)
	tris, _, _ := b.Finish()
	return tris, nil
}
