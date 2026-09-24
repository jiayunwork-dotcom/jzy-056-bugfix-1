package voronoi

import (
	"math"
	"math/rand"
	"testing"

	"github.com/example/delaunay/internal/geom"
	"github.com/example/delaunay/internal/triangulate"
)

func randomPts(n int, seed int64) []geom.Point {
	r := rand.New(rand.NewSource(seed))
	pts := make([]geom.Point, n)
	used := map[geom.Point]struct{}{}
	for i := range pts {
		for {
			p := geom.Point{X: r.Float64() * 50, Y: r.Float64() * 50}
			if _, ok := used[p]; !ok {
				used[p] = struct{}{}
				pts[i] = p
				break
			}
		}
	}
	return pts
}

// distance returns the distance between a vertex and the points of its
// triangle; a genuine circumcentre is equidistant from all three.
func circumradiusResidual(pts []geom.Point, v Vertex) float64 {
	t := v.Triangle
	d := func(i int) float64 { return math.Hypot(v.X-pts[i].X, v.Y-pts[i].Y) }
	d0 := d(t[0])
	d1 := d(t[1])
	d2 := d(t[2])
	return math.Max(math.Abs(d0-d1), math.Abs(d0-d2))
}

func TestVerticesAreCircumcenters(t *testing.T) {
	pts, err := geom.ParsePoints(randomPts(40, 5))
	if err != nil {
		t.Fatal(err)
	}
	tris, err := triangulate.Triangulate(pts)
	if err != nil {
		t.Fatal(err)
	}
	d := Build(pts, tris)
	if len(d.Vertices) != len(tris) {
		t.Fatalf("got %d Voronoi vertices, want %d (one per triangle)", len(d.Vertices), len(tris))
	}
	for i, v := range d.Vertices {
		// Must agree with the independent Circumcenter helper.
		want := geom.Circumcenter(pts[v.Triangle[0]], pts[v.Triangle[1]], pts[v.Triangle[2]])
		if math.Abs(v.X-want.X) > 1e-9*math.Max(1, math.Abs(want.X)) ||
			math.Abs(v.Y-want.Y) > 1e-9*math.Max(1, math.Abs(want.Y)) {
			t.Errorf("vertex %d = (%v,%v), want exact circumcentre (%v,%v)", i, v.X, v.Y, want.X, want.Y)
		}
		if res := circumradiusResidual(pts, v); res > 1e-8 {
			t.Errorf("vertex %d is not equidistant from its triangle (residual %.2e)", i, res)
		}
	}
}

func TestDualTopology(t *testing.T) {
	pts, err := geom.ParsePoints(randomPts(35, 9))
	if err != nil {
		t.Fatal(err)
	}
	tris, err := triangulate.Triangulate(pts)
	if err != nil {
		t.Fatal(err)
	}
	d := Build(pts, tris)

	// Count Delaunay edges and hull edges independently.
	allEdges := map[[2]int]struct{}{}
	for _, tr := range tris {
		for _, e := range [3][2]int{{tr.A, tr.B}, {tr.B, tr.C}, {tr.C, tr.A}} {
			if e[0] > e[1] {
				e[0], e[1] = e[1], e[0]
			}
			allEdges[e] = struct{}{}
		}
	}
	incidence := map[[2]int]int{}
	for _, tr := range tris {
		for _, e := range [3][2]int{{tr.A, tr.B}, {tr.B, tr.C}, {tr.C, tr.A}} {
			if e[0] > e[1] {
				e[0], e[1] = e[1], e[0]
			}
			incidence[e]++
		}
	}
	hullEdges := 0
	for _, n := range incidence {
		if n == 1 {
			hullEdges++
		}
	}
	if len(d.Rays) != hullEdges {
		t.Errorf("got %d rays, want %d (one per hull edge)", len(d.Rays), hullEdges)
	}
	if len(d.Edges) != len(allEdges)-hullEdges {
		t.Errorf("got %d Voronoi edges, want %d interior Delaunay edges",
			len(d.Edges), len(allEdges)-hullEdges)
	}

	// Rays must be unit vectors and point outward.
	for _, ry := range d.Rays {
		n := math.Hypot(ry.DX, ry.DY)
		if math.Abs(n-1) > 1e-12 {
			t.Errorf("ray %v is not unit length (got %.12f)", ry, n)
		}
	}
}

// asymmetricPts is a fixed, deliberately lopsided point set: hull edges of
// very uneven length and a centroid pulled well off-centre. On regular,
// symmetric inputs (square, regular polygon) an inward-pointing hull ray is
// easy to miss by eye; this set makes the outward/inward distinction
// unambiguous.
func asymmetricPts() []geom.Point {
	return []geom.Point{
		{X: 0.0, Y: 0.0},
		{X: 3.7, Y: 0.4},
		{X: 9.2, Y: 0.1}, // long stretched bottom edge
		{X: 11.8, Y: 2.6},
		{X: 10.4, Y: 6.3}, // far right outlier
		{X: 6.1, Y: 4.9},
		{X: 2.2, Y: 5.6},
		{X: 0.6, Y: 2.9},
		{X: 4.4, Y: 2.2}, // interior
		{X: 7.3, Y: 2.8}, // interior
		{X: 5.2, Y: 3.6}, // interior
	}
}

// Every unbounded ray dual to a convex-hull Delaunay edge must point AWAY
// from the point set: walking along it moves further from the cloud, never
// back through it. The outward side of a hull edge is the side opposite the
// adjacent triangle's third vertex.
func TestHullRaysPointOutward(t *testing.T) {
	pts, err := geom.ParsePoints(asymmetricPts())
	if err != nil {
		t.Fatal(err)
	}
	tris, err := triangulate.Triangulate(pts)
	if err != nil {
		t.Fatal(err)
	}
	d := Build(pts, tris)

	// One ray per hull edge; the hull corner count equals its edge count.
	if want := len(geom.ConvexHull(pts)); len(d.Rays) != want {
		t.Fatalf("got %d rays, want %d (one per hull edge)", len(d.Rays), want)
	}

	// The centroid is a strictly positive convex combination of the points,
	// so it lies strictly inside the hull; centroid -> edge midpoint is a
	// reliable outward reference direction for every hull edge.
	var cx, cy float64
	for _, p := range pts {
		cx += p.X
		cy += p.Y
	}
	cx /= float64(len(pts))
	cy /= float64(len(pts))

	for _, ry := range d.Rays {
		u := pts[ry.DelaunayUV[0]]
		v := pts[ry.DelaunayUV[1]]
		ex, ey := v.X-u.X, v.Y-u.Y
		lenE := math.Hypot(ex, ey)
		midX, midY := (u.X+v.X)/2, (u.Y+v.Y)/2

		if n := math.Hypot(ry.DX, ry.DY); math.Abs(n-1) > 1e-12 {
			t.Errorf("ray %v is not unit length (got %.12f)", ry, n)
		}
		// The dual of a Delaunay edge is perpendicular to it.
		if dot := ex*ry.DX + ey*ry.DY; math.Abs(dot) > 1e-9*lenE {
			t.Errorf("ray %v is not perpendicular to its hull edge (dot %.2e)", ry, dot)
		}

		// Acceptance check: same sign as the centroid-outward reference.
		if dot := (midX-cx)*ry.DX + (midY-cy)*ry.DY; dot <= 0 {
			t.Errorf("ray %v points inward: dot with centroid->midpoint reference = %v, want > 0",
				ry, dot)
		}

		// Geometric definition: the ray must leave the edge on the side
		// OPPOSITE the adjacent triangle's third vertex.
		w := -1
		for _, idx := range d.Vertices[ry.A].Triangle {
			if idx != ry.DelaunayUV[0] && idx != ry.DelaunayUV[1] {
				w = idx
			}
		}
		if w < 0 {
			t.Fatalf("ray %v: source triangle has no third vertex off the hull edge", ry)
		}
		sideThird := ex*(pts[w].Y-u.Y) - ey*(pts[w].X-u.X)
		sideRay := ex*ry.DY - ey*ry.DX
		if sideThird*sideRay > 0 {
			t.Errorf("ray %v exits toward the same side as the triangle's third vertex %d (inward)",
				ry, w)
		}
	}
}

func TestTranslationMovesVertices(t *testing.T) {
	base, err := geom.ParsePoints(randomPts(30, 13))
	if err != nil {
		t.Fatal(err)
	}
	sh := geom.Point{X: 1e7 + 0.5, Y: -(1e7 - 0.25)}
	shifted := make([]geom.Point, len(base))
	for i, p := range base {
		shifted[i] = p.Add(sh)
	}

	tb, _ := triangulate.Triangulate(base)
	ts, _ := triangulate.Triangulate(shifted)
	db := Build(base, tb)
	ds := Build(shifted, ts)

	if len(db.Vertices) != len(ds.Vertices) {
		t.Fatalf("vertex count changed under translation: %d -> %d",
			len(db.Vertices), len(ds.Vertices))
	}
	for i := range db.Vertices {
		gx := ds.Vertices[i].X - db.Vertices[i].X
		gy := ds.Vertices[i].Y - db.Vertices[i].Y
		tol := 1e-7 * math.Max(1, math.Max(math.Abs(sh.X), math.Abs(sh.Y)))
		if math.Abs(gx-sh.X) > tol || math.Abs(gy-sh.Y) > tol {
			t.Errorf("vertex %d moved by (%v,%v), want (%v,%v)", i, gx, gy, sh.X, sh.Y)
		}
		if db.Vertices[i].Triangle != ds.Vertices[i].Triangle {
			t.Errorf("vertex %d dual triangle changed under translation", i)
		}
	}
	if len(db.Edges) != len(ds.Edges) {
		t.Errorf("Voronoi edge count changed under translation")
	}
}
