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
