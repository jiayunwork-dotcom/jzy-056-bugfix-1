package triangulate

import (
	"math"
	"math/rand"
	"testing"

	"github.com/example/delaunay/internal/geom"
)

// squareSet returns the four corners of an axis-aligned square.
func squareSet() []geom.Point {
	return []geom.Point{{X: 0, Y: 0}, {X: 1, Y: 0}, {X: 1, Y: 1}, {X: 0, Y: 1}}
}

// randomSet builds a reproducible pseudo-random point cloud.
func randomSet(n int, seed int64) []geom.Point {
	r := rand.New(rand.NewSource(seed))
	pts := make([]geom.Point, n)
	used := make(map[geom.Point]struct{}, n)
	for i := range pts {
		for {
			p := geom.Point{X: r.Float64() * 100, Y: r.Float64() * 100}
			if _, ok := used[p]; !ok {
				used[p] = struct{}{}
				pts[i] = p
				break
			}
		}
	}
	return pts
}

// assertValidMesh runs every geometric invariant on the finished mesh.
func assertValidMesh(t *testing.T, pts []geom.Point, tris []Tri) {
	t.Helper()
	if len(tris) == 0 {
		t.Fatal("mesh contains no triangles")
	}
	rep := Verify(pts, tris)
	if !rep.AllCCW {
		t.Error("not every triangle is counter-clockwise")
	}
	if !rep.EmptyCircle {
		t.Error("empty circumcircle property violated")
	}
	if !rep.NonOverlapping {
		t.Error("triangles contain crossing/intersecting edges")
	}

	// Compare areas in a recentred local frame so that large rigid
	// translations do not inflate float64 round-off in the area sums.
	anchor := pts[0]
	local := make([]geom.Point, len(pts))
	for i, p := range pts {
		local[i] = p.Sub(anchor)
	}
	hull := geom.ConvexHull(local)
	hullArea := geom.PolygonArea(hull)
	tol := 1e-9 * math.Max(1, hullArea)
	if math.Abs(rep.AreaSum-hullArea) > tol {
		t.Errorf("area conservation broken: triangle area sum %.12f != hull area %.12f (delta %.3e)",
			rep.AreaSum, hullArea, rep.AreaSum-hullArea)
	}

	// Euler: T = 2N - H - 2 where H is the number of hull edges.
	h := len(geom.ConvexHull(pts))
	wantT := 2*len(pts) - h - 2
	if len(tris) != wantT {
		t.Errorf("Euler count broken: got %d triangles, want %d (N=%d,H=%d)",
			len(tris), wantT, len(pts), h)
	}

	// No vertex index may reference a super vertex.
	for _, tr := range tris {
		for _, v := range [3]int{tr.A, tr.B, tr.C} {
			if v < 0 || v >= len(pts) {
				t.Errorf("triangle %v references out-of-range/ghost vertex %d", tr, v)
			}
		}
	}
}

func TestSquare(t *testing.T) {
	pts, err := geom.ParsePoints(squareSet())
	if err != nil {
		t.Fatal(err)
	}
	tris, err := Triangulate(pts)
	if err != nil {
		t.Fatal(err)
	}
	if len(tris) != 2 {
		t.Fatalf("square triangulated into %d triangles, want 2", len(tris))
	}
	assertValidMesh(t, pts, tris)
}

func TestRandomClouds_EmptyCircleAndArea(t *testing.T) {
	for _, seed := range []int64{1, 2, 3, 42, 100} {
		pts, err := geom.ParsePoints(randomSet(60, seed))
		if err != nil {
			t.Fatalf("seed %d: %v", seed, err)
		}
		tris, err := Triangulate(pts)
		if err != nil {
			t.Fatalf("seed %d: %v", seed, err)
		}
		t.Run("seed", func(t *testing.T) { assertValidMesh(t, pts, tris) })
	}
}

func TestTranslationInvariance(t *testing.T) {
	base, err := geom.ParsePoints(randomSet(40, 7))
	if err != nil {
		t.Fatal(err)
	}
	tris0, err := Triangulate(base)
	if err != nil {
		t.Fatal(err)
	}
	baseSet := triSetOf(tris0)

	shifts := []geom.Point{
		{X: 1e6, Y: 0}, {X: 0, Y: 1e6}, {X: 1e9, Y: -1e9},
		{X: -123456789.5, Y: 987654321.25}, {X: 3.14159, Y: -2.71828},
	}
	for _, sh := range shifts {
		shifted := make([]geom.Point, len(base))
		for i, p := range base {
			shifted[i] = p.Add(sh)
		}
		sp, err := geom.ParsePoints(shifted)
		if err != nil {
			t.Fatalf("shift %v: %v", sh, err)
		}
		trisS, err := Triangulate(sp)
		if err != nil {
			t.Fatalf("shift %v: %v", sh, err)
		}
		if len(trisS) != len(tris0) {
			t.Fatalf("shift %v: triangle count changed %d -> %d", sh, len(tris0), len(trisS))
		}
		if got := triSetOf(trisS); !got.equal(baseSet) {
			t.Fatalf("shift %v changed topology", sh)
		}
		// Geometry invariants must also hold in the translated frame.
		assertValidMesh(t, sp, trisS)
	}
}

func TestInsertInteriorPoint_AddsTwoTriangles(t *testing.T) {
	pts, err := geom.ParsePoints(randomSet(30, 11))
	if err != nil {
		t.Fatal(err)
	}
	b := NewBuilder(pts)
	tris0, _, _ := b.Finish()

	// Pick a point safely inside the hull: average of three hull-adjacent
	// points (convex combination of input points is inside the hull).
	hull := geom.ConvexHull(pts)
	inner := hull[0].Mul(0.34).Add(hull[1].Mul(0.33)).Add(hull[2].Mul(0.33))
	// Guard against coinciding with an existing point.
	inner = inner.Add(geom.Point{X: 1e-6, Y: 2e-6})

	b.InsertPoint(inner)
	tris1, _, _ := b.Finish()
	if got := len(tris1) - len(tris0); got != 2 {
		t.Fatalf("inserting an interior point changed triangle count by %d, want +2", got)
	}
	all := make([]geom.Point, 0, len(pts)+1)
	all = append(all, pts...)
	all = append(all, inner)
	assertValidMesh(t, all, tris1)
}

// triSet is an unordered set of triangles for topology comparison.
type triSet map[[3]int]struct{}

func triSetOf(tris []Tri) triSet {
	s := make(triSet, len(tris))
	for _, t := range tris {
		s[[3]int{t.A, t.B, t.C}] = struct{}{}
	}
	return s
}

func (s triSet) equal(o triSet) bool {
	if len(s) != len(o) {
		return false
	}
	for k := range s {
		if _, ok := o[k]; !ok {
			return false
		}
	}
	return true
}
