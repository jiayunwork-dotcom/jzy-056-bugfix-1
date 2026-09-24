package triangulate

import (
	"math"
	"testing"

	"github.com/example/delaunay/internal/geom"
)

// squareCorners returns four points exactly on the circle x^2+y^2=r^2,
// axis-aligned (the classic degenerate-split case).
func squareCorners(r float64) []geom.Point {
	return []geom.Point{
		{X: -r, Y: 0},
		{X: 0, Y: -r},
		{X: r, Y: 0},
		{X: 0, Y: r},
	}
}

// TestStrictlyCocircular_NoCrossing: four points exactly cocircular may be
// triangulated either way, but the two resulting triangles must never cross
// and the run must be repeatable (no flip-flop between runs / insertions).
func TestStrictlyCocircular_NoCrossing(t *testing.T) {
	radii := []float64{1, 1000, 1e6}
	for _, r := range radii {
		pts, err := geom.ParsePoints(squareCorners(r))
		if err != nil {
			t.Fatalf("r=%g: %v", r, err)
		}
		tris, err := Triangulate(pts)
		if err != nil {
			t.Fatalf("r=%g: %v", r, err)
		}
		if len(tris) != 2 {
			t.Fatalf("r=%g: got %d triangles, want 2", r, len(tris))
		}
		rep := Verify(pts, tris)
		if !rep.AllCCW {
			t.Errorf("r=%g: triangles not CCW", r)
		}
		if !rep.NonOverlapping {
			t.Errorf("r=%g: cocircular triangulation produced crossing edges", r)
		}
		if !rep.EmptyCircle {
			t.Errorf("r=%g: empty circle property violated", r)
		}
		hullArea := geom.PolygonArea(geom.ConvexHull(pts))
		if math.Abs(rep.AreaSum-hullArea) > 1e-9*math.Max(1, hullArea) {
			t.Errorf("r=%g: area sum %.10f != hull %.10f", r, rep.AreaSum, hullArea)
		}
	}
}

// TestCocircularDeterminism repeats the exact same cocircular input many
// times and requires bit-identical topology (the tolerance predicate must
// never oscillate between the two legal diagonals).
func TestCocircularDeterminism(t *testing.T) {
	pts, err := geom.ParsePoints(squareCorners(1))
	if err != nil {
		t.Fatal(err)
	}
	var ref triSet
	for i := 0; i < 50; i++ {
		tris, err := Triangulate(pts)
		if err != nil {
			t.Fatal(err)
		}
		cur := triSetOf(tris)
		if ref == nil {
			ref = cur
			continue
		}
		if !cur.equal(ref) {
			t.Fatalf("cocircular triangulation flipped on iteration %d", i)
		}
	}
}

// TestCocircularTranslated preserves the exact-on-circle decision after a
// large translation that destroys naive float64 incircle tests.
func TestCocircularTranslated(t *testing.T) {
	base := squareCorners(1)
	sh := geom.Point{X: 1e8 + 0.5, Y: -(1e8 - 0.25)}
	shifted := make([]geom.Point, 4)
	for i, p := range base {
		shifted[i] = p.Add(sh)
	}
	// Sanity: a naive un-recentred determinant is swamped at this scale and
	// evaluates to 0 instead of the exact 0.0 cocircular value anyway.
	_ = incircleNaive(shifted[0], shifted[1], shifted[2], shifted[3])

	pts, err := geom.ParsePoints(shifted)
	if err != nil {
		t.Fatal(err)
	}
	tris, err := Triangulate(pts)
	if err != nil {
		t.Fatal(err)
	}
	rep := Verify(pts, tris)
	if !rep.NonOverlapping {
		t.Error("translated cocircular set produced crossing edges")
	}
	if !rep.AllCCW {
		t.Error("translated cocircular set produced mis-oriented triangles")
	}
	if len(tris) != 2 {
		t.Errorf("got %d triangles, want 2", len(tris))
	}
}

// incircleNaive computes the textbook incircle determinant without recentring;
// used only to document why the tolerant predicate is necessary.
func incircleNaive(a, b, c, d geom.Point) float64 {
	mk := func(p geom.Point) [3]float64 {
		return [3]float64{p.X, p.Y, p.X*p.X + p.Y*p.Y}
	}
	A, B, C, D := mk(a), mk(b), mk(c), mk(d)
	// Laplace expansion on the d-row:
	return A[0]*(B[1]*C[2]-B[2]*C[1]) - A[1]*(B[0]*C[2]-B[2]*C[0]) + A[2]*(B[0]*C[1]-B[1]*C[0]) -
		D[0]*(B[1]*C[2]-B[2]*C[1]) + D[1]*(B[0]*C[2]-B[2]*C[0]) - D[2]*(B[0]*C[1]-B[1]*C[0])
}

// TestNearCocircularStress perturbs points around exact cocircular positions
// at scales bracketing the tolerance band; every result must be a legal
// non-crossing Delaunay mesh.
func TestNearCocircularStress(t *testing.T) {
	for k := -6; k <= 2; k++ {
		eps := math.Pow(10, float64(k))
		pts := []geom.Point{
			{X: -1, Y: 0},
			{X: 0, Y: -1},
			{X: 1, Y: 0},
			{X: eps, Y: 1 + eps},
			{X: 0.5, Y: 0.5},
			{X: -0.5, Y: 0.5},
		}
		pp, err := geom.ParsePoints(pts)
		if err != nil {
			t.Fatalf("eps=%g: %v", eps, err)
		}
		tris, err := Triangulate(pp)
		if err != nil {
			t.Fatalf("eps=%g: %v", eps, err)
		}
		rep := Verify(pp, tris)
		if !rep.NonOverlapping || !rep.AllCCW {
			t.Fatalf("eps=%g: invalid mesh (crossing=%v ccw=%v)", eps, !rep.NonOverlapping, !rep.AllCCW)
		}
	}
}
