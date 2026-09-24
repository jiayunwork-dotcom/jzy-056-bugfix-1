package triangulate

import (
	"math"
	"math/rand"
	"testing"

	"github.com/example/delaunay/internal/geom"
)

func TestStressManySeeds(t *testing.T) {
	for seed := int64(0); seed < 200; seed++ {
		n := 5 + int(seed%60)
		r := rand.New(rand.NewSource(seed))
		pts := make([]geom.Point, 0, n)
		used := map[geom.Point]struct{}{}
		for len(pts) < n {
			p := geom.Point{X: r.Float64() * 200, Y: r.Float64() * 200}
			if _, ok := used[p]; !ok {
				used[p] = struct{}{}
				pts = append(pts, p)
			}
		}
		pp, err := geom.ParsePoints(pts)
		if err != nil {
			continue // reject collinear draws
		}
		tris, err := Triangulate(pp)
		if err != nil {
			t.Fatalf("seed %d: %v", seed, err)
		}
		rep := Verify(pp, tris)
		hull := geom.ConvexHull(pp)
		var anchor geom.Point = pp[0]
		local := make([]geom.Point, len(pp))
		for i := range pp {
			local[i] = pp[i].Sub(anchor)
		}
		ha := geom.PolygonArea(geom.ConvexHull(local))
		if !rep.EmptyCircle || !rep.AllCCW || !rep.NonOverlapping {
			t.Fatalf("seed %d: invalid mesh (ec=%v ccw=%v cross=%v)", seed, rep.EmptyCircle, rep.AllCCW, rep.NonOverlapping)
		}
		if math.Abs(rep.AreaSum-ha) > 1e-7*math.Max(1, ha) {
			t.Fatalf("seed %d: area %.6f != hull %.6f", seed, rep.AreaSum, ha)
		}
		if len(tris) != 2*len(pp)-len(hull)-2 {
			t.Fatalf("seed %d: T=%d want %d", seed, len(tris), 2*len(pp)-len(hull)-2)
		}
	}
}

func TestStressAnisotropic(t *testing.T) {
	// Elongated / thin but numerically well-conditioned sets (ratios that
	// genuinely arise in terrain / mesh preprocessing, well below the
	// ~1e6:1 regime where double precision cannot separate near-cocircular
	// diagonals at all).
	cases := [][]geom.Point{
		{{X: 0, Y: 0}, {X: 1000, Y: 0.01}, {X: 1000, Y: 1}, {X: 1, Y: 1}},
		{{X: 0, Y: 0}, {X: 100, Y: 0}, {X: 100, Y: 0.1}, {X: 0, Y: 0.1}, {X: 50, Y: 0.05}},
		{{X: 0, Y: 0}, {X: 50, Y: 1}, {X: 100, Y: 0}, {X: 100, Y: 20}, {X: 0, Y: 20}},
	}
	for i, pts := range cases {
		pp, err := geom.ParsePoints(pts)
		if err != nil {
			t.Fatalf("case %d: %v", i, err)
		}
		tris, _ := Triangulate(pp)
		rep := Verify(pp, tris)
		hull := geom.ConvexHull(pp)
		if !rep.EmptyCircle || !rep.AllCCW || !rep.NonOverlapping {
			t.Fatalf("case %d: invalid mesh (ec=%v ccw=%v cross=%v)",
				i, rep.EmptyCircle, rep.AllCCW, rep.NonOverlapping)
		}
		if len(tris) != 2*len(pp)-len(hull)-2 {
			t.Fatalf("case %d: T=%d want %d", i, len(tris), 2*len(pp)-len(hull)-2)
		}
	}
}
