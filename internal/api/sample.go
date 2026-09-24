package api

import (
	"math"

	"github.com/example/delaunay/internal/geom"
)

// SampleGridSize is the side count of the built-in regular grid.
const SampleGridSize = 5 // 5x5 points

// SampleJitter bounds the interior perturbation applied to the grid.
const SampleJitter = 1e-3

// SamplePoints returns the built-in review case: a 5x5 regular grid with a
// small deterministic perturbation. Interior points receive a bounded jitter;
// boundary points are pushed slightly INWARD so all 4 corners remain the
// extreme points, keeping the convex hull a 4-gon regardless of jitter sign.
//
// With N=25 points and a 4-corner hull (H=4) the planar Euler relation
// V - E + F = 2 gives T = 2N - H - 2 = 44 triangles, which the test suite
// pins exactly.
func SamplePoints() []geom.Point {
	n := SampleGridSize
	pts := make([]geom.Point, 0, n*n)
	j := SampleJitter
	for r := 0; r < n; r++ {
		for c := 0; c < n; c++ {
			x := float64(c)
			y := float64(r)

			interior := r > 0 && r < n-1 && c > 0 && c < n-1
			switch {
			case interior:
				x += j * deterministicPerturb(r*n+c, 0)
				y += j * deterministicPerturb(r*n+c, 1)
			case (r == 0 || r == n-1) && (c == 0 || c == n-1):
				// corners stay exactly on the outer square
			case r == 0: // bottom edge: push up
				y += j * math.Abs(deterministicPerturb(r*n+c, 3))
			case r == n-1: // top edge: push down
				y -= j * math.Abs(deterministicPerturb(r*n+c, 3))
			case c == 0: // left edge: push right
				x += j * math.Abs(deterministicPerturb(r*n+c, 4))
			case c == n-1: // right edge: push left
				x -= j * math.Abs(deterministicPerturb(r*n+c, 4))
			}
			pts = append(pts, geom.Point{X: x, Y: y})
		}
	}
	return pts
}

// deterministicPerturb is a small, reproducible pseudo-random value in [-1,1]
// derived from two integer seeds (no global RNG state, stable across runs).
func deterministicPerturb(i, salt int) float64 {
	h := uint32(i*374761393 + salt*668265263)
	h = (h ^ (h >> 13)) * 1274126177
	h ^= h >> 16
	// Map [0,1) to [-1,1).
	v := float64(h%10000) / 10000.0
	return 2*v - 1
}

// ExpectedSampleTriangles is the Euler-predicted count for the sample.
const ExpectedSampleTriangles = 2*(SampleGridSize*SampleGridSize) - 4 - 2
