package geom

import (
	"errors"
	"math"
	"testing"
)

func TestParsePoints_RejectsBadInput(t *testing.T) {
	cases := []struct {
		name string
		pts  []Point
		want error
	}{
		{"empty", nil, ErrTooFewPoints},
		{"one", []Point{{0, 0}}, ErrTooFewPoints},
		{"two", []Point{{0, 0}, {1, 1}}, ErrTooFewPoints},
		{"nan x", []Point{{0, 0}, {1, 0}, {math.NaN(), 1}}, ErrInvalidCoord},
		{"nan y", []Point{{0, 0}, {1, 0}, {0, math.NaN()}}, ErrInvalidCoord},
		{"+inf", []Point{{0, 0}, {1, 0}, {0, math.Inf(1)}}, ErrInvalidCoord},
		{"-inf", []Point{{0, 0}, {math.Inf(-1), 0}, {0, 1}}, ErrInvalidCoord},
		{
			"duplicate",
			[]Point{{0, 0}, {1, 0}, {0, 1}, {1, 0}},
			ErrDuplicate,
		},
		{
			"duplicate adjacent",
			[]Point{{2.5, 2.5}, {2.5, 2.5}, {3, 3}},
			ErrDuplicate,
		},
		{
			"collinear horizontal",
			[]Point{{0, 0}, {1, 0}, {2, 0}, {3, 0}},
			ErrCollinear,
		},
		{
			"collinear diagonal",
			[]Point{{0, 0}, {1, 1}, {2, 2}, {-3, -3}},
			ErrCollinear,
		},
		{
			"collinear far translated",
			[]Point{{1e12, 1e12 + 1}, {1e12 + 1, 1e12 + 2}, {1e12 + 5, 1e12 + 6}},
			ErrCollinear,
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			_, err := ParsePoints(tc.pts)
			if !errors.Is(err, tc.want) {
				t.Fatalf("ParsePoints() err = %v, want %v", err, tc.want)
			}
		})
	}
}

func TestParsePoints_AcceptsValidSet(t *testing.T) {
	pts := []Point{{0, 0}, {4, 0}, {0, 3}, {4, 3}}
	got, err := ParsePoints(pts)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(got) != 4 {
		t.Fatalf("got %d points, want 4", len(got))
	}
}

func TestConvexHull_Square(t *testing.T) {
	pts := []Point{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}}
	hull := ConvexHull(pts)
	if len(hull) != 4 {
		t.Fatalf("hull has %d corners, want 4: %v", len(hull), hull)
	}
	if got := PolygonArea(hull); math.Abs(got-1) > 1e-12 {
		t.Fatalf("hull area = %v, want 1", got)
	}
}

func TestConvexHull_TriangleWithEdgePoints(t *testing.T) {
	pts := []Point{{0, 0}, {2, 0}, {4, 0}, {0, 3}}
	hull := ConvexHull(pts)
	if len(hull) != 3 {
		t.Fatalf("hull has %d corners, want 3: %v", len(hull), hull)
	}
	if got := PolygonArea(hull); math.Abs(got-6) > 1e-9 {
		t.Fatalf("hull area = %v, want 6", got)
	}
}

func TestOrient2d(t *testing.T) {
	c := PredCtx{}
	a, b := Point{0, 0}, Point{1, 0}
	if c.Orient2d(a, b, Point{0, 1}) != Positive {
		t.Fatal("(0,0),(1,0),(0,1) must be CCW/Positive")
	}
	if c.Orient2d(a, b, Point{0, -1}) != Negative {
		t.Fatal("(0,0),(1,0),(0,-1) must be CW/Negative")
	}
	if c.Orient2d(a, b, Point{0.5, 0}) != On {
		t.Fatal("point on the line must be On")
	}
}

func TestInCircle_UnitSquare(t *testing.T) {
	c := PredCtx{}
	// CCW triangle: three corners of the unit square.
	a := Point{0, 0}
	b := Point{1, 0}
	cc := Point{0, 1}
	// The circle through these three right-triangle vertices has centre
	// (0.5,0.5) radius ~0.707: the fourth square corner (1,1) lies ON it.
	if s := c.InCircle(a, b, cc, Point{1, 1}); s != On {
		t.Fatalf("square corner incircle = %v, want On", s)
	}
	// The geometric centre of the circle lies strictly inside the disk.
	if s := c.InCircle(a, b, cc, Point{0.5, 0.5}); s != Positive {
		t.Fatalf("circle-centre incircle = %v, want Positive", s)
	}
	// A point near the right-angle corner is strictly inside.
	if s := c.InCircle(a, b, cc, Point{0.1, 0.1}); s != Positive {
		t.Fatalf("near-corner incircle = %v, want Positive", s)
	}
	// A far away point must be outside.
	if s := c.InCircle(a, b, cc, Point{100, 100}); s != Negative {
		t.Fatalf("far point incircle = %v, want Negative", s)
	}
}

func TestCircumcenter(t *testing.T) {
	// Any right triangle: circumcentre is the midpoint of the hypotenuse.
	got := Circumcenter(Point{0, 0}, Point{4, 0}, Point{0, 3})
	want := Point{2, 1.5}
	if math.Abs(got.X-want.X) > 1e-12 || math.Abs(got.Y-want.Y) > 1e-12 {
		t.Fatalf("circumcentre = %v, want %v", got, want)
	}
}
