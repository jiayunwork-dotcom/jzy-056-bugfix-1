// Package geom holds the geometry primitives, validated point-set parsing,
// numeric predicates with tolerances and convex-hull helpers used by the
// Delaunay/Voronoi service.
package geom

import (
	"errors"
	"math"
)

// Point is a 2D point. Index is assigned by the parser; geometry-only code
// treats points as un-indexed coordinates.
type Point struct {
	X float64 `json:"x"`
	Y float64 `json:"y"`
}

// Add returns p+q.
func (p Point) Add(q Point) Point { return Point{p.X + q.X, p.Y + q.Y} }

// Sub returns p-q.
func (p Point) Sub(q Point) Point { return Point{p.X - q.X, p.Y - q.Y} }

// Mul returns p scaled by s.
func (p Point) Mul(s float64) Point { return Point{p.X * s, p.Y * s} }

// Input validation error codes. The API layer maps these onto HTTP statuses.
var (
	ErrTooFewPoints = errors.New("AT_LEAST_THREE_POINTS: at least three points are required")
	ErrDuplicate    = errors.New("DUPLICATE_POINT: the point set contains duplicate coordinates")
	ErrInvalidCoord = errors.New("INVALID_COORDINATE: coordinates must be finite (NaN and +/-Inf are rejected)")
	ErrCollinear    = errors.New("DEGENERATE_POINT_SET: all points are collinear and enclose no area")
)

// ParsePoints validates a raw point set. Checks are ordered so that the most
// basic problem (non-finite coordinate) is reported first.
//
// Rules:
//   - every coordinate must be finite (no NaN, no +/-Inf);
//   - at least three points are required;
//   - duplicate (unordered) coordinates are rejected;
//   - strictly collinear sets are rejected (they enclose no area).
func ParsePoints(raw []Point) ([]Point, error) {
	for _, p := range raw {
		if math.IsNaN(p.X) || math.IsNaN(p.Y) || math.IsInf(p.X, 0) || math.IsInf(p.Y, 0) {
			return nil, ErrInvalidCoord
		}
	}
	if len(raw) < 3 {
		return nil, ErrTooFewPoints
	}

	pts := make([]Point, len(raw))
	copy(pts, raw)
	if err := checkDuplicates(pts); err != nil {
		return nil, err
	}
	if IsCollinear(pts) {
		return nil, ErrCollinear
	}
	return pts, nil
}

// checkDuplicates reports whether any two points share exactly the same
// coordinates (the input contract says duplicates must be declared).
func checkDuplicates(pts []Point) error {
	seen := make(map[Point]struct{}, len(pts))
	for _, p := range pts {
		if _, ok := seen[p]; ok {
			return ErrDuplicate
		}
		seen[p] = struct{}{}
	}
	return nil
}

// Bounds returns the min/max corner and the diagonal span of the point set.
type Bounds struct {
	Min Point
	Max Point
}

// Span returns the width, height and diameter (diagonal) of the bounds.
func (b Bounds) Span() (w, h, d float64) {
	w = b.Max.X - b.Min.X
	h = b.Max.Y - b.Min.Y
	d = math.Hypot(w, h)
	return
}

// BoundsOf computes the axis aligned bounding box.
func BoundsOf(pts []Point) Bounds {
	b := Bounds{Min: Point{math.Inf(1), math.Inf(1)}, Max: Point{math.Inf(-1), math.Inf(-1)}}
	for _, p := range pts {
		if p.X < b.Min.X {
			b.Min.X = p.X
		}
		if p.Y < b.Min.Y {
			b.Min.Y = p.Y
		}
		if p.X > b.Max.X {
			b.Max.X = p.X
		}
		if p.Y > b.Max.Y {
			b.Max.Y = p.Y
		}
	}
	return b
}
