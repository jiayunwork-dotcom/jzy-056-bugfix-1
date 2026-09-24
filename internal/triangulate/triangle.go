// Package triangulate implements an incremental Bowyer-Watson Delaunay
// triangulator. All public triangles are returned with their vertices in
// counter-clockwise order; the same winding convention is reused while
// stitching the re-triangulated cavity, which keeps boundary edges
// consistently oriented.
package triangulate

import (
	"github.com/example/delaunay/internal/geom"
)

// Tri is a triangle as three indices into the triangulation's coordinate
// table. Order is guaranteed counter-clockwise.
type Tri struct {
	A, B, C int
}

// verts returns the three points of t.
func (t Tri) verts(pts []geom.Point) (geom.Point, geom.Point, geom.Point) {
	return pts[t.A], pts[t.B], pts[t.C]
}

// superTriangle builds the virtual enclosing triangle.
//
// A fixed, merely "large" super triangle is not enough: the circumcircle of a
// ghost fan triangle bulges a finite distance INWARD across its hull edge, so
// a point inserted just inside that hull edge can fall inside the ghost's
// circumdisk and corrupt the cavity. For a hull edge of length ~D with a
// neighbouring point at perpendicular clearance delta, the super vertices
// must sit at roughly D^2/(2 delta) beyond the edge. The recentred coordinate
// frame plus the interval-error predicates keep the large absolute
// coordinates numerically sound.
//
// superCircumradius bounds how far the virtual super vertices must be so a
// ghost circumcircle cannot bulge inward past a hull edge whose nearest real
// point sits at perpendicular clearance delta.
//
// For an equilateral super triangle, a ghost circle is the circumcircle of a
// fan triangle spanning a hull edge (length <= D) and one super vertex. In the
// limit where the super vertex is far, that circle behaves near the hull edge
// like the line itself with a quadratic bulge; placing the vertex at roughly
// D^2/(2 delta) beyond the edge keeps the bulge below delta. A factor of 4 is
// used for margin and to cover the non-limiting geometry.
func superCircumradius(diag, clearance float64) float64 {
	if clearance <= 0 {
		clearance = 1
	}
	r := 4 * (diag*diag/(2*clearance) + diag)
	return r
}

func superTriangle(real []geom.Point, b geom.Bounds, clearance float64) (sp []geom.Point, t Tri, nReal int) {
	nReal = len(real)
	_, _, diag := b.Span()
	if diag == 0 {
		diag = 1
	}
	r := superCircumradius(diag, clearance)

	cx := (b.Min.X + b.Max.X) / 2
	cy := (b.Min.Y + b.Max.Y) / 2

	// Equilateral triangle on circumcircle radius r, one vertex straight up.
	p0 := geom.Point{X: cx, Y: cy + r}
	p1 := geom.Point{X: cx - r*cos30, Y: cy - r*0.5}
	p2 := geom.Point{X: cx + r*cos30, Y: cy - r*0.5}
	// p0 -> p1 -> p2 is counter-clockwise in x-right/y-up.
	return []geom.Point{p0, p1, p2}, Tri{nReal, nReal + 1, nReal + 2}, nReal
}

const cos30 = 0.86602540378443864676372317075293618 // sqrt(3)/2
