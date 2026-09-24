// Package voronoi derives the Voronoi diagram dual to a Delaunay
// triangulation: every Delaunay triangle contributes one Voronoi vertex at
// its exact circumcentre, and every interior Delaunay edge shared by two
// triangles contributes one Voronoi edge joining their circumcentres.
// Delaunay edges on the convex hull dualise to unbounded Voronoi rays.
package voronoi

import (
	"math"
	"sort"

	"github.com/example/delaunay/internal/geom"
	"github.com/example/delaunay/internal/triangulate"
)

// Vertex is one Voronoi vertex (the circumcentre of a Delaunay triangle),
// given in the original point-set coordinate frame.
type Vertex struct {
	geom.Point
	Triangle [3]int `json:"triangle"` // Delaunay triangle producing it
}

// Edge joins two Voronoi vertices by their indices into Vertices.
type Edge struct {
	A int `json:"a"`
	B int `json:"b"`
}

// Ray is an unbounded Voronoi edge starting at vertex A and extending along
// the unit direction (DX, DY). It is the dual of a convex-hull Delaunay edge.
type Ray struct {
	A          int     `json:"a"`
	DX         float64 `json:"dx"`
	DY         float64 `json:"dy"`
	DelaunayUV [2]int  `json:"delaunay_uv"`
}

// Diagram is the full dual structure.
type Diagram struct {
	Vertices []Vertex `json:"vertices"`
	Edges    []Edge   `json:"edges"`
	Rays     []Ray    `json:"rays"`
}

type edgeKey struct{ u, v int }

func canon(u, v int) edgeKey {
	if u > v {
		u, v = v, u
	}
	return edgeKey{u, v}
}

// Build computes the Voronoi diagram. pts are the input points in their
// original frame and tris the ghost-free Delaunay triangulation. Voronoi
// vertices follow the order of tris, so vertex i is the circumcentre of
// tris[i].
func Build(pts []geom.Point, tris []triangulate.Tri) *Diagram {
	d := &Diagram{Vertices: make([]Vertex, 0, len(tris))}

	// For every undirected Delaunay edge, record the (at most two) incident
	// triangle/vertex indices.
	incident := make(map[edgeKey][]int)

	for ti, t := range tris {
		a, b, c := pts[t.A], pts[t.B], pts[t.C]
		cc := geom.Circumcenter(a, b, c)
		d.Vertices = append(d.Vertices, Vertex{
			Point:    cc,
			Triangle: [3]int{t.A, t.B, t.C},
		})
		for _, k := range [3]edgeKey{
			canon(t.A, t.B), canon(t.B, t.C), canon(t.C, t.A),
		} {
			incident[k] = append(incident[k], ti)
		}
	}

	// Deterministic key order.
	keys := make([]edgeKey, 0, len(incident))
	for k := range incident {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool {
		if keys[i].u != keys[j].u {
			return keys[i].u < keys[j].u
		}
		return keys[i].v < keys[j].v
	})

	for _, k := range keys {
		slots := incident[k]
		if len(slots) == 2 {
			d.Edges = append(d.Edges, Edge{A: slots[0], B: slots[1]})
			continue
		}
		// Hull edge: one incident triangle; emit a ray pointing away from
		// the triangle's third vertex.
		t := tris[slots[0]]
		w, _ := thirdVertex(t, k.u, k.v)
		d.Rays = append(d.Rays, makeRay(pts, k, slots[0], w))
	}
	return d
}

// thirdVertex returns t's vertex that is neither u nor v.
func thirdVertex(t triangulate.Tri, u, v int) (int, bool) {
	has := func(x int) bool { return x == u || x == v }
	switch {
	case has(t.A) && has(t.B):
		return t.C, true
	case has(t.B) && has(t.C):
		return t.A, true
	case has(t.C) && has(t.A):
		return t.B, true
	default:
		return 0, false
	}
}

// makeRay builds the unbounded ray dual to hull edge k=(u,v), starting at the
// circumcentre vertex `id`, pointing outward (away from the adjacent
// triangle's third vertex w).
func makeRay(pts []geom.Point, k edgeKey, id, w int) Ray {
	u, v, third := pts[k.u], pts[k.v], pts[w]
	ex, ey := v.X-u.X, v.Y-u.Y
	lenE := math.Hypot(ex, ey)
	ex, ey = ex/lenE, ey/lenE

	// Unit normal to the hull edge.
	nx, ny := -ey, ex
	mid := geom.Point{X: (u.X + v.X) / 2, Y: (u.Y + v.Y) / 2}
	// Side of the third vertex relative to directed edge u->v...
	sideThird := (v.X-u.X)*(third.Y-u.Y) - (v.Y-u.Y)*(third.X-u.X)
	// ...and side of a probe one normal-length past the edge midpoint.
	probe := geom.Point{X: mid.X + nx, Y: mid.Y + ny}
	sideProbe := (v.X-u.X)*(probe.Y-u.Y) - (v.Y-u.Y)*(probe.X-u.X)
	if sideThird*sideProbe > 0 {
		nx, ny = -nx, -ny // probe is on the interior side: flip outward
	}
	return Ray{A: id, DX: nx, DY: ny, DelaunayUV: [2]int{k.u, k.v}}
}
