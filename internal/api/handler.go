// Package api is the HTTP transport layer. It owns request binding, JSON
// serialisation and status codes; all geometry lives in the geom,
// triangulate and voronoi packages.
package api

import (
	"bytes"
	"encoding/json"
	"errors"
	"math"
	"net/http"
	"regexp"
	"strconv"

	"github.com/gin-gonic/gin"

	"github.com/example/delaunay/internal/geom"
	"github.com/example/delaunay/internal/triangulate"
	"github.com/example/delaunay/internal/voronoi"
)

var errBadJSON = errors.New("INVALID_JSON")

// nonFiniteToken matches the JavaScript-style non-finite number literals that
// strict JSON rejects (NaN, Infinity, -Infinity). They are rewritten to a
// finite-looking number that overflows float64, so the geometry layer can
// classify them as INVALID_COORDINATE instead of a generic parse failure.
var nonFiniteToken = regexp.MustCompile(`[-+]?(?:NaN|Infinity)`)

// decodePoints parses a request body into validated geometry points.
func decodePoints(body []byte) ([]geom.Point, error) {
	cleaned := nonFiniteToken.ReplaceAllFunc(body, func(m []byte) []byte {
		if bytes.HasPrefix(m, []byte("-")) {
			return []byte("-1e9999") // parses to -Inf in float64
		}
		return []byte("1e9999") // parses to +Inf (NaN has no ordering either)
	})

	var req struct {
		Points []struct {
			X any `json:"x"`
			Y any `json:"y"`
		} `json:"points"`
	}
	dec := json.NewDecoder(bytes.NewReader(cleaned))
	dec.UseNumber()
	if err := dec.Decode(&req); err != nil {
		return nil, errBadJSON
	}
	if len(req.Points) == 0 {
		return nil, errBadJSON
	}

	out := make([]geom.Point, len(req.Points))
	for i, rp := range req.Points {
		x, errX := asCoord(rp.X)
		y, errY := asCoord(rp.Y)
		if errors.Is(errX, geom.ErrInvalidCoord) || errors.Is(errY, geom.ErrInvalidCoord) {
			return nil, geom.ErrInvalidCoord
		}
		if errX != nil || errY != nil {
			return nil, errBadJSON
		}
		out[i] = geom.Point{X: x, Y: y}
	}
	return geom.ParsePoints(out)
}

// asCoord converts one decoded JSON scalar to a finite float64.
func asCoord(v any) (float64, error) {
	n, ok := v.(json.Number)
	if !ok {
		return 0, errBadJSON // missing or non-numeric coordinate field
	}
	f, err := strconv.ParseFloat(string(n), 64)
	if math.IsNaN(f) || math.IsInf(f, 0) {
		return 0, geom.ErrInvalidCoord
	}
	if err != nil {
		return 0, errBadJSON
	}
	return f, nil
}

// PointDTO is the wire representation of an input point.
type PointDTO struct {
	X float64 `json:"x"`
	Y float64 `json:"y"`
}

// TriangleDTO is an index triple into the request's point array.
type TriangleDTO struct {
	A int `json:"a"`
	B int `json:"b"`
	C int `json:"c"`
}

// TriangulateResponse is returned by POST /delaunay.
type TriangulateResponse struct {
	Count     int           `json:"count"`
	Triangles []TriangleDTO `json:"triangles"`
	Delaunay  bool          `json:"delaunay"`
	HullEdges int           `json:"hull_edges"`
	HullArea  float64       `json:"hull_area"`
	AreaSum   float64       `json:"triangle_area_sum"`
	AreaDelta float64       `json:"area_delta"`
}

// VoronoiVertexDTO is a circumcentre plus its source triangle.
type VoronoiVertexDTO struct {
	X        float64 `json:"x"`
	Y        float64 `json:"y"`
	Triangle [3]int  `json:"triangle"`
}

// VoronoiEdgeDTO joins two Voronoi vertices.
type VoronoiEdgeDTO struct {
	A int `json:"a"`
	B int `json:"b"`
}

// VoronoiRayDTO is an unbounded Voronoi edge on the hull.
type VoronoiRayDTO struct {
	A          int     `json:"a"`
	DX         float64 `json:"dx"`
	DY         float64 `json:"dy"`
	DelaunayUV [2]int  `json:"delaunay_uv"`
}

// VoronoiResponse is returned by POST /voronoi.
type VoronoiResponse struct {
	Vertices []VoronoiVertexDTO `json:"vertices"`
	Edges    []VoronoiEdgeDTO   `json:"edges"`
	Rays     []VoronoiRayDTO    `json:"rays"`
}

// ErrorResponse is the single error envelope.
type ErrorResponse struct {
	Error ErrorBody `json:"error"`
}

// ErrorBody carries a stable code plus a human-readable message.
type ErrorBody struct {
	Code    string `json:"code"`
	Message string `json:"message"`
}

// Router wires the service.
func Router() *gin.Engine {
	gin.SetMode(gin.ReleaseMode)
	r := gin.New()
	r.Use(gin.Recovery())

	r.GET("/health", func(c *gin.Context) { c.JSON(http.StatusOK, gin.H{"status": "ok"}) })
	r.POST("/delaunay", handleDelaunay)
	r.POST("/voronoi", handleVoronoi)
	r.GET("/sample", handleSample)
	return r
}

func errorStatus(err error) (int, string) {
	switch {
	case errors.Is(err, errBadJSON):
		return http.StatusBadRequest, "INVALID_JSON"
	case errors.Is(err, geom.ErrTooFewPoints):
		return http.StatusBadRequest, "TOO_FEW_POINTS"
	case errors.Is(err, geom.ErrDuplicate):
		return http.StatusBadRequest, "DUPLICATE_POINT"
	case errors.Is(err, geom.ErrInvalidCoord):
		return http.StatusBadRequest, "INVALID_COORDINATE"
	case errors.Is(err, geom.ErrCollinear):
		return http.StatusUnprocessableEntity, "DEGENERATE_POINT_SET"
	default:
		return http.StatusBadRequest, "BAD_REQUEST"
	}
}

func fail(c *gin.Context, err error) {
	status, code := errorStatus(err)
	msg := err.Error()
	if code == "INVALID_JSON" {
		msg = "request body must be JSON with a non-empty 'points' array of finite numbers"
	}
	c.JSON(status, ErrorResponse{Error: ErrorBody{Code: code, Message: msg}})
}

func handleDelaunay(c *gin.Context) {
	body, err := c.GetRawData()
	if err != nil {
		fail(c, errBadJSON)
		return
	}
	pts, err := decodePoints(body)
	if err != nil {
		fail(c, err)
		return
	}

	tris, err := triangulate.Triangulate(pts)
	if err != nil {
		fail(c, err)
		return
	}
	report := triangulate.Verify(pts, tris)
	hull := geom.ConvexHull(pts)
	hullArea := geom.PolygonArea(hull)

	// Area conservation reported in the caller's own coordinate frame: sum
	// the unsigned triangle areas and compare against the hull area.
	var areaSum float64
	for _, t := range tris {
		areaSum += 0.5 * math.Abs(geom.SignedArea2(pts[t.A], pts[t.B], pts[t.C]))
	}

	out := make([]TriangleDTO, 0, len(tris))
	for _, t := range tris {
		out = append(out, TriangleDTO{A: t.A, B: t.B, C: t.C})
	}
	c.JSON(http.StatusOK, TriangulateResponse{
		Count:     len(tris),
		Triangles: out,
		Delaunay:  report.EmptyCircle && report.AllCCW && report.NonOverlapping,
		HullEdges: len(hull),
		HullArea:  hullArea,
		AreaSum:   areaSum,
		AreaDelta: areaSum - hullArea,
	})
}

func handleVoronoi(c *gin.Context) {
	body, err := c.GetRawData()
	if err != nil {
		fail(c, errBadJSON)
		return
	}
	pts, err := decodePoints(body)
	if err != nil {
		fail(c, err)
		return
	}

	tris, err := triangulate.Triangulate(pts)
	if err != nil {
		fail(c, err)
		return
	}
	d := voronoi.Build(pts, tris)

	vs := make([]VoronoiVertexDTO, 0, len(d.Vertices))
	for _, v := range d.Vertices {
		vs = append(vs, VoronoiVertexDTO{X: v.X, Y: v.Y, Triangle: v.Triangle})
	}
	es := make([]VoronoiEdgeDTO, 0, len(d.Edges))
	for _, e := range d.Edges {
		es = append(es, VoronoiEdgeDTO{A: e.A, B: e.B})
	}
	rs := make([]VoronoiRayDTO, 0, len(d.Rays))
	for _, ry := range d.Rays {
		rs = append(rs, VoronoiRayDTO{A: ry.A, DX: ry.DX, DY: ry.DY, DelaunayUV: ry.DelaunayUV})
	}
	c.JSON(http.StatusOK, VoronoiResponse{Vertices: vs, Edges: es, Rays: rs})
}

// handleSample returns the built-in review case alongside its expected
// topology so a client can re-check the service at any time.
func handleSample(c *gin.Context) {
	pts := SamplePoints()
	dtos := make([]PointDTO, 0, len(pts))
	for _, p := range pts {
		dtos = append(dtos, PointDTO{X: p.X, Y: p.Y})
	}
	hull := geom.ConvexHull(pts)
	c.JSON(http.StatusOK, gin.H{
		"points":              dtos,
		"expected_triangles":  ExpectedSampleTriangles,
		"expected_hull_edges": len(hull),
		"euler_formula":       "T = 2N - H - 2",
	})
}
