package api

import (
	"bytes"
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"testing"

	"github.com/example/delaunay/internal/geom"
)

func doPost(t *testing.T, r http.Handler, path string, body any) (int, map[string]any) {
	t.Helper()
	var raw []byte
	var err error
	switch v := body.(type) {
	case string:
		raw = []byte(v)
	default:
		raw, err = json.Marshal(body)
		if err != nil {
			t.Fatal(err)
		}
	}
	req := httptest.NewRequest(http.MethodPost, path, bytes.NewReader(raw))
	req.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	r.ServeHTTP(w, req)

	var out map[string]any
	if w.Body.Len() > 0 {
		if err := json.Unmarshal(w.Body.Bytes(), &out); err != nil {
			t.Fatalf("response is not JSON (%d): %s", w.Code, w.Body.String())
		}
	}
	return w.Code, out
}

func bodyFor(pts []geom.Point) any {
	type p struct {
		X float64 `json:"x"`
		Y float64 `json:"y"`
	}
	ps := make([]p, len(pts))
	for i, q := range pts {
		ps[i] = p{q.X, q.Y}
	}
	return map[string]any{"points": ps}
}

func TestSampleEndpoint_EulerCount(t *testing.T) {
	r := Router()
	req := httptest.NewRequest(http.MethodGet, "/sample", nil)
	w := httptest.NewRecorder()
	r.ServeHTTP(w, req)
	if w.Code != http.StatusOK {
		t.Fatalf("status = %d", w.Code)
	}
	var resp struct {
		Points            []PointDTO `json:"points"`
		ExpectedTriangles int        `json:"expected_triangles"`
		ExpectedHullEdges int        `json:"expected_hull_edges"`
	}
	if err := json.Unmarshal(w.Body.Bytes(), &resp); err != nil {
		t.Fatal(err)
	}
	if len(resp.Points) != SampleGridSize*SampleGridSize {
		t.Fatalf("sample has %d points, want %d", len(resp.Points), SampleGridSize*SampleGridSize)
	}
	if resp.ExpectedTriangles != ExpectedSampleTriangles || ExpectedSampleTriangles != 44 {
		t.Fatalf("expected triangle pin = %d, want 44", resp.ExpectedTriangles)
	}
	if resp.ExpectedHullEdges != 4 {
		t.Fatalf("expected hull edges = %d, want 4", resp.ExpectedHullEdges)
	}

	// Feed the sample back into /delaunay: the service must reproduce the
	// Euler-predicted count exactly, satisfy Delaunay, and conserve area.
	status, out := doPost(t, r, "/delaunay", bodyFor(SamplePoints()))
	if status != http.StatusOK {
		t.Fatalf("status = %d body = %v", status, out)
	}
	if int(out["count"].(float64)) != 44 {
		t.Fatalf("sample triangulated to %v triangles, want 44", out["count"])
	}
	if out["delaunay"] != true {
		t.Fatalf("delaunay = %v, want true (report: %v)", out["delaunay"], out)
	}
	areaSum, hullArea := out["triangle_area_sum"].(float64), out["hull_area"].(float64)
	if d := areaSum - hullArea; d > 1e-8 || d < -1e-8 {
		t.Fatalf("area delta = %v, want ~0", d)
	}
}

func TestDelaunayEndpoint_Square(t *testing.T) {
	r := Router()
	status, out := doPost(t, r, "/delaunay", bodyFor([]geom.Point{
		{X: 0, Y: 0}, {X: 1, Y: 0}, {X: 1, Y: 1}, {X: 0, Y: 1},
	}))
	if status != http.StatusOK {
		t.Fatalf("status = %d body = %v", status, out)
	}
	if int(out["count"].(float64)) != 2 {
		t.Fatalf("count = %v, want 2", out["count"])
	}
	if out["delaunay"] != true {
		t.Fatalf("delaunay = %v", out["delaunay"])
	}
	if int(out["hull_edges"].(float64)) != 4 {
		t.Fatalf("hull_edges = %v, want 4", out["hull_edges"])
	}
}

func TestVoronoiEndpoint_Square(t *testing.T) {
	r := Router()
	status, out := doPost(t, r, "/voronoi", bodyFor([]geom.Point{
		{X: 0, Y: 0}, {X: 2, Y: 0}, {X: 2, Y: 2}, {X: 0, Y: 2},
	}))
	if status != http.StatusOK {
		t.Fatalf("status = %d body = %v", status, out)
	}
	verts := out["vertices"].([]any)
	edges := out["edges"].([]any)
	rays := out["rays"].([]any)
	if len(verts) != 2 {
		t.Fatalf("vertices = %d, want 2", len(verts))
	}
	// Two triangles share one interior edge -> one Voronoi edge; the other
	// four Delaunay edges are hull edges -> four rays.
	if len(edges) != 1 {
		t.Fatalf("edges = %d, want 1", len(edges))
	}
	if len(rays) != 4 {
		t.Fatalf("rays = %d, want 4", len(rays))
	}
	// Both circumcentres of the two right triangles on the square split
	// coincide at (1,1).
	for _, vv := range verts {
		v := vv.(map[string]any)
		if d := v["x"].(float64) - 1; d > 1e-9 || d < -1e-9 {
			t.Fatalf("voronoi vertex x = %v, want 1", v["x"])
		}
		if d := v["y"].(float64) - 1; d > 1e-9 || d < -1e-9 {
			t.Fatalf("voronoi vertex y = %v, want 1", v["y"])
		}
	}
}

func TestRejectionScenarios(t *testing.T) {
	r := Router()
	cases := []struct {
		name     string
		rawBody  string
		wantCode string
	}{
		{
			"too few points",
			`{"points":[{"x":0,"y":0},{"x":1,"y":1}]}`,
			"TOO_FEW_POINTS",
		},
		{
			"duplicate point",
			`{"points":[{"x":0,"y":0},{"x":1,"y":0},{"x":0,"y":1},{"x":1,"y":0}]}`,
			"DUPLICATE_POINT",
		},
		{
			"collinear",
			`{"points":[{"x":0,"y":0},{"x":1,"y":1},{"x":2,"y":2},{"x":3,"y":3}]}`,
			"DEGENERATE_POINT_SET",
		},
		{
			"nan coordinate",
			`{"points":[{"x":0,"y":0},{"x":1,"y":0},{"x":0,"y":NaN}]}`,
			"INVALID_COORDINATE",
		},
		{
			"infinite coordinate",
			`{"points":[{"x":0,"y":0},{"x":1,"y":0},{"x":0,"y":1e999}]}`,
			"INVALID_COORDINATE",
		},
		{
			"malformed json",
			`{"points": [`,
			"INVALID_JSON",
		},
		{
			"missing points field",
			`{}`,
			"INVALID_JSON",
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			for _, path := range []string{"/delaunay", "/voronoi"} {
				status, out := doPost(t, r, path, tc.rawBody)
				if status == http.StatusOK {
					t.Fatalf("%s %s: expected error status, got 200", path, tc.name)
				}
				errBody, ok := out["error"].(map[string]any)
				if !ok {
					t.Fatalf("%s %s: missing error envelope: %v", path, tc.name, out)
				}
				if errBody["code"] != tc.wantCode {
					t.Fatalf("%s %s: code = %v, want %s", path, tc.name, errBody["code"], tc.wantCode)
				}
				if errBody["message"] == nil || errBody["message"] == "" {
					t.Fatalf("%s %s: empty error message", path, tc.name)
				}
			}
		})
	}
}
