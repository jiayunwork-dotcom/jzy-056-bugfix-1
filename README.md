# Delaunay / Voronoi Geometry Service

A correctness-focused HTTP service that turns a scattered set of 2D points into
a **Delaunay triangulation** (empty circumcircle property) and its dual
**Voronoi diagram**. No rendering, no UI — only the geometry kernel and a JSON
HTTP API.

The triangulator is an incremental **Bowyer–Watson** inserter over a virtual
super triangle, with numeric predicates that stay stable near cocircular
configurations and under rigid translation.

## Layout

```
cmd/server              HTTP entry point (Gin)
internal/geom           point parsing/validation, tolerant predicates, hull
  point.go              input validation & error codes
  predicates.go         orient2d / incircle via conservative interval
                        arithmetic + a translation-noise radius; strict and
                        fuzzy incircle variants; exact circumcentre
  hull.go               monotone-chain convex hull, polygon area, hull clearance
internal/triangulate
  triangle.go           triangle type + adaptive super triangle
  builder.go            Bowyer–Watson insertion, cavity extraction, CCW re-fan
  verify.go             independent mesh verifier (empty circle, winding,
                        overlap, area sum, Euler count)
internal/voronoi        Voronoi dual export (circumcentre vertices, dual
  voronoi.go            edges, unbounded hull rays)
internal/api            Gin handlers, request decoding, JSON DTOs
  sample.go             built-in 5x5 jittered grid review case
```

## Key correctness decisions

- **Counter-clockwise winding everywhere.** Cavity boundary edges are
  re-oriented explicitly so the inserted point lies on their left; the same
  convention is reused at every insertion, so boundary orientation can never
  flip.
- **Tolerant incircle predicate.** Determinants are evaluated through
  hand-rolled *conservative interval arithmetic* (every float operation widens
  a validated error interval). A determinant whose interval straddles zero is
  classified `ON` the circle. For real triangles an `ON` (cocircular) triangle
  is deleted only when the new point is inside it — this keeps the existing
  diagonal for exactly cocircular quadruples instead of oscillating between
  the two legal diagonals.
- **Translation invariance.** The point set is recentred on its first point
  before any predicate is run, and each recentred coordinate carries an
  explicit noise radius (`2·eps·|anchor|`) from float subtraction. Topology is
  therefore identical under arbitrary rigid translation; Voronoi vertices move
  with the same vector.
- **Strict predicates for the super triangle.** The virtual super vertices are
  extremely far away (their distance is chosen from the closest approach of
  any interior point to the convex hull so ghost circumcircles cannot bulge
  inward across a hull edge). Ghost triangles use the interval predicate
  *without* the relative fuzzy band, which would otherwise be swamped at that
  scale.
- **Ghost triangles are fully stripped** before any result is returned; the
  response also carries hull area vs. summed triangle area so containment can
  be checked on every call.

## API

`POST /delaunay`
```json
{ "points": [ {"x":0,"y":0}, {"x":1,"y":0}, {"x":1,"y":1}, {"x":0,"y":1} ] }
```
Returns index triples (CCW), the independent `delaunay` verdict, hull edge
count, hull area, summed triangle area and their delta.

`POST /voronoi` — same body. Returns:
- `vertices`: one per Delaunay triangle, located at its **exact circumcentre**
  (plus the source triangle);
- `edges`: dual edges joining the circumcentres of the two triangles sharing an
  interior Delaunay edge;
- `rays`: unbounded Voronoi edges dual to convex-hull Delaunay edges, as an
  origin vertex plus a unit outward direction.

`GET /sample` — the built-in review case: a 5×5 regular grid with a small
deterministic perturbation, with the Euler-predicted triangle count
(`T = 2N - H - 2 = 44`) and hull edge count pinned.

### Errors

All rejections return `{ "error": { "code", "message" } }` with an HTTP status:

| code | status | cause |
| --- | --- | --- |
| `TOO_FEW_POINTS` | 400 | fewer than 3 points |
| `DUPLICATE_POINT` | 400 | two points share exact coordinates |
| `INVALID_COORDINATE` | 400 | NaN / ±Infinity coordinate (or non-numeric) |
| `DEGENERATE_POINT_SET` | 422 | all points strictly collinear |
| `INVALID_JSON` | 400 | malformed body / missing `points` array |

## Build & run (Go 1.22)

```bash
go test ./...                 # full property + rejection test suite
go run ./cmd/server           # serves :8080 (override with PORT)
```

## Docker

```bash
docker build -t delaunay .
docker run --rm -p 8080:8080 delaunay
```

The image builds with a pinned `golang:1.22` toolchain and runs a static binary
on a distroless base.

## Tests

- empty circumcircle for every triangle (random clouds and the exact
  cocircular cases);
- triangle-area sum equals convex-hull area (checked in a recentred frame);
- identical topology under large translations, Voronoi vertices translating by
  the same vector;
- inserting one interior point increases the triangle count by exactly two;
- exactly cocircular four-point inputs never cross or oscillate (repeated
  runs), including after a large translation;
- the 5×5 jittered grid yields exactly 44 triangles (Euler relation);
- every rejection path returns the right code (too few, collinear, duplicate,
  NaN/Infinity, malformed JSON), for both geometry endpoints.
