package geom

import "math"

// Numeric predicate support.
//
// The orient2d / incircle determinants are evaluated through conservative
// interval arithmetic: every coordinate is an interval, and each floating
// point operation widens the interval by a rigorous (slightly pessimistic)
// rounding bound. If the resulting interval straddles zero the configuration
// is declared ON the line / circumcircle. Bowyer-Watson then treats ON as
// "bad", which is what makes nearly cocircular quadruples stable instead of
// flip-flopping with floating point jitter.
//
// In addition to pure round-off, coordinates carry an explicit noise radius:
// after recentring the whole set on the first input point, each local
// coordinate inherits at most 2*eps*|anchor| of error from the way the
// (possibly far-translated) original coordinates are represented. Feeding
// that radius into the intervals makes the topology invariant under rigid
// translation as long as the translation does not destroy point
// distinctness itself.
const (
	// epsFloat is the float64 unit round-off (2^-52).
	epsFloat = 2.220446049250313e-16
	// tolRel is the geometric "fuzzy zero" band: a determinant smaller than
	// tolRel times its natural coordinate scale counts as ON. 1e-12 is far
	// above round-off for well-scaled data but far below any real geometry.
	tolRel = 1e-12
	// opGuard is extra pessimism per fused operation (a few ulps) so the
	// hand-rolled error bounds stay conservative.
	opGuard = 4.0
)

// TranslationNoiseRadius bounds the per-coordinate error introduced when a
// point set far from the origin is recentred by subtracting anchor: for
// float64 subtraction the residual is at most 2*eps*|anchor|.
func TranslationNoiseRadius(anchor Point) float64 {
	return 2 * epsFloat * math.Max(math.Abs(anchor.X), math.Abs(anchor.Y))
}

// Sign is the sign of a predicate value.
type Sign int8

const (
	Negative Sign = -1 // strictly outside / clockwise
	On       Sign = 0  // on the line / circumcircle within tolerance
	Positive Sign = 1  // strictly inside / counter-clockwise
)

// interval is a validated bound [lo, hi] on a real value.
type interval struct {
	lo, hi float64
}

// exactInt wraps x with an explicit uncertainty radius.
func exactInt(x, radius float64) interval {
	return interval{x - radius, x + radius}
}

func (a interval) maxMag() float64 {
	return math.Max(math.Abs(a.lo), math.Abs(a.hi))
}

// widen grows the interval by beta times its own magnitude.
func (a interval) widen(beta float64) interval {
	m := a.maxMag()
	return interval{a.lo - beta*m, a.hi + beta*m}
}

// addInterval: floating point addition with round-off bound.
func addInterval(a, b interval) interval {
	lo, hi := a.lo+b.lo, a.hi+b.hi
	r := opGuard * epsFloat * math.Max(math.Abs(lo), math.Abs(hi))
	return interval{lo - r, hi + r}
}

// subInterval: floating point subtraction with round-off bound.
func subInterval(a, b interval) interval {
	lo, hi := a.lo-b.hi, a.hi-b.lo
	r := opGuard * epsFloat * math.Max(math.Abs(lo), math.Abs(hi))
	return interval{lo - r, hi + r}
}

// mulInterval: interval product, accounting for float64 round-off.
func mulInterval(a, b interval) interval {
	p := [4]float64{a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi}
	lo, hi := p[0], p[0]
	m := math.Abs(p[0])
	for _, v := range p[1:] {
		if v < lo {
			lo = v
		}
		if v > hi {
			hi = v
		}
		if av := math.Abs(v); av > m {
			m = av
		}
	}
	r := opGuard * epsFloat * m
	return interval{lo - r, hi + r}
}

// signOf classifies a final determinant interval.
func signOf(v interval, tolBand float64) Sign {
	if v.hi < -tolBand {
		return Negative
	}
	if v.lo > tolBand {
		return Positive
	}
	return On
}

// PredCtx evaluates predicates for a point set that has been recentred on a
// common anchor. Radius is the per-coordinate uncertainty inherited from the
// original (pre-recentre) coordinates.
type PredCtx struct {
	// Radius bounds the difference between a stored local coordinate and the
	// exact translated coordinate of the original input point.
	Radius float64
}

// iv turns a stored local coordinate into an interval.
func (c PredCtx) iv(x float64) interval {
	return exactInt(x, c.Radius+opGuard*epsFloat*math.Abs(x))
}

func (c PredCtx) pt(p Point) (x, y interval) {
	return c.iv(p.X), c.iv(p.Y)
}

// Orient2d returns the sign of the signed doubled area of (a,b,c).
// Positive means counter-clockwise.
func (c PredCtx) Orient2d(a, b, p Point) Sign {
	ax, ay := c.pt(a)
	bx, by := c.pt(b)
	px, py := c.pt(p)

	dx1 := subInterval(bx, ax)
	dy1 := subInterval(by, ay)
	dx2 := subInterval(px, ax)
	dy2 := subInterval(py, ay)

	det := subInterval(mulInterval(dx1, dy2), mulInterval(dy1, dx2))

	// Natural scale of an orient determinant is coordinate^2.
	m := max4(dx1.maxMag(), dy1.maxMag(), dx2.maxMag(), dy2.maxMag())
	return signOf(det, tolRel*m*m)
}

// InCircle returns the sign of the incircle determinant for triangle (a,b,c)
// and point p. The triangle MUST be counter-clockwise: then Positive means
// p is strictly inside the circumcircle, On means cocircular (within
// tolerance), Negative means strictly outside.
func (c PredCtx) InCircle(a, b, cc, p Point) Sign {
	return c.incircle(a, b, cc, p, true)
}

// InCircleStrict is InCircle WITHOUT the extra geometric "fuzzy" band: only
// the rigorous interval round-off/translation-error bound is applied. It is
// used for virtual super (ghost) triangles, whose enormous coordinates would
// otherwise be swamped by the relative fuzzy band intended for real cocircular
// quadruples.
func (c PredCtx) InCircleStrict(a, b, cc, p Point) Sign {
	return c.incircle(a, b, cc, p, false)
}

func (c PredCtx) incircle(a, b, cc, p Point, fuzzy bool) Sign {
	// Work relative to a: this keeps every term bounded by the local spread
	// instead of by absolute coordinates.
	bx := b.X - a.X
	by := b.Y - a.Y
	cx := cc.X - a.X
	cy := cc.Y - a.Y
	px := p.X - a.X
	py := p.Y - a.Y

	// Each difference is exact up to the inherited coordinate noise; build
	// intervals for them directly.
	ibx := c.iv(bx)
	iby := c.iv(by)
	icx := c.iv(cx)
	icy := c.iv(cy)
	ipx := c.iv(px)
	ipy := c.iv(py)

	bb := addInterval(mulInterval(ibx, ibx), mulInterval(iby, iby)) // |b-a|^2
	cq := addInterval(mulInterval(icx, icx), mulInterval(icy, icy)) // |c-a|^2
	pq := addInterval(mulInterval(ipx, ipx), mulInterval(ipy, ipy)) // |p-a|^2

	// Expanding the 4x4 incircle determinant after translating a to the
	// origin yields (CCW triangle):
	//   det = -|b|^2*cross(c,p) + |c|^2*cross(b,p) - |p|^2*cross(b,c)
	// Positive <=> p strictly inside the circumcircle.
	t1 := mulInterval(bb, subInterval(mulInterval(icy, ipx), mulInterval(icx, ipy)))
	t2 := mulInterval(cq, subInterval(mulInterval(ibx, ipy), mulInterval(iby, ipx)))
	t3 := mulInterval(pq, subInterval(mulInterval(iby, icx), mulInterval(ibx, icy)))

	det := addInterval(addInterval(t1, t2), t3)

	// Natural scale is coordinate^4 (three coordinate factors and one
	// squared-norm factor). In fuzzy mode add the geometric tolerance band so
	// near-cocircular real quadruples classify as ON; strict mode uses only
	// the rigorous round-off interval (band 0).
	m := max6(ibx.maxMag(), iby.maxMag(), icx.maxMag(), icy.maxMag(), ipx.maxMag(), ipy.maxMag())
	band := 0.0
	if fuzzy {
		band = tolRel * m * m * m * m
	}
	return signOf(det, band)
}

func max4(a, b, c, d float64) float64 {
	m := math.Max(a, b)
	m = math.Max(m, c)
	return math.Max(m, d)
}

func max6(a, b, c, d, e, f float64) float64 {
	m := max4(a, b, c, d)
	m = math.Max(m, e)
	return math.Max(m, f)
}

// Circumcenter returns the circumcentre of triangle (a,b,c). The result is in
// the same coordinate system as the inputs. Computation is done relative to a
// so it is stable for triangles far from the origin. Callers must ensure the
// triangle is non-degenerate (non-zero signed area).
func Circumcenter(a, b, c Point) Point {
	bx, by := b.X-a.X, b.Y-a.Y
	cx, cy := c.X-a.X, c.Y-a.Y
	bLen := bx*bx + by*by
	cLen := cx*cx + cy*cy
	d := 2 * (bx*cy - by*cx)
	ux := (bLen*cy - cLen*by) / d
	uy := (cLen*bx - bLen*cx) / d
	return Point{a.X + ux, a.Y + uy}
}

// OrientRaw is a plain float64 orient2d with a relative tolerance, used on
// original (non-recentred) coordinates during input validation and convex
// hull walks. Positive is counter-clockwise.
func OrientRaw(a, b, c Point) Sign {
	return rawOrient(a, b, c)
}

// rawOrient is the implementation.
func rawOrient(a, b, c Point) Sign {
	det := (b.X-a.X)*(c.Y-a.Y) - (b.Y-a.Y)*(c.X-a.X)
	m := math.Max(math.Max(math.Abs(b.X-a.X), math.Abs(b.Y-a.Y)),
		math.Max(math.Abs(c.X-a.X), math.Abs(c.Y-a.Y)))
	band := tolRel * m * m
	switch {
	case det > band:
		return Positive
	case det < -band:
		return Negative
	default:
		return On
	}
}
