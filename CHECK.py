#!/usr/bin/env python3
"""CHECK.py -- independent certificate checker (spec v4, FROZEN).
Pure standard library. Validates every load-bearing numerical claim of
the paper from data/witness.json + data/certificate.json.
Exit 0 CERTIFIED; exit 1 a CERT obligation verified false; exit 2
could-not-verify (schema, grammar, NaN/Inf, domain).
"""
import sys, json, math
print("CHECK.py v1.4 (spec v4)")
def _crash_hook(tp, val, tb):
    print(f"EXIT2 (could not verify): unexpected exception {tp.__name__}: {val}")
    sys.exit(2)
sys.excepthook = _crash_hook
from fractions import Fraction
from decimal import Decimal

def die2(msg):
    print(f"EXIT2 (could not verify): {msg}"); sys.exit(2)
FAILED = []
def fail(ob, msg):
    print(f"CERT-FAIL {ob}: {msg}"); FAILED.append(ob)

# ---------- exact decimal parsing (spec: JSON numeric semantics) ----------
def jload(path):
    try:
        return json.load(open(path), parse_float=lambda s: Fraction(Decimal(s)),
                         parse_int=lambda s: Fraction(int(s)),
                         parse_constant=lambda s: die2(f"non-finite JSON constant {s}"))
    except Exception as e:
        die2(f"parse {path}: {e}")

# ---------- interval layer (outward-rounded binary64) ----------
INF = math.inf
def _chk(lo, hi):
    if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi):
        die2(f"interval invariant violated: [{lo},{hi}]")
    return (lo, hi)
def down(x): return math.nextafter(x, -INF)
def up(x):   return math.nextafter(x,  INF)
def rat_to_iv(q):
    f = float(q)
    if not math.isfinite(f): die2("rational has no finite binary64 bound")
    qf = Fraction(f)
    if qf > q:   return _chk(down(f), f)
    if qf < q:   return _chk(f, up(f))
    return _chk(f, f)
def iadd(a, b): return _chk(down(a[0]+b[0]), up(a[1]+b[1]))
def isub(a, b): return _chk(down(a[0]-b[1]), up(a[1]-b[0]))
def imul(a, b):
    ps = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    return _chk(down(min(ps)), up(max(ps)))
def idiv(a, b):
    if b[0] <= 0.0 <= b[1]: die2("division by interval containing zero")
    ps = (a[0]/b[0], a[0]/b[1], a[1]/b[0], a[1]/b[1])
    return _chk(down(min(ps)), up(max(ps)))
def guard2(x):  # 2-ulp libm guard band
    return _chk(down(down(x)), up(up(x)))
def iexp(t):    # interval exp with underflow rule (spec)
    los = 0.0 if t[0] <= -744 else guard2(math.exp(t[0]))[0]
    his = math.pow(2.0, -1073) if t[1] <= -744 else guard2(math.exp(t[1]))[1]
    return _chk(max(los, 0.0), his)
def ilog_rat(q):  # log_interval for exact rational q > 0
    if q <= 0: die2("log of nonpositive rational")
    lo, hi = rat_to_iv(q)
    return _chk(guard2(math.log(lo))[0], guard2(math.log(hi))[1])
def isqrt_rat(q):
    if q < 0: die2("sqrt of negative rational")
    lo, hi = rat_to_iv(q)
    return _chk(max(guard2(math.sqrt(lo))[0], 0.0), guard2(math.sqrt(hi))[1])
def itrig_rat(theta, fn):
    """cos/sin at exact rational theta: input enclosure + Lipschitz-1
    widening + 2-ulp guard, clipped to [-1,1] (spec)."""
    lo, hi = rat_to_iv(theta)
    w = hi - lo
    v = fn(lo)
    a = down(down(down(v) - w)); b = up(up(up(v) + w))
    return _chk(max(a, -1.0), min(b, 1.0))
def iabs(a):
    if a[0] >= 0: return a
    if a[1] <= 0: return _chk(-a[1], -a[0])
    return _chk(0.0, max(-a[0], a[1]))
def sigp_iv(t):
    """sigma'-interval on t-interval; overflow-safe form; 1/4 branch."""
    def sp_lo(x):  # lower bound of sigma' at point |x| maximal
        ax = abs(x); e = iexp((-ax, -ax))
        den = imul(iadd((1.0, 1.0), e), iadd((1.0, 1.0), e))
        return idiv(e, den)[0]
    def sp_hi(x):
        ax = abs(x); e = iexp((-ax, -ax))
        den = imul(iadd((1.0, 1.0), e), iadd((1.0, 1.0), e))
        return idiv(e, den)[1]
    far = t[0] if abs(t[0]) > abs(t[1]) else t[1]
    lo = sp_lo(far)
    if t[0] <= 0.0 <= t[1]: hi = 0.25
    else:
        near = t[0] if abs(t[0]) < abs(t[1]) else t[1]
        hi = min(sp_hi(near), 0.25)
    return _chk(max(lo, 0.0), hi)

TAU_LO = Fraction(Decimal("6.283185307179586"))
TAU_HI = Fraction(Decimal("6.283185307179587"))

# ---------- load and schema ----------
w = jload("data/witness.json"); cert = jload("data/certificate.json")
try:
    A = w["A"]; B = w["B"]; c = w["c"]
    p = w["witness_pair"]["p"]; q = w["witness_pair"]["q"]
    R0 = cert["meta"]["R0"]; raw_total = cert["meta"]["total_leaves"]
    theta = cert["tail"]["theta"]; tree = cert["bnb"]["tree"]
    kills = cert["sectors"]["kills"]; survs = cert["sectors"]["survivors"]
except Exception as e:
    die2(f"schema: {e}")
def req_num(x, where):
    if isinstance(x, bool) or not isinstance(x, Fraction):
        die2(f"non-numeric or boolean where number required: {where}")
    return x
for r in A:
    for v in r: req_num(v, "A")
for r in B:
    for v in r: req_num(v, "B")
for v in c: req_num(v, "c")
for v in p + q: req_num(v, "witness_pair")
if len(A) != 2 or any(len(r) != 4 for r in A) or len(B) != 4 or \
   any(len(r) != 2 for r in B) or len(c) != 4 or len(p) != 2 or len(q) != 2:
    die2("dimensions")
req_num(R0, "meta.R0")
if not R0 > 0: die2("R0 not positive")
req_num(raw_total, "meta.total_leaves")
if raw_total.denominator != 1 or raw_total <= 0:
    die2("meta.total_leaves must be a positive integer")
total = raw_total.numerator
theta2 = []
for k, t in enumerate(theta):
    if isinstance(t, str):
        try: t = Fraction(Decimal(t))
        except Exception: die2(f"tail.theta[{k}] not a decimal")
    req_num(t, f"tail.theta[{k}]")
    theta2.append(t)
theta = theta2
if theta[0] < 0 or any(theta[i] >= theta[i+1] for i in range(len(theta)-1)):
    die2("theta grid not strictly increasing from >= 0")
if theta[-1] >= TAU_LO:
    die2("theta_last not provably below 2*pi (inside or above tau band)")

# ---------- O1 minors (RATIONAL) ----------
import itertools
mI = {}
for I in itertools.combinations(range(4), 2):
    dA = A[0][I[0]]*A[1][I[1]] - A[0][I[1]]*A[1][I[0]]
    dB = B[I[0]][0]*B[I[1]][1] - B[I[1]][0]*B[I[0]][1]
    mI[I] = dA*dB
negs = [I for I, v in mI.items() if v < 0]
if len(negs) != 1 or negs[0] != (2, 3):
    fail("O1", f"negative-minor structure: {negs}")
else:
    cl = w.get("claimed", {})
    ok = True
    if "negative_minor_I" in cl and [int(x) for x in cl["negative_minor_I"]] != [2, 3]:
        ok = False; fail("O1", "claimed negative_minor_I mismatch")
    if "negative_minor_value" in cl:
        ex = mI[(2, 3)]; lo, hi = rat_to_iv(ex)
        for _ in range(4): lo, hi = down(lo), up(hi)   # 4-ulp containment
        if not (Fraction(lo) <= req_num(cl["negative_minor_value"], "claim") <= Fraction(hi)):
            ok = False; fail("O1", "claimed negative minor outside 4-ulp containment")
    if "minors" not in cl:
        die2("claimed.minors missing: all six minor products must be claimed")
    for I in itertools.combinations(range(4), 2):
        key = f"{I[0]}{I[1]}"
        if key not in cl["minors"]: die2(f"claimed minor {key} missing")
        claimed = req_num(cl["minors"][key], f"claimed minor {key}")
        lo, hi = rat_to_iv(mI[I])
        for _ in range(4): lo, hi = down(lo), up(hi)
        if not (Fraction(lo) <= claimed <= Fraction(hi)):
            ok = False; fail("O1", f"claimed minor {key} outside 4-ulp containment")
    if ok: print(f"CERT-PASS O1: unique negative minor (2,3) = {float(mI[(2,3)]):.12e}")
In, POS = (2, 3), [I for I in mI if mI[I] > 0]
if len(negs) != 1 or not POS:
    fail("O1", "minor sign structure precludes dependent obligations")
    print("NOT CERTIFIED: failed obligations ['O1'] (dependent obligations skipped)")
    sys.exit(1)

# ---------- O2 separation (INTERVAL) ----------
def t_of(x):   # exact rational preactivations
    return [B[i][0]*x[0] + B[i][1]*x[1] + c[i] for i in range(4)]
tp, tq = t_of(p), t_of(q)
def sig_iv_rat(t):  # sigma at exact rational t, interval
    lo, hi = rat_to_iv(t)
    def s(x):
        if x >= 0:
            e = iexp((-x, -x)); return idiv((1.0, 1.0), iadd((1.0, 1.0), e))
        e = iexp((x, x)); return idiv(e, iadd((1.0, 1.0), e))
    return _chk(s(lo)[0], s(hi)[1])
dd = []
for i in range(4):
    num = isub(sig_iv_rat(tq[i]), sig_iv_rat(tp[i]))
    den = rat_to_iv(tq[i] - tp[i])
    dd.append(idiv(num, den))
acc = (0.0, 0.0)
for I in mI:
    term = imul(imul(rat_to_iv(mI[I]), dd[I[0]]), dd[I[1]])
    acc = iadd(acc, term)
if acc[1] < 0:
    print(f"CERT-PASS O2: seg-avg det in [{acc[0]:.6e}, {acc[1]:.6e}] < 0")
else:
    fail("O2", f"upper bound {acc[1]:.3e} not < 0")
O2_iv = acc

# ---------- O4 tail (INTERVAL, normative formulas) ----------
KJ = {}
for J in POS:
    KJ[J] = iadd(ilog_rat(16*abs(mI[In])/mI[J]),
                 rat_to_iv(abs(c[J[0]]) + abs(c[J[1]]) + abs(c[In[0]]) + abs(c[In[1]])))
norms = [isqrt_rat(B[i][0]**2 + B[i][1]**2) for i in range(4)]
Lstar = iadd((max(iadd(norms[J[0]], norms[J[1]])[0] for J in POS),
              max(iadd(norms[J[0]], norms[J[1]])[1] for J in POS)),
             iadd(norms[In[0]], norms[In[1]]))
M = len(theta)
sp_next = [theta[(j+1) % M] - theta[j] for j in range(M-1)]
s_wrap_up = rat_to_iv(TAU_HI + theta[0] - theta[-1])[1]
sp_up = [rat_to_iv(s)[1] for s in sp_next] + [s_wrap_up]
def sj(j):  # upper bound of max adjacent spacing of node j
    return max(sp_up[j-1] if j > 0 else sp_up[-1], sp_up[j] if j < M else sp_up[-1])
bad = None
R0_lo = rat_to_iv(R0)[0]
for j in range(M):
    co = itrig_rat(theta[j], math.cos); si = itrig_rat(theta[j], math.sin)
    rho = []
    for i in range(4):
        d0 = imul(rat_to_iv(B[i][0]), co); d1 = imul(rat_to_iv(B[i][1]), si)
        rho.append(iabs(iadd(d0, d1)))
    rn = iadd(rho[In[0]], rho[In[1]])
    half = up(Lstar[1] * sj(j) / 2.0)
    feas = False
    for J in POS:
        rj = iadd(rho[J[0]], rho[J[1]])
        g = isub(rn, rj)
        slack = down(g[0] - half)
        if slack > 0 and KJ[J][1] <= down(slack * R0_lo):
            feas = True; break
    if not feas:
        bad = j; break
node_margin_lo = []
if bad is None and ("claimed_gap_margin" in cert["tail"] or
                    "claimed_node_margins" in cert["tail"]):
    for j in range(M):
        co = itrig_rat(theta[j], math.cos); si = itrig_rat(theta[j], math.sin)
        rho = []
        for i in range(4):
            rho.append(iabs(iadd(imul(rat_to_iv(B[i][0]), co),
                                 imul(rat_to_iv(B[i][1]), si))))
        rn = iadd(rho[In[0]], rho[In[1]])
        node_margin_lo.append(max(isub(rn, iadd(rho[J[0]], rho[J[1]]))[0] for J in POS))
if "claimed_global_margin" in cert["tail"]:
    die2("field renamed: use claimed_gap_margin (raw domination gap margin)")
if "claimed_gap_margin" in cert["tail"]:
    cg = cert["tail"]["claimed_gap_margin"]
    cg = Fraction(Decimal(cg)) if isinstance(cg, str) else req_num(cg, "claimed_gap_margin")
    if node_margin_lo and cg > Fraction(min(node_margin_lo)):
        fail("O4", f"claimed_gap_margin {float(cg):.6e} exceeds recomputed min raw gap margin {min(node_margin_lo):.6e}")
if "claimed_node_margins" in cert["tail"]:
    cns = cert["tail"]["claimed_node_margins"]
    if len(cns) != M: die2("claimed_node_margins length mismatch")
    for j, cv in enumerate(cns):
        cvf = Fraction(Decimal(cv)) if isinstance(cv, str) else req_num(cv, "node margin")
        if cvf > Fraction(node_margin_lo[j]):
            fail("O4", f"claimed node margin {j} exceeds recomputed"); break
if "claimed_L" in cert["tail"]:
    cL = cert["tail"]["claimed_L"]
    cL = Fraction(Decimal(cL)) if isinstance(cL, str) else req_num(cL, "claimed_L")
    if cL < Fraction(Lstar[0]):
        print(f"WARNING: claimed_L {float(cL):.6f} below checker L* lower bound {Lstar[0]:.6f} (non-verdict-bearing)")
if bad is None:
    if "O4" not in FAILED:
        print(f"CERT-PASS O4: tail feasible at all {M} nodes; L* <= {Lstar[1]:.4f}; R0 = {float(R0):.0f}")
else:
    fail("O4", f"node {bad} (theta={float(theta[bad]):.6f}) infeasible")

# ---------- O5 boxes (INTERVAL; normative replay geometry) ----------
nodes = iter(tree); nleaf = 0; mu = INF; worst = None
Bf = [[rat_to_iv(B[i][k]) for k in range(2)] for i in range(4)]
cf = [rat_to_iv(c[i]) for i in range(4)]
mIf = {I: rat_to_iv(v) for I, v in mI.items()}
stack = [(Fraction(-R0), Fraction(-R0), Fraction(R0), Fraction(R0))]
try:
    while stack:
        x1, y1, x2, y2 = stack.pop()
        node = next(nodes)
        if node == "X":
            mx = (x1 + x2) / 2
            stack.append((mx, y1, x2, y2)); stack.append((x1, y1, mx, y2))
            continue
        if node == "Y":
            my = (y1 + y2) / 2
            stack.append((x1, my, x2, y2)); stack.append((x1, y1, x2, my))
            continue
        try:
            claimed = Fraction(Decimal(node))  # leaf bound, preorder position
        except Exception:
            die2(f"invalid tree node {node!r}")
        xi = (rat_to_iv(x1)[0], rat_to_iv(x2)[1]); yi = (rat_to_iv(y1)[0], rat_to_iv(y2)[1])
        sps = []
        for i in range(4):
            ti = iadd(iadd(imul(Bf[i][0], xi), imul(Bf[i][1], yi)), cf[i])
            sps.append(sigp_iv(ti))
        lo = (0.0, 0.0)
        for I in mI:
            lo = iadd(lo, imul(imul(mIf[I], sps[I[0]]), sps[I[1]]))
        lb = lo[0]
        if not (lb > 0):
            fail("O5", f"leaf {nleaf} box=({float(x1):.4f},{float(y1):.4f},{float(x2):.4f},{float(y2):.4f}) bound {lb:.3e} not positive"); break
        if Fraction(lb) < claimed:
            fail("O5", f"leaf {nleaf}: claimed {float(claimed):.3e} > recomputed {lb:.3e}"); break
        if lb < mu: mu, worst = lb, (x1, y1, x2, y2)
        nleaf += 1
    else:
        try:
            next(nodes); die2("trailing tree nodes")
        except StopIteration:
            pass
        if nleaf != total: die2(f"leaf count {nleaf} != meta {total}")
        if "O5" not in FAILED:
            print(f"CERT-PASS O5: {nleaf} leaves verified; mu = {mu:.6e}")
except StopIteration:
    die2("tree exhausted before stack")

# ---------- O6 sectors (RATIONAL, exact) ----------
def sgn(x): return (x > 0) - (x < 0)
# exact 16-pattern base
dets = {}
for I in itertools.combinations(range(4), 2):
    dets[I] = B[I[0]][0]*B[I[1]][1] - B[I[1]][0]*B[I[0]][1]
if any(v == 0 for v in dets.values()) or any(B[i][0] == 0 and B[i][1] == 0 for i in range(4)):
    fail("O6", "degenerate arrangement")
dirs = []
for i in range(4):
    for s in (1, -1):
        dirs.append((s*B[i][1], -s*B[i][0], i))
def half(v):  # upper half-plane first; exact quadrant + cross ordering
    x, y, _ = v
    return 0 if (y > 0 or (y == 0 and x > 0)) else 1
dirs.sort(key=lambda v: (half(v),), reverse=False)
def cross(u, v): return u[0]*v[1] - u[1]*v[0]
import functools
def cmp(u, v):
    if half(u) != half(v): return -1 if half(u) < half(v) else 1
    cr = cross(u, v)
    return -1 if cr > 0 else (1 if cr < 0 else 0)
dirs = sorted(dirs, key=functools.cmp_to_key(cmp))
base = set()
for k in range(8):
    x, y, i = dirs[k]
    pat = tuple(0 if j == i else sgn(B[j][0]*x + B[j][1]*y) for j in range(4))
    base.add(pat)
    nx, ny, _ = dirs[(k+1) % 8]
    if (k+1) % 8 == 0: nx, ny = -nx, -ny     # wrap into next half-turn
    mxd, myd = x + nx, y + ny                # exact midpoint direction
    if mxd == 0 and myd == 0: fail("O6", "midpoint degenerate")
    base.add(tuple(sgn(B[j][0]*mxd + B[j][1]*myd) for j in range(4)))
    base.add(tuple(0 if j == i else -sgn(B[j][0]*x + B[j][1]*y) for j in range(4)))
    base.add(tuple(-sgn(B[j][0]*mxd + B[j][1]*myd) for j in range(4)))
if len(base) != 16:
    fail("O6", f"pattern base size {len(base)} != 16")
seen = set(); ok6 = "O6" not in FAILED
if len(kills) != 14 or len(survs) != 2:
    ok6 = False; fail("O6", f"sector counts {len(kills)} kills, {len(survs)} survivors (need 14+2)")
for kc in kills:
    sg = tuple(int(x) for x in kc["pattern"])
    if sg in seen:
        ok6 = False; fail("O6", f"duplicate pattern {sg}"); continue
    seen.add(sg)
    lam = [Fraction(Decimal(x)) for x in kc["lambda"]]
    Atl = [A[0][i]*lam[0] + A[1][i]*lam[1] for i in range(4)]
    raw = None
    for i in range(4):
        if sg[i] == 0: continue
        v = sg[i]*Atl[i]
        if v <= 0: ok6 = False; fail("O6", f"kill {sg}: sign condition fails at {i}"); break
        raw = v if raw is None else min(raw, v)
    else:
        cm = Fraction(Decimal(kc["claimed_margin"]))
        if cm*cm*(lam[0]**2 + lam[1]**2) > raw*raw:
            ok6 = False; fail("O6", f"kill {sg}: claimed margin too large")
det23 = A[0][2]*A[1][3] - A[0][3]*A[1][2]
for sv in survs:
    sg = tuple(int(x) for x in sv["pattern"])
    if sg in seen:
        ok6 = False; fail("O6", f"duplicate pattern {sg}"); continue
    seen.add(sg)
    g0, g1 = [Fraction(Decimal(x)) for x in sv["gamma01"]]
    g2 = (-(A[0][0]*g0 + A[0][1]*g1)*A[1][3] + (A[1][0]*g0 + A[1][1]*g1)*A[0][3]) / det23
    g3 = ( (A[0][0]*g0 + A[0][1]*g1)*A[1][2] - (A[1][0]*g0 + A[1][1]*g1)*A[0][2]) / det23
    gam = (g0, g1, g2, g3)
    if all(x == 0 for x in gam) or \
       A[0][0]*g0+A[0][1]*g1+A[0][2]*g2+A[0][3]*g3 != 0 or \
       A[1][0]*g0+A[1][1]*g1+A[1][2]*g2+A[1][3]*g3 != 0 or \
       tuple(sgn(x) for x in gam) != sg:
        ok6 = False; fail("O6", f"survivor {sg}: kernel realization fails")
if seen != base:
    ok6 = False; fail("O6", f"patterns mismatch: missing {base-seen}, extra {seen-base}")
if ok6 and "O6" not in FAILED:
    print("CERT-PASS O6: 14 kills exact, 2 survivors realized, 16-pattern base exact")

# ---------- O3 INFO ----------
mins = INF; minsU = INF
for k in range(201):
    a = Fraction(k, 200)
    x = (p[0] + a*(q[0]-p[0]), p[1] + a*(q[1]-p[1]))
    sps = []
    for i in range(4):
        ti = rat_to_iv(B[i][0]*x[0] + B[i][1]*x[1] + c[i])
        sps.append(sigp_iv(ti))
    lo = (0.0, 0.0)
    for I in mI:
        lo = iadd(lo, imul(imul(mIf[I], sps[I[0]]), sps[I[1]]))
    if lo[0] < mins: mins = lo[0]
    if lo[1] < minsU: minsU = lo[1]
print(f"INFO O3: min sampled det DF lower bound = {mins:.6e}  [SAMPLING IS NOT A CERTIFICATE]")
if "min_detDF_on_segment" in w.get("claimed", {}) or "item3" in w.get("claimed", {}):
    cm = w["claimed"].get("min_detDF_on_segment")
    if cm is not None:
        lo4, hi4 = mins, minsU
        for _ in range(4): lo4, hi4 = down(lo4), up(hi4)
        if not (Fraction(lo4) <= req_num(cm, "claimed min_detDF") <= Fraction(hi4)):
            fail("O3CLAIM", "claimed min_detDF outside sampled containment")
if "seg_avg_det" in w.get("claimed", {}):
    ca = req_num(w["claimed"]["seg_avg_det"], "claimed seg_avg")
    lo4, hi4 = O2_iv
    for _ in range(4): lo4, hi4 = down(lo4), up(hi4)
    if not (Fraction(lo4) <= ca <= Fraction(hi4)):
        fail("O2CLAIM", "claimed seg_avg_det outside containment")

# ---------- verdict ----------
if FAILED:
    print(f"NOT CERTIFIED: failed obligations {sorted(set(FAILED))}"); sys.exit(1)
print("CERTIFIED: all CERT obligations and claim audits verified")
sys.exit(0)
