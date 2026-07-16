#!/usr/bin/env python3
"""segment_bound.py -- certified lower bound for det DF along the witness
segment p + t(q-p), t in [0,1]: one-dimensional branch and bound in
outward-rounded binary64 interval arithmetic, with exact decimal parsing
of the witness data, 2-ulp guards on libm exp including underflow
handling, and the exact 1/4 bound where a preactivation interval
straddles zero. The acceptance threshold is verified at runtime to be at
least the exact decimal 10^-6 stated by the theorem. Reads A, B, c, p, q
from data/witness.json; the file's "claimed" block is audited by
CHECK.py, not by this script. Any failure exits 1 before a certified
line can print.
Run from the repository root: python code/segment_bound.py

Version 1.3 retains the validated arithmetic path introduced in version
1.2, including the corrected negative-term accumulation and runtime
verification of the acceptance threshold.
"""
import json, math, sys, time
from fractions import Fraction
from decimal import Decimal

print("segment_bound.py v1.3")

INF = math.inf
def down(x): return math.nextafter(x, -INF)
def up(x): return math.nextafter(x, INF)
def g2lo(x): return down(down(x))
def g2hi(x): return up(up(x))
def rat_iv(q):
    f = float(q); qf = Fraction(f)
    if qf > q: return (down(f), f)
    if qf < q: return (f, up(f))
    return (f, f)
def iadd(a, b): return (down(a[0]+b[0]), up(a[1]+b[1]))
def imul(a, b):
    ps = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    return (down(min(ps)), up(max(ps)))
def iexp(t):
    lo = 0.0 if t[0] <= -744 else g2lo(math.exp(t[0]))
    hi = math.pow(2.0, -1073) if t[1] <= -744 else g2hi(math.exp(t[1]))
    return (max(lo, 0.0), hi)
def sigp_iv(t):
    def at(x, side):
        ax = abs(x); e = iexp((-ax, -ax))
        den = imul(iadd((1.0, 1.0), e), iadd((1.0, 1.0), e))
        v = (e[0]/den[1], e[1]/den[0])
        return down(v[0]) if side == 0 else up(v[1])
    far = t[0] if abs(t[0]) > abs(t[1]) else t[1]
    lo = max(at(far, 0), 0.0)
    if t[0] <= 0.0 <= t[1]:
        hi = 0.25
    else:
        near = t[0] if abs(t[0]) < abs(t[1]) else t[1]
        hi = min(at(near, 1), 0.25)
    return (lo, hi)

w = json.load(open("data/witness.json"),
              parse_float=lambda s: Fraction(Decimal(s)),
              parse_int=lambda s: Fraction(int(s)))
A, B, c = w["A"], w["B"], w["c"]
p, q = w["witness_pair"]["p"], w["witness_pair"]["q"]
import itertools
mI = {}
for I in itertools.combinations(range(4), 2):
    dA = A[0][I[0]]*A[1][I[1]] - A[0][I[1]]*A[1][I[0]]
    dB = B[I[0]][0]*B[I[1]][1] - B[I[1]][0]*B[I[0]][1]
    mI[I] = dA*dB
mIv = {I: rat_iv(v) for I, v in mI.items()}
# t_i(s) = alpha_i s + beta_i, exact rationals, then 1-2 ulp float enclosures
alpha = [B[i][0]*(q[0]-p[0]) + B[i][1]*(q[1]-p[1]) for i in range(4)]
beta = [B[i][0]*p[0] + B[i][1]*p[1] + c[i] for i in range(4)]
aiv = [rat_iv(a) for a in alpha]; biv = [rat_iv(b) for b in beta]

# single computation path: every +,* outward-rounded per operation
# (iadd/imul), libm exp 2-ulp guarded with underflow handling (iexp),
# sigma' by unimodal near/far with the 1/4 straddle branch (sigp_iv);
# interval multiplication handles the minor signs by its corner cases.
def det_lower(sv):
    sp = [sigp_iv(iadd(imul(aiv[i], sv), biv[i])) for i in range(4)]
    acc = (0.0, 0.0)
    for I, miv in mIv.items():
        acc = iadd(acc, imul(miv, imul(sp[I[0]], sp[I[1]])))
    return acc[0]

# acceptance threshold: must be at least the exact decimal 10^-6 the
# theorem states. The binary64 literal 1e-6 lies below that rational,
# so bump one float upward and verify exactly at runtime.
TARGET = 1e-6
if Fraction(TARGET) < Fraction(1, 10**6):
    TARGET = math.nextafter(TARGET, INF)
if Fraction(TARGET) < Fraction(1, 10**6):
    print("FAIL: acceptance threshold below exact decimal 10^-6"); sys.exit(1)
WMIN = 5e-9       # hard refinement floor; reaching it below TARGET is failure
t0 = time.time()
stack = [(0.0, 1.0)]
nleaf = 0
while stack:
    a, b = stack.pop()
    lb = det_lower((a, b))
    if lb >= TARGET:
        nleaf += 1
        continue
    if (b - a) <= WMIN:
        print(f"FAIL: bound {lb!r} < TARGET at floor width on [{a},{b}]")
        sys.exit(1)
    m = (a + b) / 2
    stack.append((m, b)); stack.append((a, m))
print("certified: det DF(p+t(q-p)) >= 1e-06 for all t in [0,1]")
print(f"  acceptance threshold {TARGET!r}, runtime-verified >= exact decimal 10^-6")
print(f"leaves = {nleaf}, refinement floor {WMIN}, {time.time()-t0:.1f}s")
print("VERDICT: outward-rounded interval certificate over the full segment")
