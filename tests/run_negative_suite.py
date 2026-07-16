#!/usr/bin/env python3
"""tests/run_negative_suite.py -- spec-mandated negative tests N1-N27.

Every test corrupts a copy of the real inputs (or constructs a synthetic
input where the spec requires one), runs CHECK.py as a SUBPROCESS in a
temporary directory, and asserts the genuine exit code plus, where
specified, an output marker. The suite never imports the checker's
verification logic; checker primitives are used only to CONSTRUCT
knife-edge inputs, never to judge them.

Status: 27/27 implemented. N21, N22, N27 are self-calibrating rounding
knife-edges: each locates the verdict threshold with replicas of the
checker's primitives, isolates ONE conservatism mechanism by simulating
its naive counterpart, tunes a decimal-literal input into the sliver
between the naive and guarded thresholds, and asserts the spec-compliant
checker refuses what the naive simulation would accept.

Usage: python tests/run_negative_suite.py [--only N1,N5,...]
       (run from the repository root)
"""
import sys, os, json, math, shutil, subprocess, tempfile, copy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def load(name):
    return json.load(open(os.path.join(ROOT, "data", name)))

W = load("witness.json")
C = load("certificate.json")

def run_check(w, c):
    d = tempfile.mkdtemp()
    os.mkdir(os.path.join(d, "data"))
    json.dump(w, open(os.path.join(d, "data", "witness.json"), "w"))
    json.dump(c, open(os.path.join(d, "data", "certificate.json"), "w"))
    shutil.copy(os.path.join(ROOT, "CHECK.py"), d)
    r = subprocess.run([sys.executable, "CHECK.py"], cwd=d,
                       capture_output=True, text=True, timeout=1800)
    shutil.rmtree(d, ignore_errors=True)
    LAST["out"] = r.stdout + r.stderr
    return r.returncode, LAST["out"]

def run_raw(wtxt, ctxt):
    d = tempfile.mkdtemp()
    os.mkdir(os.path.join(d, "data"))
    open(os.path.join(d, "data", "witness.json"), "w").write(wtxt)
    open(os.path.join(d, "data", "certificate.json"), "w").write(ctxt)
    shutil.copy(os.path.join(ROOT, "CHECK.py"), d)
    r = subprocess.run([sys.executable, "CHECK.py"], cwd=d,
                       capture_output=True, text=True, timeout=1800)
    shutil.rmtree(d, ignore_errors=True)
    LAST["out"] = r.stdout + r.stderr
    return r.returncode, LAST["out"]

BROKEN_PATCHES = {
    "log":  ("return _chk(guard2(math.log(lo))[0], guard2(math.log(hi))[1])",
             "return _chk(math.log(lo), math.log(hi))"),
    "trig": ("""    w = hi - lo
    v = fn(lo)
    a = down(down(down(v) - w)); b = up(up(up(v) + w))
    return _chk(max(a, -1.0), min(b, 1.0))""",
             """    v = fn(lo)
    return _chk(max(down(down(v)), -1.0), min(up(up(v)), 1.0))"""),
    "tau":  ('TAU_HI = Fraction(Decimal("6.283185307179587"))',
             'TAU_HI = Fraction(math.tau)'),
    "b64":  ("parse_float=lambda s: Fraction(Decimal(s)),",
             "parse_float=lambda s: Fraction(float(s)),"),
    "uexp0": ("his = math.pow(2.0, -1073) if t[1] <= -744 else guard2(math.exp(t[1]))[1]",
              "his = 0.0 if t[1] <= -744 else guard2(math.exp(t[1]))[1]"),
    "useL": ("half = up(Lstar[1] * sj(j) / 2.0)",
             'half = up((float(Fraction(Decimal(cert["tail"]["claimed_L"]))) if "claimed_L" in cert["tail"] else Lstar[1]) * sj(j) / 2.0)'),
}
def run_check_variant(wtxt, ctxt, variant):
    """Run a deliberately broken CHECK.py (one guard removed) on the same
    inputs; two-sided tests require the broken variant to reach the
    opposite verdict from the correct checker."""
    src = open(os.path.join(ROOT, "CHECK.py")).read()
    old, new = BROKEN_PATCHES[variant]
    if src.count(old) != 1:
        raise RuntimeError(f"broken-variant patch '{variant}' did not match exactly once")
    src = src.replace(old, new)
    d = tempfile.mkdtemp()
    os.mkdir(os.path.join(d, "data"))
    open(os.path.join(d, "data", "witness.json"), "w").write(wtxt)
    open(os.path.join(d, "data", "certificate.json"), "w").write(ctxt)
    open(os.path.join(d, "CHECK.py"), "w").write(src)
    r = subprocess.run([sys.executable, "CHECK.py"], cwd=d,
                       capture_output=True, text=True, timeout=1800)
    shutil.rmtree(d, ignore_errors=True)
    LAST["out"] = "[broken variant] " + r.stdout + r.stderr
    return r.returncode, r.stdout + r.stderr

def regen_tree(w, c, R0_float):
    """Regenerate the B&B tree for a retuned radius with emitter-replica
    logic, so O5 stays valid for both correct and broken variants; the
    discrimination must come from O4 alone."""
    B = w["B"]; cv = w["c"]; A = w["A"]
    import itertools as _it
    mp = {}
    for I in _it.combinations(range(4), 2):
        dA = A[0][I[0]]*A[1][I[1]] - A[0][I[1]]*A[1][I[0]]
        dB = B[I[0]][0]*B[I[1]][1] - B[I[1]][0]*B[I[0]][1]
        mp[I] = dA * dB
    pos_l = [(i, j, m) for (i, j), m in mp.items() if m > 0]
    In = [I for I, m in mp.items() if m < 0][0]; mneg = abs(mp[In])
    def sp(t):
        t = abs(t); e = math.exp(-min(t, 700)); return e / (1 + e) ** 2
    def bound(x1, y1, x2, y2):
        smin = [0.0]*4; smax = [0.0]*4
        for i in range(4):
            b0, b1 = B[i][0], B[i][1]
            tlo = b0*(x1 if b0 > 0 else x2) + b1*(y1 if b1 > 0 else y2) + cv[i]
            thi = b0*(x2 if b0 > 0 else x1) + b1*(y2 if b1 > 0 else y1) + cv[i]
            smin[i] = min(sp(tlo), sp(thi))
            smax[i] = 0.25 if tlo <= 0 <= thi else max(sp(tlo), sp(thi))
        lop = sum(m*smin[i]*smin[j] for (i, j, m) in pos_l)
        return lop - mneg*smax[In[0]]*smax[In[1]]
    tree = []; nleaf = 0
    stack = [(-R0_float, -R0_float, R0_float, R0_float)]
    while stack:
        x1, y1, x2, y2 = stack.pop()
        lb = bound(x1, y1, x2, y2)
        if lb > 0:
            tree.append(repr(0.0 if lb < 1e-100 else lb*1e-3)); nleaf += 1; continue
        wd, hd = x2-x1, y2-y1
        if wd < 0.008 and hd < 0.008:
            raise RuntimeError("regen_tree: nonpositive minimal box")
        if wd >= hd:
            tree.append("X"); mx = (x1+x2)/2
            stack.append((mx, y1, x2, y2)); stack.append((x1, y1, mx, y2))
        else:
            tree.append("Y"); my = (y1+y2)/2
            stack.append((x1, my, x2, y2)); stack.append((x1, y1, x2, my))
    c["bnb"]["tree"] = tree
    c["meta"]["total_leaves"] = nleaf

def cc():  # fresh deep copies
    return copy.deepcopy(W), copy.deepcopy(C)

RESULTS = []
LAST = {"out": ""}
def record(name, ok, note=""):
    RESULTS.append((name, ok, note))
    print(f"{name:5s} {'PASS' if ok else 'FAIL':4s} {note}")
    if not ok and LAST["out"]:
        for line in LAST["out"].splitlines():
            if ("EXIT2" in line or "Traceback" in line or "Error" in line
                    or "CERT-FAIL" in line or "NOT CERTIFIED" in line):
                print(f"      diag: {line.strip()[:160]}")
                break
        if os.environ.get("SUITE_VERBOSE"):
            print("      ---- full checker output ----")
            for line in LAST["out"].splitlines():
                print(f"      | {line[:200]}")

# ---------------- N1 leaf deleted ----------------
def n1():
    w, c = cc()
    t = c["bnb"]["tree"]
    for i in range(len(t) - 1, -1, -1):
        if t[i] not in ("X", "Y"):
            del t[i]; break
    code, out = run_check(w, c)
    record("N1", code == 2, f"exit={code}")

# ---------------- N2 R0 radius perturbation ----------------
# Reframed per review: v4 derives the O5 root box from meta.R0 itself,
# so no independent root exists to mismatch; the meaningful perturbation
# is a radius too small for tail feasibility.
def n2():
    w, c = cc()
    c["meta"]["R0"] = 300
    code, out = run_check(w, c)
    record("N2", code == 1 and "O4" in out, f"exit={code}")

# ---------------- N3 margin semantics: scale invariance + inflation ----------------
# Spec erratum found by this suite on first full run (2026-06-11): the
# squared-margin test claimed^2 * |lambda|^2 <= raw^2 is homogeneous of
# degree 2 in lambda, so scaling lambda with the margin retained keeps the
# claim TRUE; the sound verdict is exit 0, not the exit 1 the spec
# predicted. The attack the test actually guards is margin inflation.
def n3():
    w, c = cc()
    k = c["sectors"]["kills"][0]
    k["lambda"] = [repr(float(x) * 1e-12) for x in k["lambda"]]
    code1, _ = run_check(w, c)
    w, c = cc()
    k = c["sectors"]["kills"][0]
    k["claimed_margin"] = repr(float(k["claimed_margin"]) * 1.5)
    code2, out2 = run_check(w, c)
    ok = code1 == 0 and code2 == 1 and "O6" in out2
    record("N3", ok, f"scale:{code1} (invariant) inflate:{code2}")

# ---------------- N4 spacing gap opened at the critical direction ----------------
# Spec phrasing rotates the grid so the wraparound exceeds feasibility; an
# exact rotation needs exact 2 pi, which the decimal semantics prohibit, so
# the equivalent stressor deletes a 3000-node window centered on the
# minimum-gap direction: the adjacent spacing there exceeds what the
# half-spacing Lipschitz inequality can absorb. First smoke run showed a
# tail-truncation mutation is too weak (large gap at the wrap absorbs huge
# spacing soundly); this version targets the binding direction.
def n4():
    w, c = cc()
    th = [float(t) for t in c["tail"]["theta"]]
    B = w["B"]
    import itertools
    A = w["A"]
    mp = {}
    for I in itertools.combinations(range(4), 2):
        dA = A[0][I[0]]*A[1][I[1]] - A[0][I[1]]*A[1][I[0]]
        dB = B[I[0]][0]*B[I[1]][1] - B[I[1]][0]*B[I[0]][1]
        mp[I] = dA * dB
    POSl = [I for I, v in mp.items() if v > 0]
    gmin, jstar = float("inf"), 0
    for j, t in enumerate(th):
        u = (math.cos(t), math.sin(t))
        rho = [abs(B[i][0]*u[0] + B[i][1]*u[1]) for i in range(4)]
        rn = rho[2] + rho[3]
        g = max(rn - (rho[J[0]] + rho[J[1]]) for J in POSl)
        if g < gmin: gmin, jstar = g, j
    lo = max(1, jstar - 1500); hi = min(len(th) - 2, jstar + 1500)
    c["tail"]["theta"] = c["tail"]["theta"][:lo] + c["tail"]["theta"][hi:]
    if "claimed_gap_margin" in c["tail"]: del c["tail"]["claimed_gap_margin"]
    code, out = run_check(w, c)
    record("N4", code == 1 and "O4" in out, f"exit={code} (hole at j*={jstar}, gap={gmin:.5f})")

# ---------------- N5 + N25 claimed_L tiny: verdict invariant, warning present ----------------
def n5_n25():
    w, c = cc()
    c["tail"]["claimed_L"] = "1e-9"
    code, out = run_check(w, c)
    record("N5", code == 0, f"exit={code} (verdict invariance)")
    record("N25", "WARNING" in out and "claimed_L" in out, "warning line present" if "WARNING" in out else "warning MISSING")
    # stressed half: sparse grid makes O4 fail under the checker's L*;
    # a broken checker consuming claimed_L=1e-9 must falsely certify.
    w, c = cc()
    c["tail"]["theta"] = c["tail"]["theta"][::40]
    c["tail"]["claimed_L"] = "1e-9"
    if "claimed_gap_margin" in c["tail"]: del c["tail"]["claimed_gap_margin"]
    wtxt, ctxt = json.dumps(w), json.dumps(c)
    code2, out2 = run_raw(wtxt, ctxt)
    bcode, bout = run_check_variant(wtxt, ctxt, "useL")
    record("N5b", code2 == 1 and "CERT-FAIL O4" in out2 and bcode == 0,
           f"correct:{code2} broken[useL]:{bcode} (sparse grid; claimed_L must not rescue O4)")

# ---------------- N6 inflated claimed gap margin ----------------
def n6():
    # (a) unknown formula-like fields must not affect a valid verdict
    w, c = cc()
    c["tail"]["O4_formula"] = "exp(-rho*r)*bogus"
    c["tail"]["derivative_formula"] = "0"
    code_a, _ = run_check(w, c)
    # (b) inflated claimed margin fails even with formula fields present
    w, c = cc()
    c["tail"]["O4_formula"] = "exp(-rho*r)*bogus"
    c["tail"]["claimed_gap_margin"] = "0.5"
    code_b, out_b = run_check(w, c)
    record("N6", code_a == 0 and code_b == 1 and "O4" in out_b,
           f"formula-ignored:{code_a} inflate:{code_b}")

# ---------------- N7 straddling box, endpoint-monotone (wrong) claim ----------------
def n7():
    w, c = cc()
    from fractions import Fraction
    A = w["A"]; B = w["B"]; cv = w["c"]
    import itertools
    mp = {}
    for I in itertools.combinations(range(4), 2):
        dA = A[0][I[0]]*A[1][I[1]] - A[0][I[1]]*A[1][I[0]]
        dB = B[I[0]][0]*B[I[1]][1] - B[I[1]][0]*B[I[0]][1]
        mp[I] = dA * dB
    def sp(t):
        t = abs(t); e = math.exp(-min(t, 700)); return e / (1 + e) ** 2
    def bound(x1, y1, x2, y2, wrong):
        smin = []; smax = []
        for i in range(4):
            b0, b1 = B[i][0], B[i][1]
            tlo = b0*(x1 if b0 > 0 else x2) + b1*(y1 if b1 > 0 else y2) + cv[i]
            thi = b0*(x2 if b0 > 0 else x1) + b1*(y2 if b1 > 0 else y1) + cv[i]
            lo_ = min(sp(tlo), sp(thi))
            if (tlo <= 0 <= thi) and not wrong:
                hi_ = 0.25
            else:
                hi_ = max(sp(tlo), sp(thi))   # wrong: endpoint monotonicity
            smin.append(lo_); smax.append(hi_)
        lop = sum(m*smin[i]*smin[j] for (i, j), m in mp.items() if m > 0)
        upn = sum(-m*smax[i]*smax[j] for (i, j), m in mp.items() if m < 0)
        return lop - upn
    # replay tree to find a straddling leaf where the wrong bound exceeds the right one
    t = c["bnb"]["tree"]; R0 = c["meta"]["R0"]
    stack = [(-float(R0), -float(R0), float(R0), float(R0))]
    idx = -1; target = None
    for node in t:
        x1, y1, x2, y2 = stack.pop()
        if node == "X":
            mx = (x1 + x2) / 2
            stack.append((mx, y1, x2, y2)); stack.append((x1, y1, mx, y2)); continue
        if node == "Y":
            my = (y1 + y2) / 2
            stack.append((x1, my, x2, y2)); stack.append((x1, y1, x2, my)); continue
        idx += 1
        straddle = False
        for i in range(4):
            b0, b1 = B[i][0], B[i][1]
            tlo = b0*(x1 if b0 > 0 else x2) + b1*(y1 if b1 > 0 else y2) + cv[i]
            thi = b0*(x2 if b0 > 0 else x1) + b1*(y2 if b1 > 0 else y1) + cv[i]
            if tlo <= 0 <= thi: straddle = True; break
        if straddle:
            wrongb = bound(x1, y1, x2, y2, wrong=True)
            rightb = bound(x1, y1, x2, y2, wrong=False)
            if wrongb > rightb * 1.0001 and wrongb > 0:
                target = (idx, wrongb); break
    if target is None:
        record("N7", False, "PREMISE-UNREALIZED: no discriminating straddling leaf found"); return
    li = -1
    for k, node in enumerate(t):
        if node not in ("X", "Y"):
            li += 1
            if li == target[0]:
                t[k] = repr(target[1]); break
    code, out = run_check(w, c)
    record("N7", code == 1 and "O5" in out, f"exit={code} leaf {target[0]}")

# ---------------- N8 / N9 positive control: overflow + underflow regimes ----------------
def n8_n9():
    w, c = cc()
    code, out = run_check(w, c)
    okmu = False
    for line in out.splitlines():
        if "mu =" in line:
            try:
                mu = float(line.split("mu =")[1].strip().split()[0])
                okmu = 0 < mu < 1e-200
            except Exception:
                pass
    record("N8", code == 0, f"exit={code} (|t| up to ~1.9e3, finite bounds)")
    # N9: the 2^-1073 rule is structurally absorbed in any certificate
    # whose tail obligation passes (tail domination forces the negative
    # pair's saturation to outrun every positive floor), so no valid
    # certificate can make the rule verdict-decisive; see the gate
    # record. The probe asserts (a) the band is genuinely traversed,
    # (b) a checker variant with the band upper set to 0 reaches the
    # same mu to the last printed digit -- documenting the absorption
    # rather than pretending a discriminating oracle exists.
    okmu2 = code == 0 and okmu
    w2, c2 = cc()
    wtxt, ctxt = json.dumps(w2), json.dumps(c2)
    bcode, bout = run_check_variant(wtxt, ctxt, "uexp0")
    mu_b = None
    for line in bout.splitlines():
        if "mu =" in line:
            try: mu_b = line.split("mu =")[1].strip().split()[0]
            except Exception: pass
    mu_c = None
    for line in out.splitlines():
        if "mu =" in line:
            try: mu_c = line.split("mu =")[1].strip().split()[0]
            except Exception: pass
    record("N9", okmu2 and bcode == 0 and mu_b == mu_c,
           f"band traversed; absorption documented (mu correct={mu_c} uexp0={mu_b})")

# ---------------- N10 NaN in matrix ----------------
def n10():
    wtxt = json.dumps(W)
    first = W["A"][0][0]
    wtxt = wtxt.replace(json.dumps(first), "NaN", 1)
    code, out = run_raw(wtxt, json.dumps(C))
    record("N10", code == 2, f"exit={code}")

# ---------------- N11 truncated json ----------------
def n11():
    ctxt = json.dumps(C)[: len(json.dumps(C)) // 2]
    code, out = run_raw(json.dumps(W), ctxt)
    record("N11", code == 2, f"exit={code}")

# ---------------- N12 survivor pattern sign flipped ----------------
def n12():
    w, c = cc()
    sv = c["sectors"]["survivors"][0]
    sv["pattern"] = [-int(x) if int(x) != 0 else 0 for x in sv["pattern"]]
    # keep counts/uniqueness violation out of the way: flip makes it equal antipode ->
    # duplicate of survivor 2 OR pattern/gamma mismatch; either way exit 1 via O6
    code, out = run_check(w, c)
    record("N12", code == 1 and "O6" in out, f"exit={code}")

# ---------------- N13 theta grid violations (three sub-cases) ----------------
def n13():
    ok = True; notes = []
    for tag, mut in [
        ("unsorted", lambda th: th[:100] + [th[101], th[100]] + th[102:]),
        ("duplicate", lambda th: th[:100] + [th[100]] + th[100:]),
        ("ge2pi", lambda th: th[:-1] + ["6.2831853071795868"]),
    ]:
        w, c = cc()
        c["tail"]["theta"] = mut(c["tail"]["theta"])
        code, out = run_check(w, c)
        good = code == 2
        ok = ok and good; notes.append(f"{tag}:{code}")
    record("N13", ok, " ".join(notes))

# ---------------- N14 duplicate sector pattern, one omitted ----------------
def n14():
    w, c = cc()
    k = c["sectors"]["kills"]
    k[1] = copy.deepcopy(k[0])
    code, out = run_check(w, c)
    record("N14", code == 1 and "O6" in out, f"exit={code}")

# ---------------- N15 survivor gamma = 0 ----------------
def n15():
    w, c = cc()
    c["sectors"]["survivors"][0]["gamma01"] = ["0.0", "0.0"]
    code, out = run_check(w, c)
    record("N15", code == 1 and "O6" in out, f"exit={code}")

# ---------------- N16 exact-decimal vs binary64 semantics ----------------
def n16():
    # Synthetic A solved against the real B's pair determinants so that
    # BINARY64 semantics sees the canonical structure (unique negative
    # minor at (2,3); the (0,1) minor is exactly 0.0 in float) while
    # EXACT-DECIMAL semantics sees a second negative at (0,1) via
    # z = 1 + 1e-22 (rounds to 1.0). Correct checker: O1 structural
    # fail. b64-parsing variant: O1 passes. Downstream obligations on
    # the synthetic witness are noise for both sides and not asserted;
    # the isolation is the O1 marker differential.
    z = "1.0000000000000000000001"
    assert float(z) == 1.0
    import itertools as _it
    from fractions import Fraction as F
    from decimal import Decimal as D
    w, c = cc()
    B = w["B"]
    dB = {}
    for I in _it.combinations(range(4), 2):
        dB[I] = B[I[0]][0]*B[I[1]][1] - B[I[1]][0]*B[I[0]][1]
    if any(v == 0 for v in dB.values()):
        record("N16", False, "PREMISE-UNREALIZED: degenerate dB"); return
    s01 = 1 if dB[(0, 1)] > 0 else -1
    # rows: A = [[1, b, cc_, d], [zr, f, g, h]] with f = b so float a01 = 0;
    # exact a01 = f - b*z = b*(1 - z) = -b*1e-22; need exact m01 < 0:
    # sign(-b)*sign(dB01) < 0  =>  b = s01.
    b = s01; f = s01
    # a02 = g - cc_*z ~ g - cc_ = v2 ; a12 = b*g - cc_*f = s01*(g - cc_) = s01*v2
    # need m02 > 0: sign(v2) = sgn(dB02); need m12 > 0: s01*v2 sign = sgn(dB12)
    s02 = 1 if dB[(0, 2)] > 0 else -1
    s12 = 1 if dB[(1, 2)] > 0 else -1
    if s01 * s02 != s12:
        record("N16", False, "PREMISE-UNREALIZED: (02,12) sign constraint unsatisfiable"); return
    s03 = 1 if dB[(0, 3)] > 0 else -1
    s13 = 1 if dB[(1, 3)] > 0 else -1
    if s01 * s03 != s13:
        record("N16", False, "PREMISE-UNREALIZED: (03,13) sign constraint unsatisfiable"); return
    s23 = 1 if dB[(2, 3)] > 0 else -1     # need m23 < 0: sign(a23) = -s23
    # choose cc_, g with v2 = s02*1, and d, h with v3 = s03*1, then pick the
    # branch making a23 = cc_*h - d*g have sign -s23.
    sol = None
    for cc_ in (1.0, 2.0, -1.0, -2.0, 3.0, -3.0):
        for d in (1.0, 2.0, -1.0, -2.0, 3.0, -3.0):
            g = cc_ + s02; h = d + s03
            a23 = cc_*h - d*g
            if a23 != 0 and (1 if a23 > 0 else -1) == -s23:
                sol = (cc_, d, g, h); break
        if sol: break
    if sol is None:
        record("N16", False, "PREMISE-UNREALIZED: no a23 branch"); return
    cc_, d, g, h = sol
    A2 = [[1.0, float(b), cc_, d], [1.0, float(f), g, h]]   # 1.0 placeholder for z
    minors = {}
    for I in _it.combinations(range(4), 2):
        dA = A2[0][I[0]]*A2[1][I[1]] - A2[0][I[1]]*A2[1][I[0]]
        minors[f"{I[0]}{I[1]}"] = dA*dB[I]
    negs_f = sorted(k for k, v in minors.items() if v < 0)
    if not (minors["01"] == 0.0 and negs_f == ["23"]):
        record("N16", False, f"PREMISE-UNREALIZED: float signs {minors}"); return
    # exact check of m01 with z
    a01_exact = F(f) - F(b)*F(D(z))
    if not a01_exact*F(D(str(dB[(0, 1)]))) < 0:
        record("N16", False, "PREMISE-UNREALIZED: exact m01 not negative"); return
    w["A"] = A2
    w["claimed"]["minors"] = minors
    w["claimed"]["negative_minor_value"] = minors["23"]
    wtxt = json.dumps(w)
    needle = json.dumps(A2)
    target = needle.replace("[1.0, ", f"[{z}, ", 1) if needle.startswith('[[1.0') else None
    # inject z into A[0][0]? z multiplies row pairing via a01 = f - b*z meaning
    # z sits at A[1][0]: rebuild explicitly
    wtxt = wtxt.replace(f'[1.0, {float(f)}, {g}, {h}]', f'[{z}, {float(f)}, {g}, {h}]', 1)
    if z not in wtxt:
        record("N16", False, "PREMISE-UNREALIZED: z injection failed"); return
    ctxt = json.dumps(c)
    code, out = run_raw(wtxt, ctxt)
    okc = code == 1 and "negative-minor structure" in out
    bcode, bout = run_check_variant(wtxt, ctxt, "b64")
    okb = "CERT-PASS O1" in bout and "negative-minor structure" not in bout
    record("N16", okc and okb,
           f"correct:O1-structural-fail broken[b64]:O1-pass (exits {code}/{bcode})")

# ---------------- N17 trailing nodes / invalid axis ----------------
def n17():
    ok = True; notes = []
    w, c = cc(); c["bnb"]["tree"] = c["bnb"]["tree"] + ["X"]
    code, _ = run_check(w, c); ok &= code == 2; notes.append(f"trailing:{code}")
    w, c = cc()
    for i, n in enumerate(c["bnb"]["tree"]):
        if n == "X": c["bnb"]["tree"][i] = "Z"; break
    code, _ = run_check(w, c); ok &= code == 2; notes.append(f"axisZ:{code}")
    record("N17", ok, " ".join(notes))

# ---------------- N18 R0 nonpositive ----------------
def n18():
    ok = True; notes = []
    for v in (0, -5):
        w, c = cc(); c["meta"]["R0"] = v
        code, _ = run_check(w, c); ok &= code == 2; notes.append(f"R0={v}:{code}")
    record("N18", ok, " ".join(notes))

# ---------------- N19 claimed minor off by 6 ulps ----------------
def n19():
    w, c = cc()
    v = float(w["claimed"]["minors"]["01"])
    w["claimed"]["minors"]["01"] = v + 6 * math.ulp(v)
    code, out = run_check(w, c)
    record("N19", code == 1 and "O1" in out, f"exit={code}")

# ---------------- N20 all-zero minor products ----------------
def n20():
    w, c = cc()
    w["A"] = [[1, 1, 1, 1], [1, 1, 1, 1]]   # all A-minors zero
    code, out = run_check(w, c)
    record("N20", code == 1 and "O1" in out, f"exit={code}")


# ================= rounding knife-edges N21 / N22 / N27 =================
# Construction kernel: faithful replicas of CHECK.py's interval
# primitives, used ONLY to locate verdict thresholds; the verdict
# itself always comes from the CHECK.py subprocess. Naive variants
# simulate an implementation missing exactly one guard, isolating the
# mechanism under test. All tuned constants are emitted as decimal
# literals and all predicates are evaluated on the exact rationals of
# those literals, so the tests self-calibrate in each environment.
from fractions import Fraction as _F
from decimal import Decimal as _D
_INF = math.inf
def _dn(x): return math.nextafter(x, -_INF)
def _up(x): return math.nextafter(x, _INF)
def _riv(q):
    f = float(q); qf = _F(f)
    if qf > q: return (_dn(f), f)
    if qf < q: return (f, _up(f))
    return (f, f)
def _ia(a, b): return (_dn(a[0]+b[0]), _up(a[1]+b[1]))
def _is(a, b): return (_dn(a[0]-b[1]), _up(a[1]-b[0]))
def _im(a, b):
    ps = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    return (_dn(min(ps)), _up(max(ps)))
def _g2(x): return (_dn(_dn(x)), _up(_up(x)))
def _ilog(q, naive=False):
    lo, hi = _riv(q)
    if naive: return (math.log(lo), math.log(hi))
    return (_g2(math.log(lo))[0], _g2(math.log(hi))[1])
def _isq(q):
    lo, hi = _riv(q)
    return (max(_g2(math.sqrt(lo))[0], 0.0), _g2(math.sqrt(hi))[1])
def _itrig(th, fn, naive=False):
    lo, hi = _riv(th)
    if naive:
        v = fn(lo); return (max(_dn(_dn(v)), -1.0), min(_up(_up(v)), 1.0))
    w = hi - lo; v = fn(lo)
    return (max(_dn(_dn(_dn(v) - w)), -1.0), min(_up(_up(_up(v) + w)), 1.0))
def _iabs(a):
    if a[0] >= 0: return a
    if a[1] <= 0: return (-a[1], -a[0])
    return (0.0, max(-a[0], a[1]))
_TAU_LO = _F(_D("6.283185307179586")); _TAU_HI = _F(_D("6.283185307179587"))

def _o4_node_ok(Wd, theta_r, sj_up, R0_r, naive_log=False, naive_trig=False,
                naive_tau_sj=None):
    """Replicates CHECK.py's per-node feasibility decision; sj_up is the
    spacing upper bound (or naive_tau_sj when simulating raw math.tau)."""
    import itertools as _it
    A = Wd["A"]; B = Wd["B"]; cv = Wd["c"]
    mI = {}
    for I in _it.combinations(range(4), 2):
        dA = _F(_D(str(A[0][I[0]])))*_F(_D(str(A[1][I[1]]))) - _F(_D(str(A[0][I[1]])))*_F(_D(str(A[1][I[0]])))
        dB = _F(_D(str(B[I[0]][0])))*_F(_D(str(B[I[1]][1]))) - _F(_D(str(B[I[1]][0])))*_F(_D(str(B[I[0]][1])))
        mI[I] = dA*dB
    Inn = (2, 3); POSl = [I for I in mI if mI[I] > 0]
    Bq = [[_F(_D(str(B[i][k]))) for k in range(2)] for i in range(4)]
    cq = [_F(_D(str(cv[i]))) for i in range(4)]
    KJ = {}
    for J in POSl:
        KJ[J] = _ia(_ilog(16*abs(mI[Inn])/mI[J], naive=naive_log),
                    _riv(abs(cq[J[0]])+abs(cq[J[1]])+abs(cq[Inn[0]])+abs(cq[Inn[1]])))
    nm = [_isq(Bq[i][0]**2 + Bq[i][1]**2) for i in range(4)]
    Ls = _ia((max(_ia(nm[J[0]], nm[J[1]])[0] for J in POSl),
              max(_ia(nm[J[0]], nm[J[1]])[1] for J in POSl)),
             _ia(nm[Inn[0]], nm[Inn[1]]))
    co = _itrig(theta_r, math.cos, naive=naive_trig)
    si = _itrig(theta_r, math.sin, naive=naive_trig)
    rho = []
    for i in range(4):
        rho.append(_iabs(_ia(_im(_riv(Bq[i][0]), co), _im(_riv(Bq[i][1]), si))))
    rn = _ia(rho[Inn[0]], rho[Inn[1]])
    sj = naive_tau_sj if naive_tau_sj is not None else sj_up
    half = _up(Ls[1] * sj / 2.0)
    R0_lo = _riv(R0_r)[0]
    for J in POSl:
        g = _is(rn, _ia(rho[J[0]], rho[J[1]]))
        slack = _dn(g[0] - half)
        if slack > 0 and KJ[J][1] <= _dn(slack * R0_lo):
            return True
    return False

def _binding_R0_window(Wd, c, naive_log=False, naive_trig=False,
                       theta_override=None):
    """Binding node by float scan, then rational 28-digit refinement of
    the naive and spec R0 thresholds at that node."""
    from decimal import getcontext
    getcontext().prec = 40
    def q28(fr):
        return _F(_D(f"{_D(fr.numerator) / _D(fr.denominator):.28f}"))
    th = c["tail"]["theta"]
    M = len(th)
    th_r = [(_F(_D(t)) if isinstance(t, str) else _F(t)) for t in th]
    if theta_override:
        j0, tr = theta_override
        th_r[j0] = tr
    sp = [_riv(th_r[j+1] - th_r[j])[1] for j in range(M-1)]
    swrap = _riv(_TAU_HI + th_r[0] - th_r[-1])[1]
    spu = sp + [swrap]
    def sj(j): return max(spu[j-1] if j > 0 else spu[-1], spu[j] if j < M-1 else spu[-1])
    def ok(j, R0_r, naive):
        return _o4_node_ok(Wd, th_r[j], sj(j), R0_r,
                           naive_log=naive and naive_log,
                           naive_trig=naive and naive_trig)
    def fthresh(j):
        lo, hi = 1.0, 1000000.0
        for _ in range(80):
            mid = (lo + hi) / 2
            if ok(j, _F(_D(repr(mid))), False): hi = mid
            else: lo = mid
        return hi
    cand = list(range(0, M, max(1, M // 400)))
    if theta_override: cand.append(theta_override[0])
    best_j, best_t = None, -1.0
    for j in cand:
        t = fthresh(j)
        if t > best_t: best_t, best_j = t, j
    def rthresh(naive):
        lo, hi = _F(1), _F(1000000)
        for _ in range(140):
            mid = q28((lo + hi) / 2)
            if ok(best_j, mid, naive): hi = mid
            else: lo = mid
        return hi                      # smallest q28 grid point that passes
    return best_j, rthresh(True), rthresh(False), q28

def _knife_run(name, naive_log=False, naive_trig=False, theta_override=None,
               extra_note="", mutate=None, edge_checks=()):
    w, c = cc()
    if mutate: mutate(w, c)
    try:
        j, t_n, t_s, q28 = _binding_R0_window(w, c, naive_log=naive_log,
                                              naive_trig=naive_trig,
                                              theta_override=theta_override)
        if not (t_n < t_s):
            record(name, False,
                   f"PREMISE-UNREALIZED: naive thresh !< spec thresh at node {j} ({float(t_n)!r} vs {float(t_s)!r})")
            return
        mid = q28((t_n + t_s) / 2)     # naive-feasible, spec-infeasible
        if edge_checks:
            th_chk = c["tail"]["theta"]
            th_chk_r = [(_F(_D(t)) if isinstance(t, str) else _F(t)) for t in th_chk]
            if theta_override: th_chk_r[theta_override[0]] = theta_override[1]
            spx = [_riv(th_chk_r[k+1]-th_chk_r[k])[1] for k in range(len(th_chk_r)-1)]
            swx = _riv(_TAU_HI + th_chk_r[0] - th_chk_r[-1])[1]
            spxu = spx + [swx]
            def sjx(k): return max(spxu[k-1] if k > 0 else spxu[-1], spxu[k] if k < len(th_chk_r)-1 else spxu[-1])
            for k in edge_checks:
                if not _o4_node_ok(w, th_chk_r[k], sjx(k), mid):
                    record(name, False, f"PREMISE-UNREALIZED: hole edge node {k} infeasible at tuned R0"); return
        from decimal import getcontext
        getcontext().prec = 40
        lit = f"{_D(mid.numerator) / _D(mid.denominator):.28f}"
        if _F(_D(lit)) != mid:
            record(name, False, "PREMISE-UNREALIZED: R0 literal not exact"); return
        if theta_override:
            j0, tr = theta_override
            from decimal import getcontext
            getcontext().prec = 90
            tlit = f"{_D(tr.numerator) / _D(tr.denominator):.70f}".rstrip("0")
            if tlit.endswith("."): tlit += "0"
            if _F(_D(tlit)) != tr:
                record(name, False, "PREMISE-UNREALIZED: theta literal not exact"); return
            c["tail"]["theta"][j0] = tlit
        regen_tree(w, c, float(_F(_D(lit))))
        c["meta"]["R0"] = 123456.5     # magic placeholder
        ctxt = json.dumps(c).replace("123456.5", lit, 1)
        wtxt = json.dumps(w)
        code, out = run_raw(wtxt, ctxt)
        okc = code == 1 and "CERT-FAIL O4" in out
        bvar = "log" if naive_log else "trig"
        bcode, bout = run_check_variant(wtxt, ctxt, bvar)
        okb = bcode == 0 and "CERTIFIED" in bout
        record(name, okc and okb,
               f"correct:{code} broken[{bvar}]:{bcode} node {j} window ({float(t_n)!r},{float(t_s)!r}) {extra_note}")
    except Exception as e:
        record(name, False, f"HARNESS ERROR: {e}")

def n21():
    # Mechanism isolation for the log guard. Constraints discovered en
    # route (recorded for the spec ledger): (i) the guard survives into
    # K only via a large ln term, i.e. a small-minor binding pair;
    # (ii) R0's effective lattice is binary64 (rat_to_iv), so the
    # comparison resolves the guard only when the binding node's slack
    # is small (window ~ 2 ln/(ln+sum|c|) floats of R0); (iii) this
    # witness has no direction where the small-minor channel is both
    # exclusive and marginal, so marginality is manufactured: a hole in
    # the grid around the small-minor-exclusive direction inflates the
    # half-spacing term until slack ~ K/R0_target. Interior nodes are
    # strictly more feasible at the larger tuned R0 with unchanged
    # spacings, so the naive-pass premise needs checking only at the
    # center and the two hole edges (edge_checks).
    w0, c0 = cc()
    B = w0["B"]; A = w0["A"]
    import itertools as _it
    mp = {}
    for I in _it.combinations(range(4), 2):
        dA = A[0][I[0]]*A[1][I[1]] - A[0][I[1]]*A[1][I[0]]
        dB = B[I[0]][0]*B[I[1]][1] - B[I[1]][0]*B[I[0]][1]
        mp[I] = dA * dB
    big = [I for I, v in mp.items() if v > 0.1]
    small = [I for I, v in mp.items() if 0 < v <= 0.1]
    def gaps(t):
        u = (math.cos(t), math.sin(t))
        rho = [abs(B[i][0]*u[0] + B[i][1]*u[1]) for i in range(4)]
        rn = rho[2] + rho[3]
        gb = max(rn - rho[J[0]] - rho[J[1]] for J in big)
        gs = max(rn - rho[J[0]] - rho[J[1]] for J in small)
        return gs, gb
    best = (None, -1)
    for k in range(200000):
        t = 2*math.pi*k/200000
        gs, gb = gaps(t)
        if gb < -1e-3 and gs > best[1]: best = (t, gs)
    tstar, gs0 = best
    if tstar is None:
        record("N21", False, "PREMISE-UNREALIZED: no small-minor-exclusive direction"); return
    Lst = 5.83                                   # rough L*, refined by tuning loop
    for slack_target in (0.012, 0.018, 0.008, 0.025):
        # one-sided hole: the left flank cannot absorb wide spacing
        # (max gap ~0.6 at distance 0.6) while the right flank holds
        # ~2.05 throughout, so the spacing is manufactured entirely on
        # the right; the center keeps its ordinary left neighbor.
        D = 2.0 * (gs0 - slack_target) / Lst     # right-side hole width
        ge2 = gaps(tstar + D)
        need = Lst * D / 2 + 0.02
        if max(ge2) < need:
            continue                             # right edge could not absorb
        def mut(w, c, D=D):
            th = c["tail"]["theta"]
            thf = [float(t) for t in th]
            kept = [th[i] for i in range(len(th))
                    if thf[i] <= tstar - 1e-9 or thf[i] >= tstar + D]
            import bisect as _b
            pos = _b.bisect([float(t) for t in kept], tstar)
            kept.insert(pos, repr(float(tstar)))
            c["tail"]["theta"] = kept
            mut.center = pos
        mut.center = None
        w1, c1 = cc(); mut(w1, c1)
        ctr = mut.center
        _knife_run("N21", naive_log=True,
                   theta_override=(ctr, _F(float(tstar))),
                   mutate=lambda w, c: mut(w, c),
                   edge_checks=(mut.center - 1, mut.center + 1),
                   extra_note=f"(right hole width {D:.3f} at small-minor direction {tstar:.4f})")
        if RESULTS and RESULTS[-1][0] == "N21":
            if RESULTS[-1][1]: return
            RESULTS.pop()
    record("N21", False, "PREMISE-UNREALIZED: no slack target realized the guard window")

def _o4_precomp(w):
    import itertools as _it
    A = w["A"]; B = w["B"]; cv = w["c"]
    mI = {}
    for I in _it.combinations(range(4), 2):
        dA = _F(_D(str(A[0][I[0]])))*_F(_D(str(A[1][I[1]]))) - _F(_D(str(A[0][I[1]])))*_F(_D(str(A[1][I[0]])))
        dB = _F(_D(str(B[I[0]][0])))*_F(_D(str(B[I[1]][1]))) - _F(_D(str(B[I[1]][0])))*_F(_D(str(B[I[0]][1])))
        mI[I] = dA*dB
    Inn = (2, 3); POSl = [I for I in mI if mI[I] > 0]
    Bq = [[_F(_D(str(B[i][k]))) for k in range(2)] for i in range(4)]
    cq = [_F(_D(str(cv[i]))) for i in range(4)]
    KJ = {J: _ia(_ilog(16*abs(mI[Inn])/mI[J]),
                 _riv(abs(cq[J[0]])+abs(cq[J[1]])+abs(cq[Inn[0]])+abs(cq[Inn[1]])))
          for J in POSl}
    nm = [_isq(Bq[i][0]**2 + Bq[i][1]**2) for i in range(4)]
    Ls = _ia((max(_ia(nm[J[0]], nm[J[1]])[0] for J in POSl),
              max(_ia(nm[J[0]], nm[J[1]])[1] for J in POSl)),
             _ia(nm[Inn[0]], nm[Inn[1]]))
    return mI, Inn, POSl, Bq, KJ, Ls

def _node_ok_fast(pc, theta_r, sj_up, R0_r, naive_trig=False):
    mI, Inn, POSl, Bq, KJ, Ls = pc
    co = _itrig(theta_r, math.cos, naive=naive_trig)
    si = _itrig(theta_r, math.sin, naive=naive_trig)
    rho = [_iabs(_ia(_im(_riv(Bq[i][0]), co), _im(_riv(Bq[i][1]), si))) for i in range(4)]
    rn = _ia(rho[Inn[0]], rho[Inn[1]])
    half = _up(Ls[1] * sj_up / 2.0)
    R0_lo = _riv(R0_r)[0]
    for J in POSl:
        g = _is(rn, _ia(rho[J[0]], rho[J[1]]))
        slack = _dn(g[0] - half)
        if slack > 0 and KJ[J][1] <= _dn(slack * R0_lo):
            return True
    return False

def _node_thresh(pc, theta_r, sj_up, naive_trig=False):
    lo, hi = 1.0, 1000000.0
    for _ in range(70):
        mid = (lo + hi) / 2
        if _node_ok_fast(pc, theta_r, sj_up, _F(_D(repr(mid))), naive_trig): hi = mid
        else: lo = mid
    while _node_ok_fast(pc, theta_r, sj_up, _F(_D(repr(_dn(hi)))), naive_trig):
        hi = _dn(hi)
    return hi

def n22():
    # Grid thetas are exact binary64 values, so the input-enclosure
    # penalty exists ONLY at a non-float-exact theta; and every
    # direction has an antipodal twin with bitwise-identical thresholds.
    # The override therefore sits mid-ulp AT the binding direction,
    # making its threshold uniquely global by exactly the enclosure
    # amplification; the window floor is the surviving antipode (plus
    # local neighbors). Two-sided: the correct checker refuses inside
    # the window, the trig-naive variant certifies.
    w, c = cc()
    pc = _o4_precomp(w)
    th = c["tail"]["theta"]; M = len(th)
    th_r = [_F(_D(t)) for t in th]
    sp = [_riv(th_r[k+1]-th_r[k])[1] for k in range(M-1)]
    swrap = _riv(_TAU_HI + th_r[0] - th_r[-1])[1]
    spu = sp + [swrap]
    def sj(k): return max(spu[k-1] if k > 0 else spu[-1], spu[k] if k < M-1 else spu[-1])
    coarse = max(range(0, M, M//400),
                 key=lambda k: _node_thresh(pc, th_r[k], sj(k)))
    jb = max(range(max(0, coarse-220), min(M, coarse+220)),
             key=lambda k: _node_thresh(pc, th_r[k], sj(k)))
    f = float(th_r[jb])
    tr = _F(f) + (_F(_up(f)) - _F(f)) / 2          # mid-ulp: w = 1 ulp
    t_high = _node_thresh(pc, tr, sj(jb))
    rivals = list(range(max(0, jb-40), min(M, jb+40)))
    anti = (jb + M//2) % M
    rivals += list(range(max(0, anti-40), min(M, anti+40)))
    t_low = max(_node_thresh(pc, th_r[k], sj(k)) for k in rivals if k != jb)
    if not (t_low < t_high):
        record("N22", False,
               f"PREMISE-UNREALIZED: enclosure did not raise binding ({t_low!r} !< {t_high!r})")
        return
    from decimal import getcontext
    getcontext().prec = 40
    def q28(fr): return _F(_D(f"{_D(fr.numerator)/_D(fr.denominator):.28f}"))
    mid = q28((_F(_D(repr(t_low))) + _F(_D(repr(t_high)))) / 2)
    lit = f"{_D(mid.numerator)/_D(mid.denominator):.28f}"
    getcontext().prec = 90
    tlit = f"{_D(tr.numerator)/_D(tr.denominator):.70f}".rstrip("0")
    if _F(_D(lit)) != mid or _F(_D(tlit)) != tr:
        record("N22", False, "PREMISE-UNREALIZED: literal not exact"); return
    c["tail"]["theta"][jb] = tlit
    regen_tree(w, c, float(mid))
    c["meta"]["R0"] = 123456.5
    ctxt = json.dumps(c).replace("123456.5", lit, 1)
    wtxt = json.dumps(w)
    code, out = run_raw(wtxt, ctxt)
    okc = code == 1 and "CERT-FAIL O4" in out
    bcode, bout = run_check_variant(wtxt, ctxt, "trig")
    record("N22", okc and bcode == 0 and "CERTIFIED" in bout,
           f"correct:{code} broken[trig]:{bcode} node {jb} window ({t_low!r},{t_high!r}) (mid-ulp theta at binding direction)")

def n27():
    w, c = cc()
    th = c["tail"]["theta"]
    keep = int(len(th) * 0.80)            # flip point must sit inside range
    base = [_F(_D(t)) for t in th[:keep]]
    R0 = _F(c["meta"]["R0"])
    TAU_N = _F(math.tau)
    from decimal import getcontext
    getcontext().prec = 40
    def q28(fr):
        return _F(_D(f"{_D(fr.numerator) / _D(fr.denominator):.28f}"))
    def both_ok(tl, tau):
        sw = _riv(tau + base[0] - tl)[1]
        return (_o4_node_ok(w, base[0], sw, R0)
                and _o4_node_ok(w, tl, sw, R0))
    lo, hi = base[-1] + _F(1, 10**6), _TAU_LO - _F(1, 10**9)
    if both_ok(lo, _TAU_HI) or not both_ok(hi, _TAU_HI):
        record("N27", False,
               f"PREMISE-UNREALIZED: no spec flip in range (lo={both_ok(lo,_TAU_HI)}, hi={both_ok(hi,_TAU_HI)})")
        return
    for _ in range(140):                  # lo infeasible, hi feasible
        mid = q28((lo + hi) / 2)
        if both_ok(mid, _TAU_HI): hi = mid
        else: lo = mid
    # lo: spec-infeasible within ~1e-28 of the threshold; naive threshold
    # sits ~7.7e-16 lower, so lo must be naive-feasible.
    if both_ok(lo, _TAU_HI) or not both_ok(lo, TAU_N):
        record("N27", False, "PREMISE-UNREALIZED: tau window did not separate"); return
    lit = f"{_D(lo.numerator) / _D(lo.denominator):.28f}"
    if _F(_D(lit)) != lo:
        record("N27", False, "PREMISE-UNREALIZED: literal not exact"); return
    # dense bridge from the truncated grid to the tuned last node, original spacing
    bridge = []
    step = base[1] - base[0]
    t = base[-1] + step
    while t < lo - step:
        bridge.append(repr(float(t))); t += step
    c["tail"]["theta"] = th[:keep] + bridge + [lit]
    if "claimed_gap_margin" in c["tail"]: del c["tail"]["claimed_gap_margin"]
    wtxt, ctxt = json.dumps(w), json.dumps(c)
    code, out = run_raw(wtxt, ctxt)
    bcode, bout = run_check_variant(wtxt, ctxt, "tau")
    record("N27", code == 1 and "CERT-FAIL O4" in out and bcode == 0 and "CERTIFIED" in bout,
           f"correct:{code} broken[tau]:{bcode} (math.tau certifies, tau_hi refuses)")


# ---------------- N23 wrong-axis replay ----------------
def n23():
    w, c = cc()
    for i, n in enumerate(c["bnb"]["tree"]):
        if n in ("X", "Y"):
            c["bnb"]["tree"][i] = "Y" if n == "X" else "X"; break
    code, out = run_check(w, c)
    record("N23", code == 1 and "O5" in out, f"exit={code}")

# ---------------- N24 permuted child interpretation ----------------
def n24():
    w, c = cc()
    t = c["bnb"]["tree"]; idxs = [i for i, n in enumerate(t) if n not in ("X", "Y")]
    a, b = None, None
    for i in idxs:
        for j in idxs:
            if j > i and t[i] != t[j] and abs(float(t[i]) - float(t[j])) > 0.5 * max(float(t[i]), float(t[j])):
                a, b = i, j; break
        if a is not None: break
    t[a], t[b] = t[b], t[a]
    code, out = run_check(w, c)
    record("N24", code == 1 and "O5" in out, f"exit={code}")

# ---------------- N26 zero divided-difference denominator ----------------
def n26():
    w, c = cc()
    w["witness_pair"]["q"] = copy.deepcopy(w["witness_pair"]["p"])
    code, out = run_check(w, c)
    record("N26", code == 2, f"exit={code}")


def _one_level_cert():
    """Real witness, R0 = 0.25, one X split, two leaves with claims computed
    by emitter-replica bounds. O4 fails at this radius for both correct
    and mutated geometry (noise); the discriminator is the O5 failure
    marker: absent for the true subdivision, present for wrong-axis or
    swapped-children replay."""
    w, c = cc()
    B = w["B"]; cv = w["c"]; A = w["A"]
    import itertools as _it
    mp = {}
    for I in _it.combinations(range(4), 2):
        dA = A[0][I[0]]*A[1][I[1]] - A[0][I[1]]*A[1][I[0]]
        dB = B[I[0]][0]*B[I[1]][1] - B[I[1]][0]*B[I[0]][1]
        mp[I] = dA * dB
    pos_l = [(i, j, m) for (i, j), m in mp.items() if m > 0]
    In = [I for I, m in mp.items() if m < 0][0]; mneg = abs(mp[In])
    def sp(t):
        t = abs(t); e = math.exp(-min(t, 700)); return e / (1 + e) ** 2
    def bound(x1, y1, x2, y2):
        smin = [0.0]*4; smax = [0.0]*4
        for i in range(4):
            b0, b1 = B[i][0], B[i][1]
            tlo = b0*(x1 if b0 > 0 else x2) + b1*(y1 if b1 > 0 else y2) + cv[i]
            thi = b0*(x2 if b0 > 0 else x1) + b1*(y2 if b1 > 0 else y1) + cv[i]
            smin[i] = min(sp(tlo), sp(thi))
            smax[i] = 0.25 if tlo <= 0 <= thi else max(sp(tlo), sp(thi))
        return (sum(m*smin[i]*smin[j] for (i, j, m) in pos_l)
                - mneg*smax[In[0]]*smax[In[1]])
    bl = bound(-0.25, -0.25, 0.0, 0.25); br = bound(0.0, -0.25, 0.25, 0.25)
    if not (bl > 0 and br > 0 and abs(bl - br) > 1e-3 * max(bl, br)):
        return None, None, None
    c["meta"]["R0"] = 0.25
    c["meta"]["total_leaves"] = 2
    c["bnb"]["tree"] = ["X", repr(bl * 0.999), repr(br * 0.999)]
    return w, c, (bl, br)

def n23x():
    w, c, bb = _one_level_cert()
    if w is None:
        record("N23x", False, "PREMISE-UNREALIZED: halves not certifiable/asymmetric"); return
    base_code, base_out = run_check(w, json.loads(json.dumps(c)))
    c["bnb"]["tree"][0] = "Y"
    mut_code, mut_out = run_check(w, c)
    record("N23x", "CERT-FAIL O5" not in base_out and "CERT-FAIL O5" in mut_out,
           f"true-axis: no O5 fail; wrong-axis: O5 fail (exits {base_code}/{mut_code}, O4 noise tolerated)")

def n24x():
    w, c, bb = _one_level_cert()
    if w is None:
        record("N24x", False, "PREMISE-UNREALIZED: halves not certifiable/asymmetric"); return
    base_code, base_out = run_check(w, json.loads(json.dumps(c)))
    c["bnb"]["tree"][1], c["bnb"]["tree"][2] = c["bnb"]["tree"][2], c["bnb"]["tree"][1]
    mut_code, mut_out = run_check(w, c)
    record("N24x", "CERT-FAIL O5" not in base_out and "CERT-FAIL O5" in mut_out,
           f"lower-first replay normative: swapped children fail O5 (exits {base_code}/{mut_code})")

TESTS = {
    "N1": n1, "N2": n2, "N3": n3, "N4": n4, "N5": n5_n25, "N6": n6,
    "N7": n7, "N8": n8_n9, "N10": n10, "N11": n11, "N12": n12,
    "N13": n13, "N14": n14, "N15": n15, "N16": n16, "N17": n17,
    "N18": n18, "N19": n19, "N20": n20,
    "N21": n21, "N22": n22,
    "N23": n23, "N23x": n23x, "N24": n24, "N24x": n24x, "N26": n26,
    "N27": n27,
}
# N5 covers N25; N8 covers N9.

if __name__ == "__main__":
    only = None
    req_complete = "--require-complete" in sys.argv
    for a in sys.argv[1:]:
        if a.startswith("--only"):
            only = set(a.split("=", 1)[1].split(",")) if "=" in a else set(sys.argv[sys.argv.index(a) + 1].split(","))
    for name, fn in TESTS.items():
        if only and name not in only: continue
        LAST["out"] = ""
        try:
            fn()
        except Exception as e:
            record(name, False, f"HARNESS ERROR: {e}")
    bad = [n for n, ok, _ in RESULTS if not ok]
    print(f"\n{len(RESULTS)} rows, {len(RESULTS) - len(bad)} pass, {len(bad)} fail")
    sys.exit(1 if bad else 0)
