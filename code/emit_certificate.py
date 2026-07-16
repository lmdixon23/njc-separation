"""emit_certificate.py -- regenerate the positivity certificate AND emit
certificate.json for CHECK.py (spec v4). Mirrors certify.py's logic with
a deterministic preorder DFS so the checker can replay geometry exactly:
root [-R0,R0]^2, split wider side (tie -> x), exact dyadic midpoints.
"""
EMIT_VERSION = "v1.2 (claimed_gap_margin)"
print(f"emit_certificate {EMIT_VERSION}")
import json, math, time, itertools, pickle
import numpy as np

d = pickle.load(open('data/e3.pkl','rb'))
A, B, c, _ = d['frontier'][0]
N = 4
mp = {I: float(np.linalg.det(A[:, I]) * np.linalg.det(B[list(I), :]))
      for I in itertools.combinations(range(N), 2)}
neg = [I for I, v in mp.items() if v < 0]; pos = [I for I, v in mp.items() if v > 0]
assert len(neg) == 1
In = neg[0]; mneg = abs(mp[In])
bn = np.linalg.norm(B, axis=1)

# ---- tail grid (M nodes; coarser than original run, margins allow it) ----
M = 50000
th = np.linspace(0, 2*np.pi, M, endpoint=False)
U = np.stack([np.cos(th), np.sin(th)], 1)
proj = np.abs(U @ B.T)
rho_neg = proj[:, In[0]] + proj[:, In[1]]
psi = np.full(M, np.inf)
for J in pos:
    psi = np.minimum(psi, (proj[:, J[0]] + proj[:, J[1]]) - rho_neg)
Lpsi = max(bn[j0] + bn[j1] for (j0, j1) in pos) + bn[In[0]] + bn[In[1]]
dth = 2*np.pi / M
print(f"grid M={M}: max psi = {psi.max():.6f}; Lpsi = {Lpsi:.4f}; L*dth/2 = {Lpsi*dth/2:.2e}")
assert psi.max() + Lpsi*dth/2 < 0

# per-direction radius and R0 (as in certify.py)
Rdir = np.full(M, np.inf)
for J in pos:
    rhoJ = proj[:, J[0]] + proj[:, J[1]]
    KJ = math.log(16*mneg/mp[J]) + abs(c[J[0]]) + abs(c[J[1]]) + abs(c[In[0]]) + abs(c[In[1]])
    gap = rho_neg - rhoJ
    rJ = np.where(gap > 1e-12, KJ/np.maximum(gap, 1e-12), np.inf)
    Rdir = np.minimum(Rdir, rJ)
R0 = float(Rdir.max())
gapmin = float((-psi).min())
R0 *= gapmin / max(gapmin - Lpsi*dth/2, 1e-12)
R0 = int(np.ceil(R0) + 1)
print(f"R0 = {R0}")

# ---- interior B&B, deterministic preorder DFS with recorded tree ----
pos_l = [(j0, j1, mp[(j0, j1)]) for (j0, j1) in pos]
Bl = B.tolist(); cl = c.tolist()
def sp_minmax(lo, hi):
    # sigma' min at farthest-from-0, max at nearest (1/4 if straddling)
    if lo <= 0.0 <= hi: near = 0.0
    else: near = lo if abs(lo) < abs(hi) else hi
    far = lo if abs(lo) > abs(hi) else hi
    def sp(t):
        t = abs(t); e = math.exp(-t); return e/(1+e)**2
    return sp(far), sp(near)
def box_bound(x1, y1, x2, y2):
    smin = [0.0]*N; smax = [0.0]*N
    for i in range(N):
        b0, b1 = Bl[i][0], Bl[i][1]
        tlo = b0*(x1 if b0 > 0 else x2) + b1*(y1 if b1 > 0 else y2) + cl[i]
        thi = b0*(x2 if b0 > 0 else x1) + b1*(y2 if b1 > 0 else y1) + cl[i]
        smin[i], smax[i] = sp_minmax(tlo, thi)
    lopos = sum(m*smin[j0]*smin[j1] for (j0, j1, m) in pos_l)
    upneg = mneg*smax[In[0]]*smax[In[1]]
    return lopos - upneg
HMIN = 0.008
tree = []           # preorder: "X" / "Y" / float bound for leaf
stack = [(-float(R0), -float(R0), float(R0), float(R0))]
mu = float('inf'); nleaf = 0; t0 = time.time()
while stack:
    x1, y1, x2, y2 = stack.pop()
    lb = box_bound(x1, y1, x2, y2)
    if lb > 0:
        tree.append(0.0 if lb < 1e-100 else lb*1e-3)  # claims: 3 orders slack vs outward-rounding widening; vacuous-zero in deep saturation (claims are diagnostic; checker recomputes mu)
        mu = min(mu, lb); nleaf += 1
        continue
    w, h = x2-x1, y2-y1
    if w < HMIN and h < HMIN:
        raise SystemExit(f"FAILURE box at ({x1},{y1},{x2},{y2}) lb={lb}")
    if w >= h:
        tree.append("X"); mx = (x1+x2)/2
        stack.append((mx, y1, x2, y2)); stack.append((x1, y1, mx, y2))
    else:
        tree.append("Y"); my = (y1+y2)/2
        stack.append((x1, my, x2, y2)); stack.append((x1, y1, x2, my))
print(f"B&B: {nleaf} leaves, mu = {mu:.6e}, {time.time()-t0:.0f}s, tree nodes {len(tree)}")

# ---- sector certificates (reuse lemma generator logic) ----
def sectors():
    bnd = []
    for i in range(4):
        a = (math.degrees(math.atan2(B[i,0], -B[i,1]))) % 180
        bnd.append((a, i)); bnd.append((a+180, i))
    bnd.sort()
    angs = [a for a,_ in bnd]
    pats = []
    for k in range(8):
        lo, hi = angs[k], angs[(k+1) % 8] + (360 if k == 7 else 0)
        u = np.array([math.cos(math.radians((lo+hi)/2)), math.sin(math.radians((lo+hi)/2))])
        pats.append(("sector", tuple(int(np.sign(x)) for x in B@u)))
    for a, i in bnd:
        u = np.array([math.cos(math.radians(a)), math.sin(math.radians(a))])
        v = B@u
        pats.append(("ray", tuple(0 if j == i else int(np.sign(v[j])) for j in range(4))))
    return pats
phis = np.linspace(0, 2*np.pi, 720000, endpoint=False)
L = np.stack([np.cos(phis), np.sin(phis)], 1); V = L @ A
kills = []; survivors = []
T = -np.linalg.solve(A[:, [2, 3]], A[:, [0, 1]])
for kind, sg in sectors():
    sv = np.array(sg); nz = sv != 0
    marg = np.min(V[:, nz]*sv[nz][None, :], axis=1); j = int(np.argmax(marg))
    if marg[j] > 0:
        kills.append({"pattern": list(map(int, sg)),
                      "lambda": [repr(float(L[j,0])), repr(float(L[j,1]))],
                      "claimed_margin": repr(float(marg[j])*(1-1e-9))})
    else:
        g01 = (-1.0, 1.0) if sg[0] < 0 else (1.0, -1.0)
        g23 = T @ np.array(g01)
        survivors.append({"pattern": list(map(int, sg)),
                          "gamma01": [repr(g01[0]), repr(g01[1])]})
assert len(kills) == 14 and len(survivors) == 2

cert = {
  "meta": {"R0": R0, "total_leaves": nleaf, "format": "v4"},
  "tail": {"theta": [repr(float(t)) for t in th],
           "claimed_L": repr(float(Lpsi)),
           "claimed_gap_margin": repr(float(-psi.max())*(1-1e-6))},
  "bnb": {"tree": [n if isinstance(n, str) else repr(n) for n in tree]},
  "sectors": {"kills": kills, "survivors": survivors},
  "mu": repr(mu)
}
json.dump(cert, open('data/certificate.json', 'w'))
import os; print("certificate.json:", os.path.getsize('data/certificate.json')//1024, "KB")
