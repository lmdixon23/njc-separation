"""Recompute the sector sign patterns and their explicit dual certificates.

Run from the repository root:
    python code/sector_lemma.py
"""
import json, itertools
import numpy as np

w = json.load(open("data/witness.json"))
A, B = np.array(w["A"]), np.array(w["B"])

print("=== Table 0: det(B_I) for all six row pairs ===")
for I in itertools.combinations(range(4), 2):
    print(f"  det B_{I} = {np.linalg.det(B[list(I), :]):+.12f}")

bnd = []
for i in range(4):
    th = (np.degrees(np.arctan2(B[i, 0], -B[i, 1]))) % 180
    bnd.append((th, i)); bnd.append((th + 180, i))
bnd.sort()
print("\n=== boundary angles (deg) ===")
for th, i in bnd:
    print(f"  {th:14.7f}   (b_{i} . u = 0)")

angles = [t for t, _ in bnd]
sectors, rays = [], []
print("\n=== Table 1a: open sectors ===")
for k in range(8):
    lo, hi = angles[k], angles[(k + 1) % 8] + (360 if k == 7 else 0)
    mid = np.radians((lo + hi) / 2)
    u = np.array([np.cos(mid), np.sin(mid)])
    s = tuple(int(np.sign(x)) for x in B @ u)
    sectors.append((lo, hi % 360, s))
    print(f"  ({lo:11.7f}, {hi % 360:11.7f})  sign = {s}")
print("\n=== Table 1b: boundary rays ===")
for th, i in bnd:
    u = np.array([np.cos(np.radians(th)), np.sin(np.radians(th))])
    v = B @ u
    s = tuple(0 if j == i else int(np.sign(v[j])) for j in range(4))
    rays.append((th, i, s))
    print(f"  {th:11.7f}  zero=b_{i}  sign = {s}")

print("\n=== U5: explicit kernel vectors for the surviving patterns ===")
T = -np.linalg.solve(A[:, [2, 3]], A[:, [0, 1]])
for g01 in [(-1.0, 1.0), (1.0, -1.0)]:
    g23 = T @ np.array(g01)
    g = np.array([g01[0], g01[1], g23[0], g23[1]])
    print(f"  gamma = ({g[0]:+.1f}, {g[1]:+.1f}, {g[2]:+.10f}, {g[3]:+.10f})"
          f"  sign = {tuple(int(np.sign(x)) for x in g)}"
          f"  |A@gamma| = {np.max(np.abs(A @ g)):.1e}")

print("\n=== Table 2: kill certificates ===")
allpats = [(f"sector ({lo:.4f},{hi:.4f})", s) for lo, hi, s in sectors] + \
          [(f"ray {th:.4f} (b_{i}=0)", s) for th, i, s in rays]
phis = np.linspace(0, 2 * np.pi, 720000, endpoint=False)
L = np.stack([np.cos(phis), np.sin(phis)], 1)
V = L @ A
for name, s in allpats:
    sv = np.array(s); nz = sv != 0
    margins = np.min(V[:, nz] * sv[nz][None, :], axis=1)
    j = int(np.argmax(margins))
    if margins[j] > 0:
        print(f"  KILLED   {name:32s} sign={s}"
              f"  lambda=({L[j,0]:+.6f},{L[j,1]:+.6f})  margin={margins[j]:.6f}")
    else:
        print(f"  SURVIVES {name:32s} sign={s}  (best margin {margins[j]:.3e})")

print("\n=== corollary: witness direction theta = 0 ===")
v0 = B @ np.array([1.0, 0.0])
print(f"  Bu_0 = {v0}")
print(f"  sign = {tuple(int(np.sign(x)) for x in v0)}")
