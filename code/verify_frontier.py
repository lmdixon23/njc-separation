import pickle, numpy as np, itertools
from scipy.special import expit
from scipy.optimize import minimize
rng=np.random.default_rng(99)
d=pickle.load(open('data/e3.pkl','rb'))
A,B,c,smin=d['frontier'][0]
np.set_printoptions(precision=8, suppress=False)
print("FRONTIER NETWORK\nA=",A,"\nB=",B,"\nc=",c,"\nreported seg-min=",smin)
def sp(t): return expit(t)*expit(-t)
def dd(a,b):
    h=b-a; small=np.abs(h)<1e-8
    return np.where(small, sp((a+b)/2), (expit(b)-expit(a))/np.where(small,1.0,h))
def detADB(A,B,dv):
    M=(A*dv[None,:])@B; return M[0,0]*M[1,1]-M[0,1]*M[1,0]
def detDF(x): return detADB(A,B,sp(B@x+c))
F=lambda x: A@expit(B@x+c)
mp={I: np.linalg.det(A[:,I])*np.linalg.det(B[list(I),:]) for I in itertools.combinations(range(4),2)}
print("minor products:", {k: f"{v:.4e}" for k,v in mp.items()})
# 1) re-verify N+ membership hard
g=np.linspace(-20,20,81); X=np.array(np.meshgrid(g,g)).reshape(2,-1).T
gridmin=min(detDF(x) for x in X)
worst=np.inf
for _ in range(200):
    r=minimize(lambda x: detDF(x), rng.uniform(-22,22,2), method='Nelder-Mead',
               bounds=[(-25,25)]*2, options=dict(maxiter=500,fatol=1e-18))
    worst=min(worst,r.fun)
th=np.linspace(0,2*np.pi,720,endpoint=False)
ringmin=min(min(detDF(np.array([r0*np.cos(t),r0*np.sin(t)])) for t in th) for r0 in (30.,50.,80.))
print(f"N+ recheck: gridmin={gridmin:.4e}, optmin={worst:.4e}, ringmin={ringmin:.4e}")
# 2) recover violating pair, characterize it
def qv(z):
    ta,tb=B@z[:2]+c,B@z[2:]+c
    return detADB(A,B,dd(ta,tb))
best=np.inf; bz=None
Z=rng.uniform(-12,12,(6000,4))
for z in Z:
    v=qv(z)
    if v<best: best,bz=v,z
for _ in range(60):
    r=minimize(qv,rng.uniform(-12,12,4),method='Nelder-Mead',
               bounds=[(-12,12)]*4,options=dict(maxiter=2000,fatol=1e-22))
    if r.fun<best: best,bz=r.fun,r.x
p,q=bz[:2],bz[2:]
print(f"seg-avg det min={best:.6e} at p={p}, q={q}, |p|={np.linalg.norm(p):.2f}, |q|={np.linalg.norm(q):.2f}, |p-q|={np.linalg.norm(p-q):.2f}")
# det DF along the segment (must be positive if N+ holds)
ts=np.linspace(0,1,201); segdets=[detDF(p+t*(q-p)) for t in ts]
print(f"det DF along segment: min={min(segdets):.4e} (positive => integrand dets all positive while averaged det is negative)")
# kernel alignment at a zero crossing of seg-avg det along the line from (p,p) to (p,q)? simpler: find s in [0,1] with qv at (p, p+s(q-p)) = 0
def qs(s):
    z=np.concatenate([p,p+s*(q-p)]); return qv(z)
lo,hi=1e-6,1.0
vlo,vhi=qs(lo),qs(hi)
print(f"qs(0+)={vlo:.3e}, qs(1)={vhi:.3e}")
if vlo*vhi<0:
    for _ in range(80):
        mid=(lo+hi)/2
        if qs(mid)*vlo>0: lo=mid
        else: hi=mid
    s0=(lo+hi)/2; q0=p+s0*(q-p)
    ta,tb=B@p+c,B@q0+c
    M=(A*dd(ta,tb)[None,:])@B
    _,_,Vt=np.linalg.svd(M); k=Vt[-1]
    dirpq=(q0-p)/np.linalg.norm(q0-p)
    print(f"singular avg-Jacobian at s0={s0:.6f}; |M v|={np.linalg.norm(M@dirpq):.3e}; angle(p-q0, ker)={np.degrees(np.arccos(min(1,abs(dirpq@k)))):.2f} deg")
# 3) brutal collision hunt
def hunt(nst, seeds=None):
    out=[]
    starts=(seeds or [])+[rng.uniform(-12,12,4) for _ in range(nst)]
    for z0 in starts:
        r=minimize(lambda z: float(np.sum((F(z[:2])-F(z[2:]))**2)+50*max(0,0.8-np.linalg.norm(z[:2]-z[2:]))**2),
                   z0, method='Nelder-Mead', bounds=[(-14,14)]*4,
                   options=dict(maxiter=4000,fatol=1e-24,xatol=1e-12))
        pp,qq=r.x[:2],r.x[2:]
        res=np.linalg.norm(F(pp)-F(qq))
        if res<1e-9 and np.linalg.norm(pp-qq)>0.4: out.append((pp,qq,res))
    return out
seeds=[np.concatenate([p,q]), np.concatenate([p,p+s0*(q-p)])] if vlo*vhi<0 else [np.concatenate([p,q])]
found=hunt(600, seeds)
print(f"collision hunt (600+ starts incl. seeded): {len(found)} candidates")
for f in found[:3]: print("   ", f)
