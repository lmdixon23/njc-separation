import itertools, time
import numpy as np
from scipy.special import expit
from scipy.optimize import minimize
rng=np.random.default_rng(7)
def sp(t): return expit(t)*expit(-t)
def dd(a,b):
    h=b-a; small=np.abs(h)<1e-8
    return np.where(small, sp((a+b)/2), (expit(b)-expit(a))/np.where(small,1.0,h))
def detADB(A,B,d):
    M=np.einsum('ij,gj,jk->gik',A,np.atleast_2d(d),B)
    return M[:,0,0]*M[:,1,1]-M[:,0,1]*M[:,1,0]
def detDF_grid(A,B,c,X): return detADB(A,B,sp(X@B.T+c))
def DF(A,B,c,x): d=sp(B@x+c); return (A*d[None,:])@B
def minors(A,B):
    return {I: np.linalg.det(A[:,I])*np.linalg.det(B[list(I),:])
            for I in itertools.combinations(range(A.shape[1]),2)}
def pattern(mp,tol=1e-10):
    v=np.array(list(mp.values())); p=np.any(v>tol); n=np.any(v<-tol)
    return 'mixed' if p and n else ('pos' if p else ('neg' if n else 'zero'))
G=np.linspace(-14,14,41); XG=np.array(np.meshgrid(G,G)).reshape(2,-1).T
def in_Nplus(A,B,c,nstarts=15):
    if detDF_grid(A,B,c,XG).min()<=0: return False
    f=lambda x: detDF_grid(A,B,c,x[None,:])[0]
    for _ in range(nstarts):
        r=minimize(f,rng.uniform(-20,20,2),method='Nelder-Mead',
                   bounds=[(-25,25)]*2,options=dict(maxiter=300,fatol=1e-16))
        if r.fun<0: return False
    th=np.linspace(0,2*np.pi,360,endpoint=False)
    for rad in (40.0,70.0):
        X=np.stack([rad*np.cos(th),rad*np.sin(th)],1)
        if np.any(detDF_grid(A,B,c,X)<0): return False
    return True
def newton_polish(A,B,c,p,q,iters=40):
    F=lambda x: A@expit(B@x+c)
    for _ in range(iters):
        J=DF(A,B,c,q)
        if np.linalg.cond(J)>1e8: return None
        try: q=q-np.linalg.solve(J,F(q)-F(p))
        except np.linalg.LinAlgError: return None
        if np.linalg.norm(q)>30: return None
    res=np.linalg.norm(F(q)-F(p))
    return (p,q,res) if (res<1e-13 and np.linalg.norm(p-q)>0.4) else None
def collision_hunt(A,B,c,nstarts=40,box=10.0,sep=0.8,seeds=None):
    F=lambda x: A@expit(B@x+c)
    def obj(z):
        p,q=z[:2],z[2:]; r=F(p)-F(q)
        return float(r@r+50*max(0.0,sep-np.linalg.norm(p-q))**2)
    starts=(list(seeds) if seeds else [])+[rng.uniform(-box,box,4) for _ in range(nstarts)]
    out=[]
    for z0 in starts:
        r=minimize(obj,z0,method='Nelder-Mead',bounds=[(-14,14)]*4,
                   options=dict(maxiter=3000,fatol=1e-22,xatol=1e-11))
        p,q=r.x[:2],r.x[2:]
        if np.linalg.norm(F(p)-F(q))<1e-9 and np.linalg.norm(p-q)>0.5*sep:
            pol=newton_polish(A,B,c,p,q)
            if pol is not None: out.append(pol)
    return out
def seg_min(A,B,c,nsamp=2500,nstarts=20,box=12.0):
    def qv(z):
        ta,tb=B@z[:2]+c,B@z[2:]+c
        return float(detADB(A,B,dd(ta,tb)[None,:])[0])
    Z=rng.uniform(-box,box,(nsamp,4)); vals=np.array([qv(z) for z in Z])
    best,bz=vals.min(),Z[np.argmin(vals)]
    for _ in range(nstarts):
        r=minimize(qv,rng.uniform(-box,box,4),method='Nelder-Mead',
                   bounds=[(-box,box)]*4,options=dict(maxiter=1200,fatol=1e-20))
        if r.fun<best: best,bz=r.fun,r.x
    return best,bz
def make_base(N0=3,eta=0.0):
    B0=rng.normal(size=(N0,2))
    R=rng.normal(size=(2,2))
    if np.linalg.det(R)<0: R[:,0]*=-1
    A0=R@B0.T+eta*rng.normal(size=(2,N0))
    return A0,B0,rng.normal(scale=1.5,size=N0)
print("=== E1: n=2, N=3 sanity ===")
acc=coll=tried=0; t0=time.time()
while acc<25 and tried<400:
    tried+=1; A0,B0,c0=make_base(3,rng.uniform(0,0.6))
    if not in_Nplus(A0,B0,c0,nstarts=10): continue
    acc+=1; coll+=len(collision_hunt(A0,B0,c0,nstarts=25))
print(f"  N+ tested: {acc}/{tried} tries; collisions: {coll} (theorem: 0) [{time.time()-t0:.0f}s]")
print("=== E3: n=2, N=4 frontier ===")
t0=time.time(); deep=[]; tried=0
while len(deep)<45 and tried<1500:
    tried+=1
    eta=rng.uniform(0.0,0.5); A0,B0,c0=make_base(3,eta)
    j=rng.integers(0,3); nu=rng.choice([0.0,0.05,0.2]); mu=rng.uniform(0.0,1.0)
    b4=B0[j]+nu*rng.normal(size=2)
    eps=10**rng.uniform(-1.7,-0.3)
    a4=-eps*(A0[:,j]+mu*rng.normal(size=2))
    A=np.concatenate([A0,a4[:,None]],1); B=np.vstack([B0,b4])
    c=np.concatenate([c0,[c0[j]+rng.uniform(-1,1)]])
    if pattern(minors(A,B))!='mixed': continue
    if not in_Nplus(A,B,c): continue
    deep.append((A,B,c,eta,mu,nu,eps))
print(f"  mixed-minor N+ networks: {len(deep)} (tried {tried}) [{time.time()-t0:.0f}s]")
t0=time.time(); cert=0; frontier=[]; collisions=[]; segmins=[]
for (A,B,c,eta,mu,nu,eps) in deep:
    smin,sz=seg_min(A,B,c); segmins.append(smin)
    if smin>0:
        cert+=1; f=collision_hunt(A,B,c,nstarts=20)
    else:
        frontier.append((A,B,c,smin)); f=collision_hunt(A,B,c,nstarts=120,seeds=[sz])
    if f: collisions.append((A,B,c,f[0]))
segmins=np.array(segmins)
print(f"  certified injective via segment-average positivity: {cert}/{len(deep)}")
print(f"  frontier (seg-avg det <=0): {len(frontier)}; collisions: {len(collisions)}")
if len(segmins): print(f"  seg-min stats: min={segmins.min():.3e} median={np.median(segmins):.3e}")
for (A,B,c,smin) in frontier[:6]: print(f"    frontier seg-min={smin:.3e}")
if collisions:
    A,B,c,(p,q,res)=collisions[0]; np.set_printoptions(precision=10)
    print("!! COUNTEREXAMPLE CANDIDATE !!\nA=",A,"\nB=",B,"\nc=",c,"\np=",p,"q=",q,"res=",res)
print(f"  [{time.time()-t0:.0f}s]")
import pickle
pickle.dump(dict(deep=deep,frontier=frontier,collisions=collisions,segmins=segmins),open('e3.pkl','wb'))
