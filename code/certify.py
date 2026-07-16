"""Certificate that the frontier witness is in N+ (det DF > 0 on all of R^2).
Structure:
  TAIL LEMMA: using (1/4)e^{-|t|} <= sigma'(t) <= e^{-|t|}, for ||x||_2 > R0
  some positive Cauchy-Binet term dominates the single negative term, where
  R0 = max over directions u of min over positive pairs J of
       K_J(u)/(rho_neg(u)-rho_J(u))  whenever rho_J(u) < rho_neg(u),
       K_J = ln(16|m_neg|/m_J) + |c_j1|+|c_j2|+|c_2|+|c_3|.
  Validity requires psi(u) = min_J rho_J - rho_neg < 0 for ALL u (checked on a
  dense angular grid with an explicit Lipschitz bound in theta).
  INTERIOR: adaptive branch-and-bound on [-R0,R0]^2 with exact per-box
  enclosures of sigma' (unimodal: max at point nearest 0, min at farthest),
  comparison done in log domain via logsumexp. Box certified iff
  lower(sum of positive terms) > upper(|negative term|).
All floating point is IEEE double (certificate is rigorous modulo rounding).
"""
import pickle, time, itertools
import numpy as np
from scipy.special import logsumexp
d=pickle.load(open('data/e3.pkl','rb'))
A,B,c,_=d['frontier'][0]
n,N=A.shape
mp={I: np.linalg.det(A[:,I])*np.linalg.det(B[list(I),:]) for I in itertools.combinations(range(N),2)}
neg=[I for I,v in mp.items() if v<0]; pos=[I for I,v in mp.items() if v>0]
assert len(neg)==1
In=neg[0]; mneg=abs(mp[In])
print("negative pair:",In,"|m|=",f"{mneg:.6e}","positive pairs:",{I:f"{mp[I]:.4e}" for I in pos})
bn=np.linalg.norm(B,axis=1)

# ---- tail lemma ----
M=200000
th=np.linspace(0,2*np.pi,M,endpoint=False)
U=np.stack([np.cos(th),np.sin(th)],1)          # (M,2)
proj=np.abs(U@B.T)                              # (M,4) = |b_i . u|
rho_neg=proj[:,In[0]]+proj[:,In[1]]
psi=np.full(M,np.inf); R_need=np.zeros(M)
for J in pos:
    rhoJ=proj[:,J[0]]+proj[:,J[1]]
    KJ=np.log(16*mneg/mp[J])+abs(c[J[0]])+abs(c[J[1]])+abs(c[In[0]])+abs(c[In[1]])
    gap=rho_neg-rhoJ
    ok=gap>0
    rJ=np.where(ok,KJ/np.maximum(gap,1e-30),np.inf)
    better=rhoJ-rho_neg<psi
    psi=np.minimum(psi,rhoJ-rho_neg)
    R_need=np.where(rJ<np.inf,np.minimum(R_need+np.where(R_need==0,np.inf,0),rJ),R_need) if False else R_need
# per-direction minimal radius over admissible J:
Rdir=np.full(M,np.inf)
for J in pos:
    rhoJ=proj[:,J[0]]+proj[:,J[1]]
    KJ=np.log(16*mneg/mp[J])+abs(c[J[0]])+abs(c[J[1]])+abs(c[In[0]])+abs(c[In[1]])
    gap=rho_neg-rhoJ
    rJ=np.where(gap>1e-12,KJ/np.maximum(gap,1e-12),np.inf)
    Rdir=np.minimum(Rdir,rJ)
psimax=psi.max()
Lpsi=max(bn[J0]+bn[J1] for (J0,J1) in pos)+bn[In[0]]+bn[In[1]]
dth=2*np.pi/M
print(f"tail lemma: max_u psi = {psimax:.6f} (need <0); Lipschitz slack {Lpsi*dth/2:.2e}")
assert psimax + Lpsi*dth/2 < 0, "tail lemma fails"
R0=float(Rdir.max())
# Lipschitz correction on R0: Rdir is K/gap; gap Lipschitz<=Lpsi => inflate
gapmin=float((-psi).min())
R0*= gapmin/max(gapmin-Lpsi*dth/2,1e-12)
R0=np.ceil(R0)+1
print(f"R0 = {R0}")

# ---- interior B&B ----
logm_pos=np.array([np.log(mp[J]) for J in pos]); pos_idx=np.array(pos)
def logsp(t): t=np.abs(t); return -t-2*np.log1p(np.exp(-t))
def lo_hi_t(boxes):
    # boxes: (K,4) x1,y1,x2,y2 ; corners (K,4,2)
    x1,y1,x2,y2=boxes.T
    cs=np.stack([np.stack([x1,y1],1),np.stack([x1,y2],1),
                 np.stack([x2,y1],1),np.stack([x2,y2],1)],1)
    T=cs@B.T+c                                  # (K,4corners,4units)
    return T.min(1),T.max(1)
def certify_batch(boxes):
    lo,hi=lo_hi_t(boxes)                        # (K,4)
    far=np.where(np.abs(lo)>np.abs(hi),lo,hi)
    inside=(lo<=0)&(hi>=0)
    near=np.where(inside,0.0,np.where(np.abs(lo)<np.abs(hi),lo,hi))
    lsmin=logsp(far); lsmax=logsp(near)
    Lpos=logsumexp(logm_pos[None,:]+lsmin[:,pos_idx[:,0]]+lsmin[:,pos_idx[:,1]],axis=1)
    Uneg=np.log(mneg)+lsmax[:,In[0]]+lsmax[:,In[1]]
    return Lpos>Uneg
t0=time.time()
root=np.array([[-R0,-R0,R0,R0]],float)
queue=[root]; ncert=0; nproc=0; fails=[]
HMIN=0.008; BUDGET=8_000_000
while queue and nproc<BUDGET:
    boxes=np.concatenate(queue); queue=[]
    # process in chunks
    for s in range(0,len(boxes),250000):
        b=boxes[s:s+250000]; nproc+=len(b)
        ok=certify_batch(b); ncert+=ok.sum()
        bad=b[~ok]
        if len(bad)==0: continue
        w=bad[:,2]-bad[:,0]; h=bad[:,3]-bad[:,1]
        tiny=(w<HMIN)&(h<HMIN)
        if tiny.any(): fails.append(bad[tiny])
        bad=bad[~tiny]; w=bad[:,2]-bad[:,0]; h=bad[:,3]-bad[:,1]
        sx=w>=h
        bx=bad[sx]; mx=(bx[:,0]+bx[:,2])/2
        c1=np.stack([bx[:,0],bx[:,1],mx,bx[:,3]],1); c2=np.stack([mx,bx[:,1],bx[:,2],bx[:,3]],1)
        by=bad[~sx]; my=(by[:,1]+by[:,3])/2
        c3=np.stack([by[:,0],by[:,1],by[:,2],my],1); c4=np.stack([by[:,0],my,by[:,2],by[:,3]],1)
        for cc in (c1,c2,c3,c4):
            if len(cc): queue.append(cc)
nf=sum(len(f) for f in fails)
print(f"B&B: processed {nproc} boxes in {time.time()-t0:.0f}s; certified {ncert}; FAILURES at h<{HMIN}: {nf}; budget_exhausted={nproc>=BUDGET and bool(queue)}")
if nf:
    F=np.concatenate(fails); print("failure box sample:",F[:3], "extent x:",F[:,0].min(),F[:,2].max(),"y:",F[:,1].min(),F[:,3].max())
