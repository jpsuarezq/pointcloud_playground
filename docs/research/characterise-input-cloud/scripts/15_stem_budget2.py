"""Refit: per-stem point count vs stem diameter, using the woody cluster's own
size and the count of all band points inside its own radius."""
import numpy as np, os, sys, json
sys.path.insert(0, "/tmp/opencode/cloud-research/scripts")
from lib import load_xyz, load_attr, SCALE, OX, OY
from scipy.spatial import cKDTree
HC="/tmp/opencode/cloud-research/out/hag"; F="/tmp/opencode/cloud-research/out/fig"
G="/tmp/opencode/cloud-research/out/grids"; O="/tmp/opencode/cloud-research/out/dem"
m5=dict(l.split("=") for l in open(f"{G}/grid_5000.txt"))
R5,X5,Y5,NX5,NY5=int(m5["res"]),int(m5["x0"]),int(m5["y0"]),int(m5["nx"]),int(m5["ny"])
cnt5=np.load(f"{G}/cnt_5000.npy"); slope5=np.load(f"{O}/slope05.npy")
usable=(cnt5>=20)&np.isfinite(slope5)&(slope5<60)
X=load_xyz()
xr=np.asarray(X["x"],dtype=np.int64)*SCALE+OX; yr=np.asarray(X["y"],dtype=np.int64)*SCALE+OY
zr=np.asarray(X["z"],dtype=np.float64)*SCALE
xraw=np.asarray(X["x"],dtype=np.int64); yraw=np.asarray(X["y"],dtype=np.int64); del X
hag=np.load(f"{HC}/hag_smrf025.npy"); okh=np.load(f"{HC}/ok_smrf025.npy")
ix5=np.clip(((xraw-X5)//R5).astype(np.int32),0,NX5-1); iy5=np.clip(((yraw-Y5)//R5).astype(np.int32),0,NY5-1)
valid=usable[iy5,ix5]&okh
blk=(xr>=675360)&(xr<675420)&(yr>=758800)&(yr<758860)
bi=np.where(blk)[0]
P=np.column_stack([xr[bi],yr[bi],zr[bi]]); h=hag[bi]
t=cKDTree(P); K=16
dd,nb=t.query(P,k=K+1,workers=-1)
Q=P[nb]; Q-=Q.mean(1,keepdims=True)
cov=np.einsum("nki,nkj->nij",Q,Q)/K
ev,evec=np.linalg.eigh(cov); ev=np.clip(ev[:,::-1],1e-12,None)
lin=(ev[:,0]-ev[:,1])/ev[:,0]; vert=np.abs(evec[:,:,2][:,2])
band=(h>=1.15)&(h<=1.45)&valid[bi]
m=band&(lin>.6)&(vert>.8); idx=np.where(m)[0]
W=np.column_stack([xr[bi][idx],yr[bi][idx]])
pairs=cKDTree(W).query_pairs(0.10,output_type="ndarray")
par=np.arange(len(W))
def find(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
for a,b in pairs:
    ra,rb=find(a),find(b)
    if ra!=rb: par[ra]=rb
uq,lab=np.unique(np.array([find(i) for i in range(len(W))]),return_inverse=True)
sizes=np.bincount(lab); cent=np.array([W[lab==u].mean(0) for u in range(len(uq))])
BP=np.column_stack([xr[bi][band],yr[bi][band]]); t3=cKDTree(BP)
rows=[]
for u in range(len(uq)):
    if sizes[u]<6: continue
    own=W[lab==u]; dq=np.hypot(*(own-cent[u]).T)
    d=2*np.percentile(dq,75)
    if d<0.03 or d>1.2: continue
    alld=np.hypot(*(BP-cent[u]).T)
    rows.append((d,sizes[u],int(np.sum(alld<=d/2+0.03)),int(np.sum(alld<=d/2+0.10))))
R=np.array(rows)
print("clusters used:",len(R))
for j,lbl in [(1,"woody pts linked <0.10 m"),(2,"all band pts in r=d/2+3cm"),(3,"all band pts in r=d/2+10cm")]:
    A=np.polyfit(R[:,0],R[:,j],1); r=np.corrcoef(R[:,0],R[:,j])[0,1]
    print(f"{lbl:28s} n = {A[0]:7.1f}*d {A[1]:+6.1f}  R={r:.3f}  median={np.median(R[:,j]):.0f}")
print("\n  d(cm)  n(woody)  n(r+3cm)  n(r+10cm)   count")
for i in range(len(R)):
    pass
bins=[0.03,0.07,0.10,0.13,0.16,0.20,0.26,0.40,1.2]
for i in range(len(bins)-1):
    k=(R[:,0]>=bins[i])&(R[:,0]<bins[i+1])
    if k.sum()>=5:
        print(f"  {bins[i]*100:5.1f}-{bins[i+1]*100:5.1f} {np.median(R[k,1]):7.0f} {np.median(R[k,2]):9.0f} {np.median(R[k,3]):10.0f} {k.sum():7d}")
A=np.polyfit(R[:,0],R[:,2],1)
rho=A[0]/(np.pi*0.30)
print(f"\ncalibration from r=d/2+3cm counts: n = {A[0]:.1f}*d {A[1]:+.1f} -> rho_s = {rho:.0f} pts/m2")
print("\n  DBH   budget (fit)")
for d in (0.10,0.15,0.20,0.30,0.50,0.80):
    print(f"  {d*100:5.0f} cm {A[0]*d+A[1]:10.1f}")
json.dump(dict(rows=R.tolist(), fit=list(map(float,np.polyfit(R[:,0],R[:,2],1))), rho_s=float(rho)),
          open(f"{F}/stem_budget2.json","w"),indent=1)
