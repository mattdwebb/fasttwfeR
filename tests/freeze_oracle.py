"""Development-only frozen historical Python oracle, never a runtime dependency."""
import os
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
import argparse, importlib, json, pathlib, sys
import numpy as np
import pandas as pd
p=argparse.ArgumentParser()
p.add_argument("--source-root",required=True)
p.add_argument("--package",required=True)
args=p.parse_args()
sys.dont_write_bytecode=True
sys.path.insert(0,str(pathlib.Path(args.source_root)/args.package/"src"))
Design=importlib.import_module(args.package).Design
out=pathlib.Path(__file__).resolve().parent/"fixtures"
out.mkdir(exist_ok=True)
rng=np.random.default_rng(482031)
rows=[]
a=np.array([0,0,0,3,3,4,4,0.])
for s in range(8):
    for t in range(1,6):
        for j in range(6+s%3):
            x=rng.normal()+.03*t
            d=float(a[s]>0 and t>=a[s])
            y=.6*s+.2*t+1.2*x+d*(1.4+.03*s+.02*x)+rng.normal(scale=.4)
            rows.append((s,t,j,x,y,a[s],d, .5+rng.uniform()))
data=pd.DataFrame(rows,columns=["state","time","person","x","y","g","d","w"])
data["unit"]=data["state"]*10+data["person"]
data.to_csv(out/"data.csv",index=False)
assignments=np.array([np.roll(a,k) for k in (1,2,3,4)])
pd.DataFrame(assignments).to_csv(out/"assignments.csv",index=False,header=False)
for cov in (False,True):
    design=Design(data,outcome="y",cluster="state",time="time",assignment=a,covariates=["x"] if cov else [])
    n=design.native_cluster(engine="reference")
    theta=design.point(engine="reference")
    deleted=np.array([design.point(omit=i,engine="reference") for i in range(8)])
    jk=np.sqrt(7/8*np.sum((deleted-deleted.mean())**2))
    pd.DataFrame({"estimate":[theta],"native_se":[n["native_cluster_se"]],"native_p":[n["native_cluster_p"]],"jackknife_se":[jk]}).to_csv(out/("expected_x.csv" if cov else "expected.csv"),index=False)
    pd.DataFrame({"deleted":deleted}).to_csv(out/("deletions_x.csv" if cov else "deletions.csv"),index=False)
    rr=[]
    for draw in assignments:
        dd=Design(data,outcome="y",cluster="state",time="time",assignment=draw,covariates=["x"] if cov else [])
        nn=dd.native_cluster(engine="reference")
        val=dd.point(engine="reference")
        dels=np.array([dd.point(omit=i,engine="reference") for i in range(8)])
        se=np.sqrt(7/8*np.sum((dels-dels.mean())**2))
        rr.append((val,val/nn["native_cluster_se"],val/se))
    pd.DataFrame(rr,columns=["ri_beta","ri_t_native"]).to_csv(out/("randomized_x.csv" if cov else "randomized.csv"),index=False)
print("frozen historical fixture generated:",args.package)
