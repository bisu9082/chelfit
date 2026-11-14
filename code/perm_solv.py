# -*- coding: utf-8 -*-
"""세트 내부 라벨 순열검정 — 용매화 기술자. QM/고전 양쪽, 동일 서브셋, 동일 파이프라인."""
import numpy as np,pandas as pd,sys,warnings,json
warnings.filterwarnings("ignore")
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import r2_score
sys.path.insert(0,"/home/claude/chelfit")
from descriptors import featurize,FEATURES
rng=np.random.default_rng(31337); NPERM=2000
df=pd.read_csv("raw_results.csv",encoding="utf-8-sig")
qs=pd.read_csv("qm_descriptors_solv.csv",encoding="utf-8-sig")
ex=pd.read_csv("exchange_energy.csv",encoding="utf-8-sig")[["chelator","metal","dE_exch_eV"]]
qs=qs.merge(ex,on=["chelator","metal"],how="left"); qs["dE_exch_eV"]=qs["dE_exch_eV"].fillna(0.0)
QMF=["dE_exch_eV","q_metal","dipole_D","CN_N","CN_O","CN_total",
     "d_MN_mean","d_MO_mean","q_abs_sum","q_min","q_max","rgyr"]
d=df.merge(qs[["chelator","metal"]+QMF],on=["chelator","metal"],how="inner")
cnt=d.groupby("set_id")["chelator"].nunique(); d=d[d.set_id.isin(cnt[cnt>=2].index)].copy()
d["y"]=d["logv"]-d.groupby("set_id")["logv"].transform("mean")
for c in ["d_MN_mean","d_MO_mean"]: d[c]=d[c].fillna(d[c].median())
X_qm=d[QMF].values.astype(float)
X_cls=pd.DataFrame([featurize(c,m) for c,m in zip(d.chelator,d.metal)])[FEATURES].values
fams=d.vector_family.values; sets=d.set_id.values; y=d["y"].values
def lopo(X,yy,a=1.0):
    p=np.full(len(yy),np.nan)
    for f in sorted(set(fams)):
        tr,te=fams!=f,fams==f
        if tr.sum()<10 or te.sum()<2: continue
        p[te]=make_pipeline(StandardScaler(),Ridge(alpha=a)).fit(X[tr],yy[tr]).predict(X[te])
    o=~np.isnan(p); return r2_score(yy[o],p[o])
obs={"QM_solv":lopo(X_qm,y),"classical":lopo(X_cls,y)}
null={k:[] for k in obs}
for _ in range(NPERM):
    yp=y.copy()
    for s in set(sets):
        m=sets==s; v=yp[m]; rng.shuffle(v); yp[m]=v
    null["QM_solv"].append(lopo(X_qm,yp)); null["classical"].append(lopo(X_cls,yp))
print(f"n={len(d)} 세트={d.set_id.nunique()} 계열={d.vector_family.nunique()} / {NPERM}회 순열")
res={}
for k,o in obs.items():
    n=np.array(null[k]); p=(np.sum(n>=o)+1)/(NPERM+1)
    print(f"{k:<11} 관측 R²={o:+.4f} | 귀무 평균 {n.mean():+.4f} "
          f"[{np.percentile(n,2.5):+.4f}, {np.percentile(n,97.5):+.4f}] | p={p:.4f}")
    res[k]=dict(observed=float(o),null_mean=float(n.mean()),
        null_ci=[float(np.percentile(n,2.5)),float(np.percentile(n,97.5))],p=float(p))
res["n_perm"]=NPERM; res["n_rows"]=int(len(d))
json.dump(res,open("perm_solv.json","w"),indent=1); print("→ perm_solv.json")
