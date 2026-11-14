# -*- coding: utf-8 -*-
"""일반화 수준별 전이 (Fig3b). 두 팔: 고전 기술자(n=95) / 고전+QM 용매화(n=82)."""
import pandas as pd,numpy as np,json,warnings,sys
warnings.filterwarnings("ignore")
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
sys.path.insert(0,"/home/claude/chelfit")
from descriptors import featurize,FEATURES
QMF=["dE_exch_eV","q_metal","dipole_D","CN_N","CN_O","CN_total",
     "d_MN_mean","d_MO_mean","q_abs_sum","q_min","q_max","rgyr"]
full=pd.read_csv("raw_results.csv",encoding="utf-8-sig")
qs=pd.read_csv("qm_descriptors_solv.csv",encoding="utf-8-sig")
ex=pd.read_csv("exchange_energy.csv",encoding="utf-8-sig")[["chelator","metal","dE_exch_eV"]]
qs=qs.merge(ex,on=["chelator","metal"],how="left"); qs["dE_exch_eV"]=qs["dE_exch_eV"].fillna(0.0)
sub=full.merge(qs[["chelator","metal"]+QMF],on=["chelator","metal"],how="inner")
cnt=sub.groupby("set_id")["chelator"].nunique(); sub=sub[sub.set_id.isin(cnt[cnt>=2].index)].copy()
sub["y"]=sub["logv"]-sub.groupby("set_id")["logv"].transform("mean")
for c in ["d_MN_mean","d_MO_mean"]: sub[c]=sub[c].fillna(sub[c].median())
LV=[("chelator","chelator"),("vector","vector"),("paper","doi"),
    ("target family","vector_family"),("metal","metal")]
def cls(d): return pd.DataFrame([featurize(c,m) for c,m in zip(d.chelator,d.metal)])[FEATURES].values
out={}
for arm,d,X in [("classical_n95",full,cls(full)),
                ("classical+QMsolv_n82",sub,np.hstack([cls(sub),sub[QMF].values.astype(float)]))]:
    y=d["y"].values; out[arm]={}
    print(f"--- {arm}")
    for lab,col in LV:
        gs=d[col].values; p=np.full(len(y),np.nan)
        for g in sorted(set(gs)):
            tr,te=gs!=g,gs==g
            if tr.sum()<10 or te.sum()<2: continue
            p[te]=make_pipeline(StandardScaler(),Ridge(alpha=1.0)).fit(X[tr],y[tr]).predict(X[te])
        o=~np.isnan(p); r2=float(r2_score(y[o],p[o]))
        out[arm][lab]=dict(r2=r2,n=int(o.sum()),n_groups=len(set(gs)))
        print(f"  held-out {lab:<14} R² = {r2:+.4f}  (n={o.sum()}, 그룹 {len(set(gs))})")
json.dump(out,open("levels.json","w"),indent=1); print("→ levels.json")
