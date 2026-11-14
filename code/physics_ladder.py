# -*- coding: utf-8 -*-
"""물리 개선의 두 축을 분리: (기체상↔용매화) × (전하 인공물 기술자↔교환에너지).
in-sample / LOPO를 2x2로 보고."""
import pandas as pd,numpy as np,json,warnings
warnings.filterwarnings("ignore")
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
from scipy import stats
BASE=["q_metal","dipole_D","CN_N","CN_O","CN_total","d_MN_mean","d_MO_mean",
      "q_abs_sum","q_min","q_max","rgyr"]
df=pd.read_csv("raw_results.csv",encoding="utf-8-sig")
ex=pd.read_csv("exchange_energy.csv",encoding="utf-8-sig")[["chelator","metal","dE_exch_eV"]]
def prep(qfile,energy):
    q=pd.read_csv(qfile,encoding="utf-8-sig").merge(ex,on=["chelator","metal"],how="left")
    q["dE_exch_eV"]=q["dE_exch_eV"].fillna(0.0)
    F=[energy]+BASE
    d=df.merge(q[["chelator","metal"]+F],on=["chelator","metal"],how="inner")
    c=d.groupby("set_id")["chelator"].nunique(); d=d[d.set_id.isin(c[c>=2].index)].copy()
    d["y"]=d["logv"]-d.groupby("set_id")["logv"].transform("mean")
    for cc in ["d_MN_mean","d_MO_mean"]: d[cc]=d[cc].fillna(d[cc].median())
    return d.reset_index(drop=True),F
def run(d,F):
    X=d[F].values.astype(float); y=d["y"].values; fams=d.vector_family.values
    ins=r2_score(y,make_pipeline(StandardScaler(),Ridge(alpha=1.0)).fit(X,y).predict(X))
    p=np.full(len(y),np.nan)
    for f in sorted(set(fams)):
        tr,te=fams!=f,fams==f
        if tr.sum()<10 or te.sum()<2: continue
        p[te]=make_pipeline(StandardScaler(),Ridge(alpha=1.0)).fit(X[tr],y[tr]).predict(X[te])
    o=~np.isnan(p)
    rz=stats.pearsonr(d[F[0]],d["q"])[0] if "q" in d else float("nan")
    return float(ins),float(r2_score(y[o],p[o])),float(rz),int(len(d))
OUT={}
for tag,qf,en in [("gas_complexation","qm_descriptors.csv","E_complexation_eV"),
                  ("gas_exchange","qm_descriptors.csv","dE_exch_eV"),
                  ("solv_complexation","qm_descriptors_solv.csv","E_complexation_eV"),
                  ("solv_exchange","qm_descriptors_solv.csv","dE_exch_eV")]:
    d,F=prep(qf,en); ins,lo,rz,n=run(d,F)
    OUT[tag]=dict(insample_r2=ins,lopo_r2=lo,energy_vs_charge_r=rz,n=n)
    print(f"{tag:<19} n={n:<4} in-sample R²={ins:+.3f}  LOPO R²={lo:+.3f}  "
          f"에너지항-전하 상관 r={rz:+.3f}")
json.dump(OUT,open("physics_ladder.json","w"),indent=1); print("→ physics_ladder.json")
