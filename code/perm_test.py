# -*- coding: utf-8 -*-
"""Permutation test: is the observed LOPO R² distinguishable from a null where the
chelator label carries no transferable information? Shuffle chelator WITHIN each set."""
import numpy as np, pandas as pd, sys, warnings, json
warnings.filterwarnings("ignore")
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import r2_score
sys.path.insert(0,"/home/claude/chelfit")
from descriptors import featurize, FEATURES
rng=np.random.default_rng(31337)
d=pd.read_csv("merged_qm.csv",encoding="utf-8-sig")
QMF=["E_complexation_eV","q_metal","dipole_D","CN_N","CN_O","CN_total",
     "d_MN_mean","d_MO_mean","q_abs_sum","q_min","q_max","rgyr"]
X_qm=d[QMF].values.astype(float)
X_cls=pd.DataFrame([featurize(c,m) for c,m in zip(d.chelator,d.metal)])[FEATURES].values
fams=d.vector_family.values; sets=d.set_id.values

def lopo_r2(X,y,alpha=1.0):
    pred=np.full(len(y),np.nan)
    for f in sorted(set(fams)):
        tr,te=fams!=f,fams==f
        if tr.sum()<10 or te.sum()<2: continue
        m=make_pipeline(StandardScaler(),Ridge(alpha=alpha)); m.fit(X[tr],y[tr]); pred[te]=m.predict(X[te])
    ok=~np.isnan(pred); return r2_score(y[ok],pred[ok])

y=d["y"].values
obs={"QM":lopo_r2(X_qm,y),"classical":lopo_r2(X_cls,y)}
NPERM=2000
null={k:[] for k in obs}
for i in range(NPERM):
    yp=y.copy()
    for s in set(sets):                      # shuffle within matched set
        m=sets==s; v=yp[m]; rng.shuffle(v); yp[m]=v
    null["QM"].append(lopo_r2(X_qm,yp)); null["classical"].append(lopo_r2(X_cls,yp))

print("="*76); print(f"순열검정 (세트 내부 라벨 셔플, {NPERM}회)"); print("="*76)
res={}
for k in obs:
    n=np.array(null[k]); o=obs[k]
    p=(np.sum(n>=o)+1)/(NPERM+1)
    print(f"{k:<12} 관측 R²={o:+.3f} | 귀무분포 평균 {n.mean():+.3f} "
          f"[{np.percentile(n,2.5):+.3f}, {np.percentile(n,97.5):+.3f}] | p={p:.3f}")
    res[k]=dict(observed=float(o),null_mean=float(n.mean()),
                null_ci=[float(np.percentile(n,2.5)),float(np.percentile(n,97.5))],p=float(p))
print("\n해석: p가 크면 '관측 성능이 무정보 귀무와 구분되지 않는다'는 뜻 —")
print("      즉 기술자가 전이 가능한 정보를 담고 있다는 증거가 없다.")
json.dump(res,open("perm_test.json","w"),indent=1)
