# -*- coding: utf-8 -*-
"""Scale-free second metric: can a model recover the ORDER of chelators within a held-out
matched set? Negative R² is hard to read; rank agreement is not."""
import numpy as np, pandas as pd, sys, warnings, json
warnings.filterwarnings("ignore")
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from scipy import stats
sys.path.insert(0,"/home/claude/chelfit")
from descriptors import featurize, FEATURES
rng=np.random.default_rng(99)
d=pd.read_csv("merged_qm.csv",encoding="utf-8-sig")
QMF=["E_complexation_eV","q_metal","dipole_D","CN_N","CN_O","CN_total",
     "d_MN_mean","d_MO_mean","q_abs_sum","q_min","q_max","rgyr"]
Xc=pd.DataFrame([featurize(c,m) for c,m in zip(d.chelator,d.metal)])[FEATURES].values.astype(float)
Xq=d[QMF].values.astype(float); Xch=pd.get_dummies(d.chelator,dtype=float).values
y=d["y"].values; fams=d.vector_family.values; sets=d.set_id.values

def lopo_pred(X):
    p=np.full(len(y),np.nan)
    for f in sorted(set(fams)):
        tr,te=fams!=f,fams==f
        if tr.sum()<10 or te.sum()<2: continue
        m=make_pipeline(StandardScaler(),Ridge(alpha=1.0)).fit(X[tr],y[tr]); p[te]=m.predict(X[te])
    return p

print("="*76); print("세트 내부 순위 재현 (척도 무관 지표)"); print("="*76)
print(f"{'모델':<26}{'평균 Spearman':>14}{'95% CI':>20}{'정순위 세트':>12}")
out={}
for name,X in [("킬레이터 범용순위",Xch),("고전 기술자",Xc),("QM 기술자",Xq),("고전+QM",np.hstack([Xc,Xq]))]:
    p=lopo_pred(X); rhos=[]; conc=0; tot=0
    for s in set(sets):
        m=(sets==s)&(~np.isnan(p))
        if m.sum()<3: continue
        r,_=stats.spearmanr(y[m],p[m])
        if np.isfinite(r): rhos.append(r); tot+=1; conc+= (r>0)
    rhos=np.array(rhos)
    bs=[np.mean(rng.choice(rhos,len(rhos),True)) for _ in range(10000)]
    lo,hi=np.percentile(bs,[2.5,97.5])
    print(f"{name:<26}{rhos.mean():>14.3f}   [{lo:+.3f}, {hi:+.3f}]{conc:>6}/{tot}")
    out[name]=dict(mean_rho=float(rhos.mean()),ci=[float(lo),float(hi)],
                   concordant=int(conc),n_sets=int(tot))
print("\n(세트당 킬레이터 3종 이상인 세트만. rho=0이면 순위 정보 없음, 1이면 완전 재현)")
json.dump(out,open("rank_metric.json","w"),ensure_ascii=False,indent=1)
