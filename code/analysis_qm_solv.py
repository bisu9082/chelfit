# -*- coding: utf-8 -*-
"""Re-run the transfer audit with solvated (ALPB/water) descriptors and the
charge-artifact-free exchange energy. Does the conclusion change?"""
import sys, math, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.linear_model import Ridge, LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import r2_score, mean_squared_error
from scipy import stats
sys.path.insert(0,"/home/claude/chelfit")
from descriptors import featurize, FEATURES
rng=np.random.default_rng(20261114); NBOOT=10000

df=pd.read_csv("raw_results.csv",encoding="utf-8-sig")
qs=pd.read_csv("qm_descriptors_solv.csv",encoding="utf-8-sig")
ex=pd.read_csv("exchange_energy.csv",encoding="utf-8-sig")[["chelator","metal","dE_exch_eV"]]
qs=qs.merge(ex,on=["chelator","metal"],how="left")
qs["dE_exch_eV"]=qs["dE_exch_eV"].fillna(0.0)          # apo rows
# drop the charge-artifact descriptor, keep the physically sound exchange energy
QMF=["dE_exch_eV","q_metal","dipole_D","CN_N","CN_O","CN_total",
     "d_MN_mean","d_MO_mean","q_abs_sum","q_min","q_max","rgyr"]
d=df.merge(qs[["chelator","metal"]+QMF],on=["chelator","metal"],how="inner")
cnt=d.groupby("set_id")["chelator"].nunique(); d=d[d.set_id.isin(cnt[cnt>=2].index)].copy()
d["y"]=d["logv"]-d.groupby("set_id")["logv"].transform("mean")
for c in ["d_MN_mean","d_MO_mean"]: d[c]=d[c].fillna(d[c].median())

print("="*78); print("용매화 기술자 재분석"); print("="*78)
print(f"행 {len(d)} / 세트 {d.set_id.nunique()} / 계열 {d.vector_family.nunique()}")
X_cls=pd.DataFrame([featurize(c,m) for c,m in zip(d.chelator,d.metal)])[FEATURES].values
X_qm=d[QMF].values.astype(float)
X_all=np.hstack([X_cls,X_qm]); X_chel=pd.get_dummies(d["chelator"],dtype=float).values
y=d["y"].values; fams=d["vector_family"].values

def lopo(X,a=1.0):
    p=np.full(len(y),np.nan)
    for f in sorted(set(fams)):
        tr,te=fams!=f,fams==f
        if tr.sum()<10 or te.sum()<2: continue
        m=make_pipeline(StandardScaler(),Ridge(alpha=a)).fit(X[tr],y[tr]); p[te]=m.predict(X[te])
    return p,~np.isnan(p)
def rep(name,X):
    if X is None: p,ok=np.zeros(len(y)),np.ones(len(y),bool)
    else: p,ok=lopo(X)
    r2=r2_score(y[ok],p[ok]); rm=math.sqrt(mean_squared_error(y[ok],p[ok]))
    idx=np.arange(ok.sum()); yo,po=y[ok],p[ok]
    bs=[r2_score(yo[i],po[i]) for i in (rng.choice(idx,len(idx),True) for _ in range(NBOOT))]
    lo,hi=np.percentile(bs,[2.5,97.5])
    print(f"{name:<34}{r2:>9.3f}{rm:>8.3f}   [{lo:+.3f}, {hi:+.3f}]")
    return dict(r2=float(r2),rmse=float(rm),ci=[float(lo),float(hi)],n=int(ok.sum()))

print(f"\n{'모델':<34}{'LOPO R²':>9}{'RMSE':>8}{'   95% CI':>22}")
R={}
R["B0"]=rep("B0 효과없음(0 예측)",None)
R["B2"]=rep("B2 킬레이터 범용순위",X_chel)
R["B3"]=rep("B3 고전 기술자",X_cls)
R["B4solv"]=rep("B4 QM 기술자 (ALPB 용매화)",X_qm)
R["B5solv"]=rep("B5 고전 + QM(용매화)",X_all)
print("\n정규화 민감도 (QM 용매화):")
for a in [0.1,1.0,10.0,100.0]:
    p,ok=lopo(X_qm,a); print(f"   alpha={a:<7} R²={r2_score(y[ok],p[ok]):+.3f}")

print("\n"+"="*78); print("기체상 결과와 직접 비교"); print("="*78)
old=json.load(open("qm_summary.json"))["models"]
print(f"{'모델':<26}{'기체상':>10}{'용매화':>10}{'변화':>10}")
for k,ok_ in [("B2 범용순위","B2 킬레이터범용"),("B3 고전","B3 고전기술자"),
              ("B4 QM","B4 QM기술자"),("B5 고전+QM","B5 고전+QM")]:
    o=old[ok_]["r2"]; n={"B2 범용순위":R["B2"],"B3 고전":R["B3"],"B4 QM":R["B4solv"],
                         "B5 고전+QM":R["B5solv"]}[k]["r2"]
    print(f"{k:<26}{o:>10.3f}{n:>10.3f}{n-o:>+10.3f}")

print("\n"+"="*78); print("in-sample 적합 상한"); print("="*78)
m=make_pipeline(StandardScaler(),Ridge(alpha=1.0)).fit(X_qm,y)
print(f"QM(용매화) in-sample R² = {r2_score(y,m.predict(X_qm)):.3f}  (기체상 0.070)")
inter=d.groupby(["chelator","vector_family"])["y"].transform("mean")
print(f"킬레이터x계열 상호작용 R² = {1-((d.y-inter)**2).sum()/((d.y-d.y.mean())**2).sum():.3f}")
json.dump({"n_rows":int(len(d)),"n_sets":int(d.set_id.nunique()),"models":R,
           "qm_features":QMF,"solvation":"ALPB/water"},
          open("qm_summary_solv.json","w"),indent=1,ensure_ascii=False)
d.to_csv("merged_qm_solv.csv",index=False,encoding="utf-8-sig")
print("\n→ qm_summary_solv.json / merged_qm_solv.csv")
