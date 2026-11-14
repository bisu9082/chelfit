# -*- coding: utf-8 -*-
"""Decisive test: do GFN2-xTB descriptors of the metal-chelate complex rescue transfer?"""
import sys, math, json, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.linear_model import Ridge, LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import r2_score, mean_squared_error
from scipy import stats
sys.path.insert(0,"/home/claude/chelfit")
from descriptors import featurize, FEATURES
rng=np.random.default_rng(20261113); NBOOT=10000; SEEDS=list(range(10))

df=pd.read_csv("raw_results.csv",encoding="utf-8-sig")
qm=pd.read_csv("qm_descriptors.csv",encoding="utf-8-sig")
QMF=["E_complexation_eV","q_metal","dipole_D","CN_N","CN_O","CN_total",
     "d_MN_mean","d_MO_mean","q_abs_sum","q_min","q_max","rgyr"]
d=df.merge(qm[["chelator","metal"]+QMF],on=["chelator","metal"],how="inner")
# re-center y within set AFTER the merge (sets may lose members)
cnt=d.groupby("set_id")["chelator"].nunique()
d=d[d.set_id.isin(cnt[cnt>=2].index)].copy()
d["y"]=d["logv"]-d.groupby("set_id")["logv"].transform("mean")
for c in ["d_MN_mean","d_MO_mean"]:
    d[c]=d[c].fillna(d[c].median())
print("="*80); print("QM 기술자 병합 결과"); print("="*80)
print(f"행 {len(d)} (원본 {len(df)}) / 매치드세트 {d.set_id.nunique()} / 계열 {d.vector_family.nunique()}")
print(f"탈락: {sorted(set(df.chelator)-set(d.chelator))} 킬레이터, AlF 금속")

X_cls=pd.DataFrame([featurize(c,m) for c,m in zip(d.chelator,d.metal)])[FEATURES].values
X_qm =d[QMF].values.astype(float)
X_all=np.hstack([X_cls,X_qm])
X_chel=pd.get_dummies(d["chelator"],dtype=float).values
y=d["y"].values; fams=d["vector_family"].values

def lopo(X,alpha=1.0):
    pred=np.full(len(y),np.nan)
    for f in sorted(set(fams)):
        tr,te=fams!=f,fams==f
        if tr.sum()<10 or te.sum()<2: continue
        m=make_pipeline(StandardScaler(),Ridge(alpha=alpha)); m.fit(X[tr],y[tr]); pred[te]=m.predict(X[te])
    ok=~np.isnan(pred); return pred,ok

def report(name,X,alpha=1.0):
    if X is None: pred,ok=np.zeros(len(y)),np.ones(len(y),bool)
    else: pred,ok=lopo(X,alpha)
    r2=r2_score(y[ok],pred[ok]); rmse=math.sqrt(mean_squared_error(y[ok],pred[ok]))
    idx=np.arange(ok.sum()); yo,po=y[ok],pred[ok]
    bs=[r2_score(yo[i],po[i]) for i in (rng.choice(idx,len(idx),True) for _ in range(NBOOT))]
    lo,hi=np.percentile(bs,[2.5,97.5])
    print(f"{name:<30}{r2:>9.3f}{rmse:>8.3f}   [{lo:+.3f}, {hi:+.3f}]")
    return dict(r2=float(r2),rmse=float(rmse),ci=[float(lo),float(hi)],n=int(ok.sum()))

print("\n"+"="*80); print("LOPO 교차검증 — 벡터계열 홀드아웃"); print("="*80)
print(f"{'모델':<30}{'R²':>9}{'RMSE':>8}{'   95% CI':>22}")
R={}
R["B0 효과없음"]      = report("B0 효과없음(0 예측)",None)
R["B2 킬레이터범용"]  = report("B2 킬레이터 범용순위",X_chel)
R["B3 고전기술자"]    = report("B3 고전 기술자",X_cls)
R["B4 QM기술자"]      = report("B4 QM 기술자 (GFN2-xTB)",X_qm)
R["B5 고전+QM"]       = report("B5 고전 + QM",X_all)
# alpha sweep for the QM model (robustness, not cherry-picking: report the whole sweep)
print("\nQM 모델 정규화 민감도 (alpha sweep):")
for a in [0.1,1.0,10.0,100.0]:
    p,ok=lopo(X_qm,a); print(f"   alpha={a:<7} R²={r2_score(y[ok],p[ok]):+.3f}")

print("\n"+"="*80); print("QM 기술자가 '세트 내부'에서는 신호를 갖는가 (상한 점검)"); print("="*80)
# in-sample fit with set fixed effects already removed -> how much CAN QM explain at best
m=make_pipeline(StandardScaler(),Ridge(alpha=1.0)).fit(X_qm,y)
print(f"QM 기술자 in-sample R² = {r2_score(y,m.predict(X_qm)):.3f}  (전이 아님, 적합 상한)")
m2=make_pipeline(StandardScaler(),Ridge(alpha=1.0)).fit(X_all,y)
print(f"고전+QM in-sample R²   = {r2_score(y,m2.predict(X_all)):.3f}")
inter=d.groupby(["chelator","vector_family"])["y"].transform("mean")
print(f"킬레이터x계열 상호작용  R² = {1-((d.y-inter)**2).sum()/((d.y-d.y.mean())**2).sum():.3f}  (설명 가능 구조의 크기)")

json.dump({"n_rows":int(len(d)),"n_sets":int(d.set_id.nunique()),"models":R,
           "qm_features":QMF},open("qm_summary.json","w"),indent=1,ensure_ascii=False)
d.to_csv("merged_qm.csv",index=False,encoding="utf-8-sig")
print("\n→ qm_summary.json / merged_qm.csv 기록")
