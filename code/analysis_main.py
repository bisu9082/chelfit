# -*- coding: utf-8 -*-
"""chelfit Step 4 — 킬레이트 친화도 손실 분석.
[C-3] 모든 입력 수치는 affinity_raw.csv(공개전문 표 수기 전사). 가상 데이터 없음."""
import csv, math, json, sys, random
import numpy as np, pandas as pd
from collections import defaultdict
from scipy import stats
from sklearn.linear_model import Ridge, LinearRegression
from sklearn.metrics import r2_score, mean_squared_error
sys.path.insert(0, "/home/claude/chelfit")
from descriptors import featurize, FEATURES, complex_charge

SEEDS = [0,1,2,3,4,5,6,7,8,9]      # 10 seeds (기준 5 이상)
NBOOT = 10000                       # 기준 1000 이상
rng_global = np.random.default_rng(12345)

# ---------------------------------------------------------------- load
df = pd.read_csv("affinity_raw.csv", encoding="utf-8-sig")
df["logv"] = np.log10(df["value_nM"].astype(float))
df["set_id"] = df["doi"]+"|"+df["vector"]+"|"+df["metal"]+"|"+df["endpoint"]

# keep only sets with >=2 chelators (matched design)
cnt = df.groupby("set_id")["chelator"].nunique()
df = df[df["set_id"].isin(cnt[cnt>=2].index)].copy()
# within-set centering: removes assay/vector/metal baseline (= matched-pair design)
df["y"] = df["logv"] - df.groupby("set_id")["logv"].transform("mean")

print("="*78); print("데이터"); print("="*78)
print(f"행 {len(df)} / 매치드세트 {df.set_id.nunique()} / 논문 {df.doi.nunique()} / "
      f"벡터계열 {df.vector_family.nunique()} / 킬레이터 {df.chelator.nunique()}")
print(f"y(세트내 중심화 log10) 범위 {df.y.min():.2f} ~ {df.y.max():.2f}, SD {df.y.std():.3f}")

# ---------------------------------------------------------------- SANITY CHECKS
print("\n"+"="*78); print("SANITY CHECK 3항목"); print("="*78)
# (1) 단조성: 착물 음전하가 커질수록 친화도 손실이 커지는가 (사전 가설 방향)
df["q"] = [complex_charge(c,m) for c,m in zip(df.chelator, df.metal)]
rho_q, p_q = stats.spearmanr(df.q, df.y)
print(f"① 단조성  : 알짜전하 q vs y  Spearman rho={rho_q:+.3f} p={p_q:.4f}")
print(f"   사전가설(q 음수↑ → y↑ 손실↑)이면 rho는 음수여야 함 → {'부합' if rho_q<0 else '불일치'}")
# (2) 베이스라인 타당성: 세트내 중심화 후 평균은 0이어야 함
m_ok = abs(df.groupby('set_id')['y'].mean()).max()
print(f"② 베이스라인: 세트내 y 평균 최대절대값 {m_ok:.2e} (0이어야 정상) → {'통과' if m_ok<1e-9 else '실패'}")
# (3) 교차조건 일관성: IC50만 vs 전체에서 같은 패턴인가
sub = df[df.endpoint=="IC50"]
rho_sub,_ = stats.spearmanr(sub.q, sub.y)
print(f"③ 교차조건: IC50 전용 rho={rho_sub:+.3f}  vs 전체 rho={rho_q:+.3f}  "
      f"→ {'일관' if rho_q*rho_sub>0 else '불일치'}")

# ---------------------------------------------------------------- 분산분해 (헤드라인)
print("\n"+"="*78); print("분산분해 — 주효과 vs 상호작용"); print("="*78)
# 킬레이터 주효과만으로 설명되는 분산 (전체 평균 킬레이터 효과)
chel_main = df.groupby("chelator")["y"].mean()
pred_main = df["chelator"].map(chel_main)
ss_tot = ((df.y - df.y.mean())**2).sum()
ss_main = ((pred_main - df.y.mean())**2).sum()
# 킬레이터 x 벡터계열 상호작용 포함
inter = df.groupby(["chelator","vector_family"])["y"].transform("mean")
pred_int = inter.values
ss_int = ((pred_int - df.y.mean())**2).sum()
print(f"킬레이터 주효과만      설명분산 R² = {ss_main/ss_tot:.3f}")
print(f"킬레이터×벡터계열 포함 설명분산 R² = {ss_int/ss_tot:.3f}")
print(f"→ 상호작용이 추가로 설명하는 몫  = {(ss_int-ss_main)/ss_tot:.3f}")
print(f"→ 주효과/상호작용 비             = {ss_main/max(ss_int-ss_main,1e-12):.3f}")

# ---------------------------------------------------------------- 모델 비교 + LOPO
print("\n"+"="*78); print("LOPO 교차검증 (벡터계열 단위 leave-one-out)"); print("="*78)
X_desc = pd.DataFrame([featurize(c,m) for c,m in zip(df.chelator, df.metal)])[FEATURES].values
X_chel = pd.get_dummies(df["chelator"], dtype=float).values      # 범용 순위 모델
X_q    = df[["q"]].values.astype(float)
y = df["y"].values
fams = df["vector_family"].values
uf = sorted(set(fams))

def lopo(X, model_fn, seed=0):
    pred = np.full(len(y), np.nan)
    for f in uf:
        tr, te = fams!=f, fams==f
        if tr.sum()<10 or te.sum()<2: continue
        mdl = model_fn(seed); mdl.fit(X[tr], y[tr]); pred[te] = mdl.predict(X[te])
    ok = ~np.isnan(pred)
    return pred, ok

MODELS = {
 "B0 효과없음(0 예측)"   : ("zero", None),
 "B1 전하 단독"          : (X_q,    lambda s: LinearRegression()),
 "B2 킬레이터 범용순위"  : (X_chel, lambda s: LinearRegression()),
 "B3 전체기술자 Ridge"   : (X_desc, lambda s: Ridge(alpha=1.0, random_state=s)),
}
print(f"{'모델':<24}{'LOPO R²':>10}{'RMSE':>9}{'  95% CI(R²)':>22}")
results={}
for name,(X,fn) in MODELS.items():
    if isinstance(X,str) and X=="zero":
        pred = np.zeros(len(y)); ok = np.ones(len(y),bool)
    else:
        preds=[]
        for s in SEEDS:
            p,ok = lopo(X, fn, s); preds.append(p)
        pred = np.nanmean(np.vstack(preds),axis=0); ok = ~np.isnan(pred)
    r2 = r2_score(y[ok], pred[ok]); rmse = math.sqrt(mean_squared_error(y[ok], pred[ok]))
    # bootstrap CI
    idx = np.arange(ok.sum()); yo, po = y[ok], pred[ok]
    bs = [r2_score(yo[i], po[i]) for i in (rng_global.choice(idx,len(idx),True) for _ in range(NBOOT))]
    lo,hi = np.percentile(bs,[2.5,97.5])
    results[name]=dict(r2=r2,rmse=rmse,ci=[lo,hi],n=int(ok.sum()))
    print(f"{name:<24}{r2:>10.3f}{rmse:>9.3f}   [{lo:+.3f}, {hi:+.3f}]")

# ---------------------------------------------------------------- 효과크기 + 정규성
print("\n"+"="*78); print("효과크기 · 검정"); print("="*78)
best = max((k for k in results if not k.startswith("B0")), key=lambda k: results[k]["r2"])
res_b = y - np.zeros(len(y))
sw = stats.shapiro(y)
print(f"정규성(Shapiro-Wilk) on y : W={sw.statistic:.4f} p={sw.pvalue:.4f} → "
      f"{'정규' if sw.pvalue>0.05 else '비정규 → 비모수 사용'}")
# DOTA vs NODAGA 짝지은 대비 (가장 많은 쌍)
pairs=defaultdict(list)
for sid,g in df.groupby("set_id"):
    d=dict(zip(g.chelator,g.logv))
    for a in d:
        for b in d:
            if a<b: pairs[(a,b)].append(d[b]-d[a])
for (a,b),v in sorted(pairs.items(), key=lambda x:-len(x[1]))[:4]:
    v=np.array(v)
    if len(v)<4: continue
    w = stats.wilcoxon(v) if len(v)>=6 else None
    d_z = v.mean()/v.std(ddof=1) if v.std(ddof=1)>0 else float('nan')
    bs=[np.mean(rng_global.choice(v,len(v),True)) for _ in range(NBOOT)]
    lo,hi=np.percentile(bs,[2.5,97.5])
    print(f"{a} vs {b:<10} n={len(v):<3} Δlog10={v.mean():+.3f} 95%CI[{lo:+.3f},{hi:+.3f}] "
          f"Cohen's d_z={d_z:+.3f}" + (f" Wilcoxon p={w.pvalue:.3f}" if w else ""))

json.dump({"n_rows":int(len(df)),"n_sets":int(df.set_id.nunique()),
           "var_main_r2":float(ss_main/ss_tot),"var_inter_r2":float(ss_int/ss_tot),
           "models":{k:{kk:(vv if not isinstance(vv,list) else [float(x) for x in vv])
                        for kk,vv in v.items()} for k,v in results.items()},
           "sanity":{"monotonicity_rho":float(rho_q),"monotonicity_p":float(p_q),
                     "baseline_max_abs_mean":float(m_ok),"cross_condition_rho":float(rho_sub)},
           "seeds":len(SEEDS),"bootstrap":NBOOT},
          open("experiment_summary.json","w"), indent=1, ensure_ascii=False)
df.to_csv("raw_results.csv", index=False, encoding="utf-8-sig")
print("\n→ experiment_summary.json / raw_results.csv 기록")
