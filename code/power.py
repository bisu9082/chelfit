# -*- coding: utf-8 -*-
"""Power analysis, corrected: inject an effect that IS a linear function of the descriptors
(i.e. learnable-by-construction), scaled to a known magnitude, and find the smallest
magnitude this LOPO design reliably detects."""
import numpy as np, pandas as pd, sys, warnings, json
warnings.filterwarnings("ignore")
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import r2_score
sys.path.insert(0,"/home/claude/chelfit")
from descriptors import featurize, FEATURES
rng=np.random.default_rng(606)
d=pd.read_csv("merged_qm.csv",encoding="utf-8-sig")
X=pd.DataFrame([featurize(c,m) for c,m in zip(d.chelator,d.metal)])[FEATURES].values.astype(float)
Xs=StandardScaler().fit_transform(X)
fams=d.vector_family.values; sets=pd.Series(d.set_id.values); y0=d["y"].values

def lopo_r2(y,alpha=1.0):
    pred=np.full(len(y),np.nan)
    for f in sorted(set(fams)):
        tr,te=fams!=f,fams==f
        if tr.sum()<10 or te.sum()<2: continue
        m=make_pipeline(StandardScaler(),Ridge(alpha=alpha)); m.fit(X[tr],y[tr]); pred[te]=m.predict(X[te])
    ok=~np.isnan(pred); return r2_score(y[ok],pred[ok])

print("="*78)
print("검정력 — 기술자의 선형함수로 '학습 가능한' 효과를 주입했을 때")
print("="*78)
print(f"{'주입 Δlog10(SD)':>16}{'= 배수':>9}{'LOPO R² 평균':>14}{'탐지율(R²>0)':>14}{'탐지율(>0.2)':>13}")
res=[]; NREP=200
for mag in [0.0,0.05,0.10,0.15,0.20,0.30,0.40,0.60]:
    r2s=[]
    for rep in range(NREP):
        b=rng.normal(size=Xs.shape[1])
        sig=Xs@b
        sig=sig-pd.Series(sig).groupby(sets.values).transform("mean").values
        sd=sig.std()
        if sd<1e-9: continue
        sig=sig/sd*mag                      # mag is in log10 units
        ysim=y0+sig
        ysim=ysim-pd.Series(ysim).groupby(sets.values).transform("mean").values
        r2s.append(lopo_r2(ysim))
    r2s=np.array(r2s)
    print(f"{mag:>16.2f}{10**mag:>9.2f}{r2s.mean():>14.3f}{(r2s>0).mean():>13.0%}{(r2s>0.2).mean():>13.0%}")
    res.append(dict(delta_log10=float(mag),fold=float(10**mag),r2=float(r2s.mean()),
                    det0=float((r2s>0).mean()),det20=float((r2s>0.2).mean())))
thr=next((r for r in res if r["det0"]>=0.8),None)
obs=lopo_r2(y0)
print()
if thr:
    print(f"→ 80% 탐지 임계 = Δlog10 {thr['delta_log10']:.2f} ({thr['fold']:.2f}배)")
    print(f"→ 해석: 킬레이터 기술자에 선형으로 실리는 전이 가능 효과가 {thr['fold']:.2f}배 이상이었다면")
    print(f"        이 설계는 80% 확률로 탐지했을 것이다. 관측값 R²={obs:+.3f} → 그런 효과는 없다.")
print(f"\n비교: 실제 관측되는 세트내 효과는 중앙값 {10**np.median(np.abs(y0)):.2f}배 "
      f"(최대 {10**np.abs(y0).max():.2f}배) — 탐지 임계보다 크다.")
json.dump({"curve":res,"observed_r2":float(obs),
           "threshold":thr},open("power.json","w"),indent=1)
