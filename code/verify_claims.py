# -*- coding: utf-8 -*-
import pandas as pd,numpy as np
from scipy import stats
np.random.seed(0)
df=pd.read_csv("raw_results.csv")
print("=== 효과 크기 ===")
g=df.groupby("set_id")["value_nM"]; fold=g.max()/g.min()
print(f"세트수={len(fold)} 중앙배수={fold.median():.3f} 최대={fold.max():.2f} 5배초과={int((fold>5).sum())}")
im=fold.idxmax(); print(df[df.set_id==im][["chelator","metal","value_nM"]].to_string(index=False))
print("\n=== 측정오차 천장 ===")
d=df.dropna(subset=["sd_nM"]).copy(); d["sl"]=d["sd_nM"]/(d["value_nM"]*np.log(10))
med=d["sl"].median(); con=med*np.sqrt(2); vy=df["y"].var(ddof=1); fr=con**2/vy
print(f"n(sd 보유)={len(d)} 단일오차중앙={med:.4f} 대비오차={con:.4f}")
print(f"세트내분산={vy:.4f} → 오차비율={fr*100:.1f}% 천장R²={1-fr:.3f}")
print("\n=== Fig4c/d 기술자거리 ===")
p=pd.read_csv("pairdist.csv"); print("n=",len(p))
r,pv=stats.pearsonr(p.ddesc,p.mag); rho,pr=stats.spearmanr(p.ddesc,p.mag)
print(f"|Δ| pooled Pearson r={r:+.3f} p={pv:.4f} | Spearman rho={rho:+.3f} p={pr:.4f}")
idx=[np.random.randint(0,len(p),len(p)) for _ in range(10000)]
bs=[stats.pearsonr(p.ddesc.values[i],p.mag.values[i])[0] for i in idx]
print(f"   bootstrap 95% CI [{np.percentile(bs,2.5):+.3f}, {np.percentile(bs,97.5):+.3f}]")
rs,ps=stats.pearsonr(p.ddesc,p.signed); print(f"signed  Pearson r={rs:+.3f} p={ps:.4f}")
print("\n=== Fig4d LOPO ===")
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score,roc_auc_score
def lopo(Xc):
    X=p[Xc].values.reshape(len(p),-1); y=p.mag.values; pr_=np.zeros(len(p))
    for f in p.fam.unique():
        m=(p.fam==f).values; sc=StandardScaler().fit(X[~m])
        pr_[m]=Ridge(alpha=1.0).fit(sc.transform(X[~m]),y[~m]).predict(sc.transform(X[m]))
    return r2_score(y,pr_)
print(f"거리만 LOPO R²={lopo(['ddesc']):+.3f}")
sgn=(p.signed>0).astype(int)
pr_=np.zeros(len(p))
for f in p.fam.unique():
    m=(p.fam==f).values; sc=StandardScaler().fit(p[['ddesc']].values[~m])
    pr_[m]=Ridge(alpha=1.0).fit(sc.transform(p[['ddesc']].values[~m]),p.signed.values[~m]).predict(sc.transform(p[['ddesc']].values[m]))
print(f"방향 AUC={roc_auc_score(sgn,pr_):.3f} 정확도={((pr_>0).astype(int)==sgn).mean():.3f}")
