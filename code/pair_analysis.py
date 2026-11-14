# -*- coding: utf-8 -*-
"""Fig4 c/d: 기술자 공간 거리 vs 효과 크기. 기체상/용매화 양쪽 전량 재계산.
쌍 순서 임의성은 반대칭 증강으로 제거(방향 예측 기회수준 = 0.500). 고정 시드."""
import pandas as pd,numpy as np,json,itertools
from scipy import stats
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score,roc_auc_score
SEED=20261114; NBOOT=10000
QMF=["dE_exch_eV","q_metal","dipole_D","CN_N","CN_O","CN_total",
     "d_MN_mean","d_MO_mean","q_abs_sum","q_min","q_max","rgyr"]
df=pd.read_csv("raw_results.csv",encoding="utf-8-sig")
ex=pd.read_csv("exchange_energy.csv",encoding="utf-8-sig")[["chelator","metal","dE_exch_eV"]]
def build(qfile):
    qs=pd.read_csv(qfile,encoding="utf-8-sig").merge(ex,on=["chelator","metal"],how="left")
    qs["dE_exch_eV"]=qs["dE_exch_eV"].fillna(0.0)
    d=df.merge(qs[["chelator","metal"]+QMF],on=["chelator","metal"],how="inner")
    cnt=d.groupby("set_id")["chelator"].nunique(); d=d[d.set_id.isin(cnt[cnt>=2].index)].copy()
    d["y"]=d["logv"]-d.groupby("set_id")["logv"].transform("mean")
    for c in ["d_MN_mean","d_MO_mean"]: d[c]=d[c].fillna(d[c].median())
    d=d.reset_index(drop=True)
    Z=StandardScaler().fit_transform(d[QMF].values.astype(float))
    rows=[]
    for sid,g in d.groupby("set_id"):
        for i,j in itertools.combinations(list(g.index),2):
            dv=Z[i]-Z[j]; dy=float(d.y[i]-d.y[j])
            rows.append(dict(set_id=sid,fam=d.vector_family[i],
                ddesc=float(np.linalg.norm(dv)),signed=dy,mag=abs(dy),
                **{f"v{k}":v for k,v in enumerate(dv)}))
    return pd.DataFrame(rows)
V=[f"v{k}" for k in range(len(QMF))]
def lopo(T,Xc,yc,auc=False):
    X=T[Xc].values.reshape(len(T),-1); y=T[yc].values; p=np.full(len(T),np.nan)
    for f in sorted(T.fam.unique()):
        tr,te=(T.fam!=f).values,(T.fam==f).values
        if tr.sum()<10 or te.sum()<2: continue
        p[te]=make_pipeline(StandardScaler(),Ridge(alpha=1.0)).fit(X[tr],y[tr]).predict(X[te])
    o=~np.isnan(p)
    if auc: return roc_auc_score((y[o]>0).astype(int),p[o]),float(((p[o]>0)==(y[o]>0)).mean())
    return r2_score(y[o],p[o])
OUT={"seed":SEED,"n_bootstrap":NBOOT,"n_descriptors":len(QMF)}
for tag,qf in [("gas","qm_descriptors.csv"),("solv","qm_descriptors_solv.csv")]:
    U=build(qf)
    A=pd.concat([U,U.assign(signed=-U.signed,**{c:-U[c] for c in V})],ignore_index=True)
    rng=np.random.default_rng(SEED)
    r,pv=stats.pearsonr(U.ddesc,U.mag); rho,prho=stats.spearmanr(U.ddesc,U.mag)
    bs=[stats.pearsonr(U.ddesc.values[i],U.mag.values[i])[0]
        for i in (rng.integers(0,len(U),len(U)) for _ in range(NBOOT))]
    ci=[float(np.percentile(bs,2.5)),float(np.percentile(bs,97.5))]
    rs,ps=stats.pearsonr(U.ddesc,U.signed)
    auc,acc=lopo(A,V,"signed",auc=True)
    o=dict(n_pairs=int(len(U)),pooled_r=float(r),pooled_p=float(pv),pooled_ci=ci,
           spearman=float(rho),spearman_p=float(prho),signed_r=float(rs),signed_p=float(ps),
           lopo_mag_dist=float(lopo(U,["ddesc"],"mag")),lopo_mag_vec=float(lopo(U,V,"mag")),
           lopo_dir_auc=float(auc),lopo_dir_acc=acc,lopo_signed=float(lopo(A,V,"signed")))
    OUT[tag]=o
    print(f"--- {tag.upper()} (n={o['n_pairs']})")
    print(f"  풀링 r={r:+.3f} p={pv:.4f} CI[{ci[0]:+.3f},{ci[1]:+.3f}] | Spearman rho={rho:+.3f} p={prho:.4f}")
    print(f"  부호있는 상관 r={rs:+.3f} p={ps:.4f}")
    print(f"  LOPO 크기(거리)={o['lopo_mag_dist']:+.3f} 크기(벡터)={o['lopo_mag_vec']:+.3f} "
          f"부호={o['lopo_signed']:+.3f} | 방향 AUC={auc:.3f} 정확도={acc:.3f}")
    if tag=="solv": U.to_csv("pairdist.csv",index=False)
json.dump(OUT,open("pair_summary.json","w"),indent=1)
print("→ pair_summary.json / pairdist.csv(용매화)")
