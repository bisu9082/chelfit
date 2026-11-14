# -*- coding: utf-8 -*-
"""Reviewer-5 test: is the chelate affinity effect a lipophilicity confound?
All logD/logP values hand-transcribed from the same open-access tables as the affinity data."""
import numpy as np, pandas as pd, itertools, warnings
warnings.filterwarnings("ignore")
from scipy import stats
L=[]
def a(doi,vec,chel,metal,val,note=""):
    L.append(dict(doi=doi,vector=vec,chelator=chel,metal=metal,logD=val,note=note))

# 10.1186/s13550-025-01241-7  RM2 derivatives, logD7.4
S="10.1186/s13550-025-01241-7"
RM2={("DOTA","RM2"):(-2.78,-2.37,-2.51),("DOTA","RM2_Hse7"):(-2.90,-2.04,-2.25),
 ("DOTA","RM2_Bta8"):(-2.32,-1.73,-1.81),("DOTA","RM2_aMeTrp8"):(-2.26,-2.17,-2.28),
 ("NOTA","RM2"):(-1.93,-2.26,None),("NOTA","RM2_Hse7"):(-1.56,-1.95,None),
 ("NOTA","RM2_Bta8"):(-1.28,-1.59,None),("NOTA","RM2_aMeTrp8"):(-1.41,-1.79,None),
 ("NODAGA","RM2"):(-1.68,-2.14,None),("NODAGA","RM2_Hse7"):(-1.51,-2.00,None),
 ("NODAGA","RM2_Bta8"):(-1.14,-1.62,None),("NODAGA","RM2_aMeTrp8"):(-1.42,-1.94,None)}
for (ch,v),(g,c,l) in RM2.items():
    for m_,x in zip(("Ga","Cu","Lu"),(g,c,l)):
        if x is not None: a(S,v,ch,m_,x)
# 10.2967/jnumed.123.267126  PSMA-I&T, logD pH7.4
a("10.2967/jnumed.123.267126","PSMA-I&T","DOTAGA","Cu",-2.69)
a("10.2967/jnumed.123.267126","PSMA-I&T","NODAGA","Cu",-2.95)
# 10.3389/fnume.2024.1481343  TOC, log D(O/PBS pH7.4)
a("10.3389/fnume.2024.1481343","TOC","DOTA","Cu",-2.81)
a("10.3389/fnume.2024.1481343","TOC","NODAGA","Cu",-2.60)
a("10.3389/fnume.2024.1481343","TOC","DOTA","Ga",-3.18)
# 10.1186/s13550-014-0063-1  KuE ligands, log P(o/w)
a("10.1186/s13550-014-0063-1","Ahx-KuE","HBED_CC","Ga",-4.1)
a("10.1186/s13550-014-0063-1","FFK_SubKuE","DOTA","Ga",-3.1)
a("10.1186/s13550-014-0063-1","FFK_SubKuE","DOTA","Lu",-2.7)
a("10.1186/s13550-014-0063-1","FFK_SubKuE","DOTAGA","Ga",-3.6,"FFK/ffk 합산 보고")
a("10.1186/s13550-014-0063-1","FFK_SubKuE","DOTAGA","Lu",-3.9,"FFK/ffk 합산 보고")
LD=pd.DataFrame(L); LD.to_csv("lipophilicity.csv",index=False,encoding="utf-8-sig")

d=pd.read_csv("raw_results.csv",encoding="utf-8-sig")
m=d.merge(LD[["doi","vector","chelator","metal","logD"]],on=["doi","vector","chelator","metal"],how="inner")
print("="*76); print("친유성 교란 검정 (리뷰어 5 지적)"); print("="*76)
print(f"전사 logD/logP {len(LD)}건 (논문 4편) · 친화도와 결합된 행 {len(m)}건")
cnt=m.groupby("set_id")["chelator"].nunique(); m=m[m.set_id.isin(cnt[cnt>=2].index)]
print(f"킬레이터 2종 이상 보유 세트 {m.set_id.nunique()}개 · 행 {len(m)}건\n")
rows=[]
for sid,g in m.groupby("set_id"):
    for (_,x),(_,y) in itertools.combinations(list(g.iterrows()),2):
        rows.append(dict(set_id=sid,fam=x.vector_family,dlogD=x.logD-y.logD,daff=x.logv-y.logv))
P=pd.DataFrame(rows)
r1,p1=stats.pearsonr(P.dlogD,P.daff); r2,p2=stats.pearsonr(P.dlogD.abs(),P.daff.abs())
print(f"쌍 {len(P)}개")
print(f"ΔlogD  vs  Δlog친화도   : r = {r1:+.3f}  p = {p1:.3f}")
print(f"|ΔlogD| vs |Δlog친화도| : r = {r2:+.3f}  p = {p2:.3f}")
rng=np.random.default_rng(5)
bs=[stats.pearsonr(P.dlogD.values[i],P.daff.values[i])[0]
    for i in (rng.choice(len(P),len(P),True) for _ in range(10000))]
print(f"  bootstrap 95% CI [{np.percentile(bs,2.5):+.3f}, {np.percentile(bs,97.5):+.3f}]")
# 친유성 통제 후 킬레이터 효과가 남는가 — 부분상관
import statsmodels.api as sm
Xs=sm.add_constant(P[["dlogD"]]); res=sm.OLS(P.daff,Xs).fit()
print(f"\nΔlog친화도 ~ ΔlogD 회귀:  R² = {res.rsquared:.3f}, 기울기 p = {res.pvalues['dlogD']:.3f}")
print(f"잔차 SD {res.resid.std():.3f} vs 원 SD {P.daff.std():.3f} "
      f"→ 친유성이 설명하는 분산 {100*res.rsquared:.1f}%")
print(f"\n→ {'친유성은 주요 교란이 아님' if res.rsquared<0.15 else '친유성 교란 유의 — 통제 필요'}")
