# -*- coding: utf-8 -*-
"""ESI 보충 그림 4종. 모든 수치는 산출물에서 직접 읽음 (리터럴 금지)."""
import json,warnings; warnings.filterwarnings("ignore")
import numpy as np,pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
import sys; sys.path.insert(0,"/home/claude/chelfit")
from descriptors import featurize,FEATURES
C={"blue":"#0072B2","orange":"#E69F00","green":"#009E73","red":"#D55E00",
   "sky":"#56B4E9","yellow":"#F0E442","purple":"#CC79A7","grey":"#999999"}
plt.rcParams.update({"font.size":8,"axes.titlesize":9,"axes.labelsize":8,
                     "xtick.labelsize":7.5,"ytick.labelsize":7.5,"axes.linewidth":.8,
                     "figure.dpi":300,"savefig.dpi":300})
def despine(ax):
    for s in ("top","right"): ax.spines[s].set_visible(False)

QMF=["dE_exch_eV","q_metal","dipole_D","CN_N","CN_O","CN_total",
     "d_MN_mean","d_MO_mean","q_abs_sum","q_min","q_max","rgyr"]
df=pd.read_csv("raw_results.csv",encoding="utf-8-sig")
ex=pd.read_csv("exchange_energy.csv",encoding="utf-8-sig")
qg=pd.read_csv("qm_descriptors.csv",encoding="utf-8-sig")
qs=pd.read_csv("qm_descriptors_solv.csv",encoding="utf-8-sig")
LAD=json.load(open("physics_ladder.json")); POW=json.load(open("power.json"))
LVL=json.load(open("levels.json"))

# ---------- Fig S1 : 전하 인공물과 그 해소 ----------
from descriptors import METAL
def prep(q):
    m=q[q.metal!="apo"].copy()
    m["z"]=[METAL[x]["z"] for x in m.metal]
    m=m.merge(ex[["chelator","metal","dE_exch_kcal"]],on=["chelator","metal"],how="left")
    m["E_kcal"]=m["E_complexation_eV"]*23.0605
    return m
G,S=prep(qg),prep(qs)
fig,axes=plt.subplots(1,2,figsize=(7.1,3.0),gridspec_kw={"wspace":.40})
ax=axes[0]
for d,mk,col,lab in [(G,"o",C["red"],"gas phase"),(S,"D",C["orange"],"solvated (ALPB/water)")]:
    r=stats.pearsonr(d.z,d.E_kcal)
    ax.scatter(d.z+(0.06 if mk=="D" else -0.06),d.E_kcal,s=38,marker=mk,c=col,
               edgecolors="k",linewidths=.55,label=f"{lab},  $r$ = {r[0]:+.3f}")
ax.set_xlabel("formal charge of the metal ion",fontsize=8)
ax.set_ylabel("complexation energy vs. bare ion\n(kcal mol$^{-1}$)",fontsize=8)
ax.set_xticks([2,3,4])
ax.set_title("a   Referenced to the bare ion, the energy is charge",
             loc="left",fontweight="bold",fontsize=8.5)
ax.legend(fontsize=6.6,frameon=False,loc="upper center",bbox_to_anchor=(.5,-.22),ncol=1); despine(ax)
ax=axes[1]
ok=S.dE_exch_kcal.notna(); r=stats.pearsonr(S.z[ok],S.dE_exch_kcal[ok])
CH={"DOTA":C["blue"],"DOTAGA":C["sky"],"NOTA":C["orange"],"NODAGA":C["green"],"DTPA":C["purple"]}
for ch,col in CH.items():
    sel=ok&(S.chelator==ch)
    if sel.sum(): ax.scatter(S.z[sel],S.dE_exch_kcal[sel],s=42,c=col,edgecolors="k",
                             linewidths=.55,label=ch)
ax.axhline(0,color=C["grey"],lw=.9,ls=":")
ax.text(.5,-.22,f"$r$ = {r[0]:+.3f},  $p$ = {r[1]:.3f},  $n$ = {int(ok.sum())}",
        transform=ax.transAxes,ha="center",va="top",fontsize=7.2)
ax.set_xlabel("formal charge of the metal ion",fontsize=8); ax.set_xticks([2,3,4])
ax.set_ylabel("ligand-exchange energy vs. DOTA\n(kcal mol$^{-1}$)",fontsize=8)
ax.set_title("b   Referenced to DOTA, the charge term cancels",
             loc="left",fontweight="bold",fontsize=8.5)
ax.legend(fontsize=6.6,frameon=False,loc="center left",bbox_to_anchor=(1.01,.5),
          handletextpad=.4,labelspacing=.9)
despine(ax)
plt.savefig("FigS1.png",bbox_inches="tight"); plt.close()
gr=stats.pearsonr(G.z,G.E_kcal); sr=stats.pearsonr(S.z,S.E_kcal)
print(f"FigS1 전하 인공물: 기체 r={gr[0]:+.3f} 용매 r={sr[0]:+.3f} 교환 r={r[0]:+.3f} p={r[1]:.3f}")
json.dump(dict(gas_complex_vs_z=float(gr[0]),solv_complex_vs_z=float(sr[0]),
               exch_vs_z=float(r[0]),exch_vs_z_p=float(r[1]),n=int(ok.sum()),
               E_gas_min=float(G.E_kcal.min()),E_gas_max=float(G.E_kcal.max()),
               exch_min=float(S.dE_exch_kcal[ok].min()),exch_max=float(S.dE_exch_kcal[ok].max())),
          open("charge_artifact.json","w"),indent=1)

# ---------- Fig S2 : 물리 사다리 ----------
order=["gas_complexation","gas_exchange","solv_complexation","solv_exchange"]
lab=["complexation\ngas","exchange\ngas","complexation\nsolvated","exchange\nsolvated"]
ins=[LAD[k]["insample_r2"] for k in order]; lo=[LAD[k]["lopo_r2"] for k in order]
fig,ax=plt.subplots(figsize=(4.6,3.1))
x=np.arange(4)
ax.plot(x,ins,"o-",color=C["blue"],lw=1.8,ms=7,mec="k",mew=.6,label="in-sample $R^2$")
ax.plot(x,lo,"s-",color=C["red"],lw=1.8,ms=7,mec="k",mew=.6,label="held-out (LOPO) $R^2$")
for xi,(a,b) in enumerate(zip(ins,lo)):
    ax.annotate(f"{a:.3f}",(xi,a),textcoords="offset points",xytext=(0,8),ha="center",fontsize=6.8)
    ax.annotate(f"{b:.3f}",(xi,b),textcoords="offset points",xytext=(0,-13),ha="center",fontsize=6.8)
ax.axhline(0,color=C["grey"],lw=.9,ls=":")
ax.set_xticks(x); ax.set_xticklabels(lab,fontsize=7)
ax.set_ylabel("$R^2$",fontsize=8); ax.set_ylim(-0.42,0.26)
ax.set_title("Every improvement in physics widens the gap",loc="left",fontweight="bold",fontsize=8.5)
ax.legend(fontsize=7,frameon=False,loc="center left")
despine(ax); plt.savefig("FigS2.png",bbox_inches="tight"); plt.close()
print("FigS2 물리 사다리")

# ---------- Fig S3 : 주입 검정력 곡선 ----------
cur=pd.DataFrame(POW["curve"])
fig,ax=plt.subplots(figsize=(4.6,3.1))
xc=[c for c in cur.columns if "fold" in c.lower()][0]
yc=[c for c in cur.columns if "det" in c.lower() or "power" in c.lower()][0]
ax.plot(cur[xc],cur[yc]*(100 if cur[yc].max()<=1.01 else 1),"o-",color=C["purple"],
        lw=1.8,ms=6,mec="k",mew=.6)
ax.axhline(80,color=C["grey"],ls="--",lw=.9); ax.text(cur[xc].min(),82,"80 %",ha="left",fontsize=7,color=C["grey"])
fold=df.groupby("set_id")["value_nM"].apply(lambda v:(v.max()/v.min()))
ax.axvline(fold.median(),color=C["orange"],ls=":",lw=1.4)
ax.text(fold.median()*1.03,8,f"observed median\n{fold.median():.2f}-fold",fontsize=6.8,color=C["orange"])
ax.set_xlabel("injected transferable effect (fold)",fontsize=8)
ax.set_ylabel("detection rate (%)",fontsize=8); ax.set_ylim(-4,104)
ax.set_title("What a transferable effect would have to be to be caught",
             loc="left",fontweight="bold",fontsize=8.5)
despine(ax); plt.savefig("FigS3.png",bbox_inches="tight"); plt.close()
print(f"FigS3 검정력 곡선 (관측 중앙 {fold.median():.2f}배)")

# ---------- Fig S4 : 수준별 전이, 두 팔 ----------
keys=["chelator","vector","paper","target family","metal"]
a=[LVL["classical_n95"][k]["r2"] for k in keys]
b=[LVL["classical+QMsolv_n82"][k]["r2"] for k in keys]
fig,ax=plt.subplots(figsize=(5.4,3.0))
y=np.arange(len(keys)); h=.36
ax.barh(y+h/2,a,h,color=C["sky"],edgecolor="k",lw=.6,label="classical descriptors ($n$ = 95)")
ax.barh(y-h/2,b,h,color=C["purple"],edgecolor="k",lw=.6,label="classical + QM ($n$ = 82)")
for yi,v in zip(y+h/2,a): ax.text(v-.02,yi,f"{v:.3f}",va="center",ha="right",fontsize=6.6)
for yi,v in zip(y-h/2,b): ax.text(v-.02,yi,f"{v:.3f}",va="center",ha="right",fontsize=6.6)
ax.axvline(0,color="k",lw=1)
ax.set_yticks(y); ax.set_yticklabels([f"held-out {k}" for k in keys],fontsize=7.5)
ax.set_xlabel("cross-validated $R^2$",fontsize=8); ax.set_xlim(-1.18,0.12)
ax.set_title("Transfer fails under every grouping",loc="left",fontweight="bold",fontsize=8.5)
ax.legend(fontsize=6.8,frameon=False,loc="upper center",bbox_to_anchor=(.5,-.20),ncol=2)
despine(ax); plt.savefig("FigS4.png",bbox_inches="tight"); plt.close()
print("FigS4 수준별 전이")
