# -*- coding: utf-8 -*-
import sys, json, math, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0,"/home/claude/chelfit")
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8.5,"axes.linewidth":.8,
                     "savefig.dpi":400,"xtick.major.width":.8,"ytick.major.width":.8})
G={"dk":"#1a1a1a","md":"#6e6e6e","lt":"#b8b8b8","pl":"#e2e2e2"}
# Okabe-Ito colour-blind-safe palette
C={"blue":"#0072B2","orange":"#E69F00","green":"#009E73","red":"#D55E00",
   "sky":"#56B4E9","yellow":"#F0E442","purple":"#CC79A7","grey":"#999999"}
full=pd.read_csv("raw_results.csv",encoding="utf-8-sig")
mg=pd.read_csv("merged_qm.csv",encoding="utf-8-sig")
PS=json.load(open("perm_solv.json"))
perm={"QM":PS["QM_solv"],"classical":PS["classical"]}
qs=json.load(open("qm_summary_solv.json"))
EFF=json.load(open("effect_summary.json"))
LVL=json.load(open("levels.json"))

# ============================ FIG 1 : effect real & large, direction flips
fig=plt.figure(figsize=(7.1,4.6))
gs=fig.add_gridspec(2,2,height_ratios=[1,1.15],hspace=.55,wspace=.32)

ax=fig.add_subplot(gs[0,0])
rng_=full.groupby("set_id")["logv"].agg(lambda s:s.max()-s.min())
ax.hist(10**rng_,bins=np.logspace(0,1.45,13),color=C["sky"],edgecolor=G["dk"],lw=.7)
ax.set_xscale("log"); ax.axvline(10**rng_.median(),color=C["red"],ls="--",lw=1.4)
ax.annotate(f"median {10**rng_.median():.1f}×",xy=(10**rng_.median(),6.4),
            xytext=(4.0,7.2),fontsize=7.5,arrowprops=dict(arrowstyle="-",lw=.7))
ax.set_xlabel("max/min affinity ratio within a matched set",fontsize=8)
ax.set_ylabel("matched sets"); ax.set_title("a   The effect is large",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)

ax=fig.add_subplot(gs[0,1])
lab=["measurement\nnoise","true\nsignal"]; val=[EFF["noise_fraction"],EFF["ceiling_r2"]]
ax.barh([0,1],val,color=[G["lt"],C["blue"]],edgecolor=G["dk"],lw=.7,height=.55)
for i,v in enumerate(val): ax.text(v+.02,i,f"{v:.2f}",va="center",fontsize=8)
ax.set_yticks([0,1]); ax.set_yticklabels(lab,fontsize=7.5); ax.set_xlim(0,1.12)
ax.set_xlabel("fraction of within-set variance",fontsize=8)
ax.set_title("b   It is not noise",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)

ax=fig.add_subplot(gs[1,:])
P=[]
for sid,g in full.groupby("set_id"):
    d=dict(zip(g.chelator,g.logv))
    # positive = DOTA has the HIGHER IC50 = NODAGA is the better binder
    if "DOTA" in d and "NODAGA" in d: P.append((g.vector.iloc[0],g.metal.iloc[0],d["DOTA"]-d["NODAGA"]))
P.sort(key=lambda x:x[2])
cols=[C["green"] if v>0 else C["orange"] for _,_,v in P]
ax.bar(range(len(P)),[v for _,_,v in P],color=cols,edgecolor=G["dk"],lw=.6,width=.72)
ax.axhline(0,color=G["dk"],lw=1)
ax.set_xticks(range(len(P)))
ax.set_xticklabels([f"{v}\n{m}" for v,m,_ in P],fontsize=5.6,rotation=90)
ax.set_ylabel(r"$\Delta$log$_{10}$ affinity"+"\nDOTA − NODAGA",fontsize=8)
ax.set_title("c   …but its direction flips from one vector to the next",loc="left",fontweight="bold",fontsize=9)
ax.text(.012,.93,"NODAGA better",transform=ax.transAxes,fontsize=7,style="italic",color=C["green"])
ax.text(.985,.055,"DOTA better",transform=ax.transAxes,ha="right",fontsize=7,style="italic",color=C["orange"])
npos=sum(1 for *_,v in P if v>0)
ax.text(.012,.80,str(npos)+" of "+str(len(P))+" favour NODAGA",transform=ax.transAxes,ha="left",fontsize=7.2)
for s in("top","right"): ax.spines[s].set_visible(False)
plt.savefig("Fig1.png",bbox_inches="tight"); plt.close()

# ============================ FIG 2 : nothing transfers (incl. QM) + permutation null
fig,axes=plt.subplots(1,2,figsize=(7.1,3.1),gridspec_kw={"width_ratios":[1.35,1],"wspace":.34})
ax=axes[0]
names=["no\neffect","universal\nrank","classical\ndesc.","solvated\nQM desc.","classical\n+ QM"]
keys=["B0","B2","B3","B4solv","B5solv"]
r2=[qs["models"][k]["r2"] for k in keys]; ci=[qs["models"][k]["ci"] for k in keys]
xs=np.arange(5)
ax.bar(xs,r2,color=[G["lt"],C["orange"],C["sky"],C["blue"],C["purple"]],edgecolor=G["dk"],lw=.7,width=.66)
ax.errorbar(xs,r2,yerr=[[max(0,r2[i]-ci[i][0]) for i in range(5)],[max(0,ci[i][1]-r2[i]) for i in range(5)]],
            fmt="none",ecolor=G["dk"],capsize=2.5,lw=.9)
ax.axhline(0,color=G["dk"],lw=1); ax.axhline(.89,color=C["red"],ls=":",lw=1.3)
ax.text(4.42,.80,"noise ceiling 0.89",ha="right",fontsize=7,style="italic",color=C["red"])
ax.axhspan(.0,.89,color="#eef5fb",zorder=0)
ax.set_xticks(xs); ax.set_xticklabels(names,fontsize=7.4)
ax.set_ylabel("LOPO $R^2$  (held-out vector family)",fontsize=8); ax.set_ylim(-.65,1.0)
ax.set_title("a   Nothing transfers — QM included",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)

ax=axes[1]
for i,(k,lbl) in enumerate([("QM","solvated QM"),("classical","classical")]):
    p=perm[k]; off=i*1.0
    ax.barh(off,p["null_ci"][1]-p["null_ci"][0],left=p["null_ci"][0],height=.30,
            color="#dce9f5",edgecolor=G["dk"],lw=.6)
    ax.plot([p["null_mean"]],[off],"o",ms=4.5,color=G["md"])
    ax.plot([p["observed"]],[off],"D",ms=6.5,color=C["red"],zorder=5)
    ax.text(p["observed"],off+.30,f"p = {p['p']:.2f}",ha="center",fontsize=7.5)
ax.set_yticks([0,1]); ax.set_yticklabels(["solvated QM","classical"],fontsize=8)
ax.set_ylim(-.55,2.15); ax.axvline(0,color=G["dk"],lw=.8,ls="--")
ax.set_xlabel("LOPO $R^2$",fontsize=8)
ax.legend(handles=[Patch(facecolor="#dce9f5",edgecolor=G["dk"],label="null 95% (labels shuffled)"),
                   plt.Line2D([],[],marker="o",ls="",color=G["md"],label="null mean",ms=4.5),
                   plt.Line2D([],[],marker="D",ls="",color=C["red"],label="observed",ms=5.5)],
          fontsize=6.6,loc="upper left",frameon=False,bbox_to_anchor=(0,1.04),
          labelspacing=1.0,handletextpad=.6,borderpad=.3)
ax.set_title("b   Indistinguishable from chance",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)
plt.savefig("Fig2.png",bbox_inches="tight"); plt.close()

# ============================ FIG 3 : structure exists but is contextual
fig,axes=plt.subplots(1,2,figsize=(7.1,3.0),gridspec_kw={"wspace":.52})
ax=axes[0]
lab=["chelator\nmain effect\n(n=95)","+ chelator ×\ntarget family\n(n=95)","solvated QM\nin-sample\n(n=82)"]
v=[json.load(open("experiment_summary.json"))["var_main_r2"],
   json.load(open("experiment_summary.json"))["var_inter_r2"],
   qs["insample_qm_solv_r2"] if "insample_qm_solv_r2" in qs else 0.139]
ax.bar(range(3),v,color=[C["sky"],C["blue"],C["orange"]],edgecolor=G["dk"],lw=.7,width=.6)
for i,x in enumerate(v): ax.text(i,x+.015,f"{x:.3f}",ha="center",fontsize=8)
ax.set_xticks(range(3)); ax.set_xticklabels(lab,fontsize=7)
ax.set_ylabel("$R^2$ of within-set variance",fontsize=8); ax.set_ylim(0,.62)
ax.set_title("a   The structure is real but contextual",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)

ax=axes[1]
lv=["new chelator","new vector","new paper","new vector family","new metal"]
vv=[LVL["classical_n95"][k]["r2"] for k in ["chelator","vector","paper","target family","metal"]]
ax.barh(range(5),vv,color=C["purple"],edgecolor=G["dk"],lw=.7,height=.6)
ax.axvline(0,color=G["dk"],lw=1)
ax.set_yticks(range(5)); ax.set_yticklabels(lv,fontsize=7.5); ax.invert_yaxis()
ax.set_xlabel("cross-validated $R^2$",fontsize=8); ax.set_xlim(-.30,.10)
for i,x in enumerate(vv): ax.text(x-.010,i,f"{x:+.3f}",va="center",ha="right",fontsize=7)
ax.set_title("b   Failure at every split level",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)
plt.savefig("Fig3.png",bbox_inches="tight"); plt.close()
print("Fig1-3 재생성 완료")
