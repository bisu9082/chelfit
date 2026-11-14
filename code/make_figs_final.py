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
full=pd.read_csv("raw_results.csv",encoding="utf-8-sig")
mg=pd.read_csv("merged_qm.csv",encoding="utf-8-sig")
perm=json.load(open("perm_test.json")); qs=json.load(open("qm_summary.json"))

# ============================ FIG 1 : effect real & large, direction flips
fig=plt.figure(figsize=(7.1,4.6))
gs=fig.add_gridspec(2,2,height_ratios=[1,1.15],hspace=.55,wspace=.32)

ax=fig.add_subplot(gs[0,0])
rng_=full.groupby("set_id")["logv"].agg(lambda s:s.max()-s.min())
ax.hist(10**rng_,bins=np.logspace(0,1.45,13),color=G["lt"],edgecolor=G["dk"],lw=.7)
ax.set_xscale("log"); ax.axvline(10**rng_.median(),color=G["dk"],ls="--",lw=1.2)
ax.annotate(f"median {10**rng_.median():.1f}×",xy=(10**rng_.median(),6.4),
            xytext=(4.0,7.2),fontsize=7.5,arrowprops=dict(arrowstyle="-",lw=.7))
ax.set_xlabel("max/min affinity ratio within a matched set",fontsize=8)
ax.set_ylabel("matched sets"); ax.set_title("a   The effect is large",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)

ax=fig.add_subplot(gs[0,1])
lab=["measurement\nnoise","true\nsignal"]; val=[0.110,0.890]
ax.barh([0,1],val,color=[G["pl"],G["md"]],edgecolor=G["dk"],lw=.7,height=.55)
for i,v in enumerate(val): ax.text(v+.02,i,f"{v:.2f}",va="center",fontsize=8)
ax.set_yticks([0,1]); ax.set_yticklabels(lab,fontsize=7.5); ax.set_xlim(0,1.12)
ax.set_xlabel("fraction of within-set variance",fontsize=8)
ax.set_title("b   It is not noise",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)

ax=fig.add_subplot(gs[1,:])
P=[]
for sid,g in full.groupby("set_id"):
    d=dict(zip(g.chelator,g.logv))
    if "DOTA" in d and "NODAGA" in d: P.append((g.vector.iloc[0],g.metal.iloc[0],d["NODAGA"]-d["DOTA"]))
P.sort(key=lambda x:x[2])
cols=[G["dk"] if v>0 else G["lt"] for _,_,v in P]
ax.bar(range(len(P)),[v for _,_,v in P],color=cols,edgecolor=G["dk"],lw=.6,width=.72)
ax.axhline(0,color=G["dk"],lw=1)
ax.set_xticks(range(len(P)))
ax.set_xticklabels([f"{v}\n{m}" for v,m,_ in P],fontsize=5.6,rotation=90)
ax.set_ylabel(r"$\Delta$log$_{10}$IC$_{50}$"+"\nNODAGA − DOTA",fontsize=8)
ax.set_title("c   …but its direction flips from one vector to the next",loc="left",fontweight="bold",fontsize=9)
ax.text(.012,.93,"DOTA better",transform=ax.transAxes,fontsize=7,style="italic")
ax.text(.012,.08,"NODAGA better",transform=ax.transAxes,fontsize=7,style="italic")
npos=sum(1 for *_,v in P if v>0)
ax.text(.985,.60,str(npos)+"/"+str(len(P))+" DOTA better",transform=ax.transAxes,ha="right",fontsize=7.2)
for s in("top","right"): ax.spines[s].set_visible(False)
plt.savefig("Fig1.png",bbox_inches="tight"); plt.close()

# ============================ FIG 2 : nothing transfers (incl. QM) + permutation null
fig,axes=plt.subplots(1,2,figsize=(7.1,3.1),gridspec_kw={"width_ratios":[1.35,1],"wspace":.34})
ax=axes[0]
names=["no\neffect","universal\nrank","classical\ndesc.","GFN2-xTB\ndesc.","classical\n+ QM"]
keys=["B0 효과없음","B2 킬레이터범용","B3 고전기술자","B4 QM기술자","B5 고전+QM"]
r2=[qs["models"][k]["r2"] for k in keys]; ci=[qs["models"][k]["ci"] for k in keys]
xs=np.arange(5)
ax.bar(xs,r2,color=[G["pl"],G["lt"],G["md"],G["dk"],G["md"]],edgecolor=G["dk"],lw=.7,width=.66)
ax.errorbar(xs,r2,yerr=[[max(0,r2[i]-ci[i][0]) for i in range(5)],[max(0,ci[i][1]-r2[i]) for i in range(5)]],
            fmt="none",ecolor=G["dk"],capsize=2.5,lw=.9)
ax.axhline(0,color=G["dk"],lw=1); ax.axhline(.89,color=G["dk"],ls=":",lw=1.1)
ax.text(4.42,.80,"noise ceiling 0.89",ha="right",fontsize=7,style="italic")
ax.axhspan(.0,.89,color="#f6f6f6",zorder=0)
ax.set_xticks(xs); ax.set_xticklabels(names,fontsize=7.4)
ax.set_ylabel("LOPO $R^2$  (held-out vector family)",fontsize=8); ax.set_ylim(-.65,1.0)
ax.set_title("a   Nothing transfers — QM included",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)

ax=axes[1]
for i,(k,lbl) in enumerate([("QM","GFN2-xTB"),("classical","classical")]):
    p=perm[k]; off=i*1.0
    ax.barh(off,p["null_ci"][1]-p["null_ci"][0],left=p["null_ci"][0],height=.30,
            color=G["pl"],edgecolor=G["dk"],lw=.6)
    ax.plot([p["null_mean"]],[off],"o",ms=4,color=G["md"])
    ax.plot([p["observed"]],[off],"D",ms=6,color=G["dk"],zorder=5)
    ax.text(p["observed"],off+.30,f"p = {p['p']:.2f}",ha="center",fontsize=7.5)
ax.set_yticks([0,1]); ax.set_yticklabels(["GFN2-xTB","classical"],fontsize=8)
ax.set_ylim(-.55,1.75); ax.axvline(0,color=G["dk"],lw=.8,ls="--")
ax.set_xlabel("LOPO $R^2$",fontsize=8)
ax.legend(handles=[Patch(facecolor=G["pl"],edgecolor=G["dk"],label="null 95% (label shuffled)"),
                   plt.Line2D([],[],marker="D",ls="",color=G["dk"],label="observed",ms=5)],
          fontsize=6.6,loc="upper left",frameon=False,bbox_to_anchor=(0,1.02))
ax.set_title("b   Indistinguishable from chance",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)
plt.savefig("Fig2.png",bbox_inches="tight"); plt.close()

# ============================ FIG 3 : structure exists but is contextual
fig,axes=plt.subplots(1,2,figsize=(7.1,3.0),gridspec_kw={"wspace":.33})
ax=axes[0]
lab=["chelator\nmain effect\n(n=95)","+ chelator ×\nvector\n(n=95)","QM desc.\nin-sample\n(n=82)"]
v=[0.133,0.498,0.070]
ax.bar(range(3),v,color=[G["lt"],G["dk"],G["md"]],edgecolor=G["dk"],lw=.7,width=.6)
for i,x in enumerate(v): ax.text(i,x+.015,f"{x:.3f}",ha="center",fontsize=8)
ax.set_xticks(range(3)); ax.set_xticklabels(lab,fontsize=7)
ax.set_ylabel("$R^2$ of within-set variance",fontsize=8); ax.set_ylim(0,.62)
ax.set_title("a   The structure is real but contextual",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)

ax=axes[1]
lv=["new chelator","new vector","new paper","new vector family","new metal"]
vv=[-0.004,-0.081,-0.081,-0.122,-0.211]
ax.barh(range(5),vv,color=G["md"],edgecolor=G["dk"],lw=.7,height=.6)
ax.axvline(0,color=G["dk"],lw=1)
ax.set_yticks(range(5)); ax.set_yticklabels(lv,fontsize=7.5); ax.invert_yaxis()
ax.set_xlabel("cross-validated $R^2$",fontsize=8); ax.set_xlim(-.30,.10)
for i,x in enumerate(vv): ax.text(x-.010,i,f"{x:+.3f}",va="center",ha="right",fontsize=7)
ax.set_title("b   Failure at every split level",loc="left",fontweight="bold",fontsize=9)
for s in("top","right"): ax.spines[s].set_visible(False)
plt.savefig("Fig3.png",bbox_inches="tight"); plt.close()
print("Fig1-3 재생성 완료")
