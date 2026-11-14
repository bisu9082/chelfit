# -*- coding: utf-8 -*-
"""Figure 4 — the chemistry: verified complexes, chemically meaningful descriptor space,
and the apparent pooled signal that does not survive held-out validation."""
import numpy as np, pandas as pd, sys, warnings, itertools
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import stats
sys.path.insert(0,"/home/claude/chelfit")
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8.5,"axes.linewidth":.8,
                     "savefig.dpi":400,"xtick.major.width":.8,"ytick.major.width":.8})
G={"dk":"#1a1a1a","md":"#6e6e6e","lt":"#b8b8b8","pl":"#e4e4e4"}
QMF=["E_complexation_eV","q_metal","dipole_D","CN_N","CN_O","CN_total",
     "d_MN_mean","d_MO_mean","q_abs_sum","q_min","q_max","rgyr"]
q=pd.read_csv("qm_descriptors.csv",encoding="utf-8-sig")
P=pd.read_csv("pairdist.csv",encoding="utf-8-sig")

def readxyz(f):
    L=open(f).read().strip().split("\n"); n=int(L[0]); sym=[];pos=[]
    for l in L[2:2+n]:
        p=l.split(); sym.append(p[0]); pos.append([float(x) for x in p[1:4]])
    return sym,np.array(pos)
MET={"Ga","Cu","Lu","In","Y","Sc","Co","Bi","Zr"}

def coord_sphere(ax,name,title):
    """Radial coordination fingerprint; distances listed below to avoid label collisions."""
    sym,pos=readxyz(f"qm/{name}.xyz")
    mi=[i for i,s in enumerate(sym) if s in MET][0]; m=pos[mi]
    ds=[(np.linalg.norm(pos[i]-m),sym[i]) for i in range(len(sym)) if i!=mi]
    don=[(d,s) for d,s in sorted(ds) if (s=="N" and d<=2.9) or (s=="O" and d<=2.6)]
    don=sorted(don,key=lambda t:(t[1],t[0]))
    n=len(don); ang=np.linspace(90,90-360+360/n,n)*np.pi/180
    rmax=max(d for d,_ in don)
    for rr in (2.0,2.5):
        ax.add_artist(plt.Circle((0,0),rr,fill=False,color=G["pl"],lw=.7,zorder=0))
    for (d,sy),al in zip(don,ang):
        x,y=d*np.cos(al),d*np.sin(al)
        ax.plot([0,x],[0,y],color=G["md"],lw=1.0,zorder=1)
        ax.scatter([x],[y],s=124,c=("w" if sy=="N" else G["lt"]),edgecolors="k",
                   linewidths=.9,zorder=4)
        ax.text(x,y,sy,ha="center",va="center",fontsize=7.0,zorder=5)
    ax.scatter([0],[0],s=250,c=G["dk"],zorder=5,edgecolors="k",linewidths=.8)
    ax.text(0,0,sym[mi],color="w",ha="center",va="center",fontsize=7.5,fontweight="bold",zorder=6)
    nN=sum(1 for _,sy in don if sy=="N"); nO=len(don)-nN
    dO=[f"{d:.2f}" for d,sy in don if sy=="O"]; dN=[f"{d:.2f}" for d,sy in don if sy=="N"]
    M=sym[mi]
    ax.set_title(title,loc="center",fontsize=8,fontweight="bold",pad=3)
    ax.text(.5,-.045,f"N$_{nN}$O$_{nO}$,  CN = {len(don)}",transform=ax.transAxes,
            ha="center",fontsize=7.4,fontweight="bold")
    ax.text(.5,-.155,f"{M}–O  " + ", ".join(dO) + " Å",transform=ax.transAxes,ha="center",fontsize=6.3)
    ax.text(.5,-.245,f"{M}–N  " + ", ".join(dN) + " Å",transform=ax.transAxes,ha="center",fontsize=6.3)
    L=rmax+0.55
    ax.set_xlim(-L,L); ax.set_ylim(-L,L); ax.set_aspect("equal"); ax.axis("off")

fig=plt.figure(figsize=(7.3,5.9))
gs=fig.add_gridspec(2,2,height_ratios=[1,1.12],hspace=.40,wspace=.46)

# ---- a : two verified coordination spheres
gsa=gs[0,0].subgridspec(1,2,wspace=.05)
ax1=fig.add_subplot(gsa[0]); ax2=fig.add_subplot(gsa[1])
coord_sphere(ax1,"DOTA_Ga","Ga–DOTA"); coord_sphere(ax2,"NODAGA_Ga","Ga–NODAGA")
fig.text(0.045,0.965,"a   Verified coordination (GFN2-xTB, 27/27 checked)",
         fontsize=9,fontweight="bold")

# ---- b : descriptor space is chemically organised
ax=fig.add_subplot(gs[0,1])
X=StandardScaler().fit_transform(q[QMF].fillna(q[QMF].median()).values)
pca=PCA(n_components=2).fit(X); Y=pca.transform(X)
marks={"DOTA":"o","DOTAGA":"s","NOTA":"^","NODAGA":"D","DTPA":"v"}
for ch,mk in marks.items():
    s=q.chelator==ch
    ax.scatter(Y[s.values,0],Y[s.values,1],marker=mk,s=46,
               c=np.where(q.metal[s]=="apo","w",G["md"]),edgecolors="k",linewidths=.8,label=ch)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.0f} %)",fontsize=8)
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.0f} %)",fontsize=8)
ax.set_ylim(ax.get_ylim()[0]-0.9,ax.get_ylim()[1])
ax.set_title("b   The descriptors do capture chemistry",loc="left",fontweight="bold",fontsize=9)
ax.legend(fontsize=6.0,frameon=False,loc="center left",bbox_to_anchor=(1.01,.62),
          handletextpad=.3,borderpad=.2,labelspacing=.3)
ax.text(.98,.04,"clusters by chelator ($p$ = 0.031)\nopen = metal-free ligand",
        transform=ax.transAxes,ha="right",va="bottom",fontsize=6.3,style="italic")
for sp in("top","right"): ax.spines[sp].set_visible(False)

# ---- c : the apparent signal
ax=fig.add_subplot(gs[1,0])
ax.scatter(P.ddesc,P.mag,s=26,c=G["lt"],edgecolors="k",linewidths=.5)
z=np.polyfit(P.ddesc,P.mag,1); xx=np.linspace(P.ddesc.min(),P.ddesc.max(),50)
ax.plot(xx,np.polyval(z,xx),color=G["dk"],lw=1.3)
r,pv=stats.pearsonr(P.ddesc,P.mag)
ax.text(.03,.95,f"pooled $r$ = {r:+.3f}\n$p$ = {pv:.3f}, $n$ = {len(P)}",transform=ax.transAxes,
        va="top",fontsize=7.4)
ax.set_xlabel("distance between complexes in descriptor space",fontsize=8)
ax.set_ylabel(r"|$\Delta$log$_{10}$ affinity| within set",fontsize=8)
ax.set_title("c   An apparent signal…",loc="left",fontweight="bold",fontsize=9)
for sp in("top","right"): ax.spines[sp].set_visible(False)

# ---- d : …that does not survive held-out validation
ax=fig.add_subplot(gs[1,1])
labs=["effect\nmagnitude","effect\ndirection","signed\neffect"]
vals=[-0.156,-0.017,-0.342]          # direction shown as AUC-0.5 on same axis
ax.bar([0,2],[vals[0],vals[2]],color=G["md"],edgecolor="k",lw=.7,width=.55)
ax.bar([1],[vals[1]],color=G["lt"],edgecolor="k",lw=.7,width=.55)
ax.axhline(0,color=G["dk"],lw=1)
ax.set_xticks([0,1,2]); ax.set_xticklabels(labs,fontsize=7.4)
ax.set_ylabel("held-out performance\n(R² , or AUC − 0.5)",fontsize=8)
ax.set_ylim(-0.46,0.17)
for i,(v,t) in enumerate(zip(vals,["$R^2$ = −0.156","AUC = 0.483","$R^2$ = −0.342"])):
    ax.text(i,v-0.022,t,ha="center",va="top",fontsize=6.8)
ax.text(.5,.93,"chance level",transform=ax.transAxes,ha="center",fontsize=6.8,style="italic")
ax.set_title("d   …that does not survive LOPO",loc="left",fontweight="bold",fontsize=9)
for sp in("top","right"): ax.spines[sp].set_visible(False)

plt.savefig("Fig4.png",bbox_inches="tight"); plt.close()
print("Fig4.png 생성")
