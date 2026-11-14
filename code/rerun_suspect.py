# -*- coding: utf-8 -*-
"""Re-run conformationally failed complexes with a wider seed sweep.
Selection criterion = lowest energy AMONG structures with plausible coordination."""
import numpy as np, pandas as pd, sys, time, warnings, os
warnings.filterwarnings("ignore")
sys.path.insert(0,"/home/claude/chelfit")
from build_complex import make_complex, optimize, coordination, MET
from ase.io import write
from ase import Atoms
from tblite.ase import TBLite

SUS=[("DOTA","Lu"),("DOTAGA","Lu"),("DOTAGA","Y"),("DTPA","Cu"),("DTPA","Lu")]
EXPN={"DOTA":4,"DOTAGA":4,"NOTA":3,"NODAGA":3,"DTPA":3}
SEEDS=[42,7,2024,1,13,99,314,777,1234,555,8080,61]
q=pd.read_csv("qm_descriptors.csv",encoding="utf-8-sig")

def ion_energy(sym,z,u):
    a=Atoms(sym,positions=[[0,0,0]]); a.calc=TBLite(method="GFN2-xTB",charge=float(z),uhf=int(u),verbosity=0)
    return a.get_potential_energy()

t0=time.time(); fixed=[]
for chel,metal in SUS:
    cands=[]
    for sd in SEEDS:
        try:
            r=make_complex(chel,metal,seed=sd)
            if r is None: continue
            at,qq,u,mi=r; at,ns=optimize(at,qq,u)
            E=at.get_potential_energy()
            if not np.isfinite(E): continue
            nN,nO,ds=coordination(at,mi)
            cands.append((nN>=EXPN[chel]-1, E, at, qq, u, mi, nN, nO, ds, sd, ns))
        except Exception: pass
    if not cands: print(f"!! {metal}-{chel} 전부 실패"); continue
    good=[c for c in cands if c[0]]
    pick=min(good,key=lambda c:c[1]) if good else min(cands,key=lambda c:c[1])
    ok,E,at,qq,u,mi,nN,nO,ds,sd,ns=pick
    lat,lq,lu,_=make_complex(chel,"apo",seed=sd); lat,_=optimize(lat,lq,lu)
    Ecplx=E-lat.get_potential_energy()-ion_energy(*MET[metal])
    ch=at.get_charges(); pos=at.get_positions()
    dN=[d for d,s in ds if s=="N" and d<=2.9]; dO=[d for d,s in ds if s=="O" and d<=2.6]
    rec=dict(chelator=chel,metal=metal,n_atoms=len(at),charge=qq,uhf=u,seed=sd,
        E_eV=float(E),E_complexation_eV=float(Ecplx),q_metal=float(ch[mi]),
        dipole_D=float(np.linalg.norm(at.get_dipole_moment())*4.8032),
        CN_N=nN,CN_O=nO,CN_total=nN+nO,
        d_MN_mean=float(np.mean(dN)) if dN else np.nan,
        d_MO_mean=float(np.mean(dO)) if dO else np.nan,
        q_abs_sum=float(np.abs(ch).sum()),q_min=float(ch.min()),q_max=float(ch.max()),
        rgyr=float(np.sqrt(((pos-pos.mean(0))**2).sum(1).mean())),opt_steps=ns)
    old=q[(q.chelator==chel)&(q.metal==metal)].iloc[0]
    print(f"{metal}-{chel:<7} 이전 N{old.CN_N}O{old.CN_O} -> 신규 N{nN}O{nO} "
          f"(타당구조 {len(good)}/{len(cands)}) {'해결' if ok else '미해결'}  ({time.time()-t0:.0f}s)",flush=True)
    q.loc[(q.chelator==chel)&(q.metal==metal),list(rec.keys())]=list(rec.values())
    write(f"/home/claude/chelfit/qm/{chel}_{metal}.xyz",at); fixed.append((chel,metal,ok))
    q.to_csv("qm_descriptors.csv",index=False,encoding="utf-8-sig")
    q.to_csv("qm_descriptors.csv",index=False,encoding="utf-8-sig")
q.to_csv("qm_descriptors.csv",index=False,encoding="utf-8-sig")
print(f"\n수정 {sum(1 for _,_,o in fixed if o)}/{len(SUS)} 해결, {time.time()-t0:.0f}s")
