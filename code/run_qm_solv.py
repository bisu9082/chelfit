# -*- coding: utf-8 -*-
"""GFN2-xTB / ALPB(water) re-computation of every metal-chelate complex.
Addresses: gas-phase charged species, conformer sampling, descriptor physicality.
Start points: the gas-phase optimum (fast re-relaxation) + fresh ETKDG seeds."""
import numpy as np, pandas as pd, os, sys, time, warnings, json
warnings.filterwarnings("ignore")
sys.path.insert(0,"/home/claude/chelfit")
from build_complex import make_complex, coordination, LIG, MET
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from tblite.ase import TBLite

SOLV=("alpb","water")
FRESH=[42,7,2024,1313]
OUT="/home/claude/chelfit/qm_solv"; os.makedirs(OUT,exist_ok=True)
EXPN={"DOTA":4,"DOTAGA":4,"NOTA":3,"NODAGA":3,"DTPA":3}

def calc(q,u): return TBLite(method="GFN2-xTB",charge=float(q),uhf=int(u),
                             solvation=SOLV,accuracy=1.0,verbosity=0,max_iterations=350)
def opt(at,q,u,fmax=0.08,steps=500):
    at.calc=calc(q,u); o=BFGS(at,logfile=None); o.run(fmax=fmax,steps=steps); return at

def ion_E(sym,z,u):
    a=Atoms(sym,positions=[[0,0,0]]); a.calc=calc(z,u); return a.get_potential_energy()

old=pd.read_csv("qm_descriptors.csv",encoding="utf-8-sig")
combos=[(r.chelator,r.metal) for _,r in old.iterrows()]
ION={m:ion_E(*MET[m]) for m in sorted({m for _,m in combos if m!="apo"})}
print(f"용매화 금속이온 참조에너지 완료 ({len(ION)}종)",flush=True)

# solvated free-ligand reference, best over seeds
LIGE={}
for chel in sorted(LIG):
    best=None
    for sd in FRESH:
        try:
            r=make_complex(chel,"apo",seed=sd)
            if r is None: continue
            at,q,u,_=r; at=opt(at,q,u); E=at.get_potential_energy()
            if np.isfinite(E) and (best is None or E<best[0]): best=(E,at,q,u)
        except Exception: pass
    gp=f"qm/{chel}_apo.xyz"
    if os.path.exists(gp):
        try:
            at=read(gp); q=LIG[chel][1]; at=opt(at,q,0); E=at.get_potential_energy()
            if np.isfinite(E) and (best is None or E<best[0]): best=(E,at,q,0)
        except Exception: pass
    LIGE[chel]=best[0]; write(f"{OUT}/{chel}_apo.xyz",best[1])
    print(f"  리간드 {chel:<7} E={best[0]:.3f} eV",flush=True)

rows=[]; t0=time.time()
for chel,metal in combos:
    cands=[]
    # (a) restart from gas-phase optimum
    gp=f"qm/{chel}_{metal}.xyz"
    if os.path.exists(gp):
        try:
            at=read(gp)
            q=old[(old.chelator==chel)&(old.metal==metal)].charge.iloc[0]
            u=old[(old.chelator==chel)&(old.metal==metal)].uhf.iloc[0]
            at=opt(at,q,u); E=at.get_potential_energy()
            if np.isfinite(E): cands.append((E,at,int(q),int(u),"restart"))
        except Exception as e: print(f"   {chel}-{metal} restart 실패 {type(e).__name__}",flush=True)
    # (b) fresh seeds
    for sd in FRESH:
        try:
            r=make_complex(chel,metal,seed=sd)
            if r is None: continue
            at,q,u,mi=r; at=opt(at,q,u); E=at.get_potential_energy()
            if np.isfinite(E): cands.append((E,at,q,u,f"seed{sd}"))
        except Exception: pass
    if not cands:
        print(f"!! {chel}-{metal} 전부 실패",flush=True); continue
    if metal!="apo":
        mi_of=lambda a: [i for i,s in enumerate(a.get_chemical_symbols()) if s==MET[metal][0]][0]
        good=[c for c in cands if coordination(c[1],mi_of(c[1]))[0]>=EXPN[chel]-1]
        pick=min(good,key=lambda c:c[0]) if good else min(cands,key=lambda c:c[0])
        ok=bool(good)
    else:
        pick=min(cands,key=lambda c:c[0]); ok=True
    E,at,q,u,src=pick
    if metal!="apo":
        mi=[i for i,s in enumerate(at.get_chemical_symbols()) if s==MET[metal][0]][0]
        nN,nO,ds=coordination(at,mi); qm_=float(at.get_charges()[mi])
        Ecplx=E-LIGE[chel]-ION[metal]
        dN=[d for d,s in ds if s=="N" and d<=2.9]; dO=[d for d,s in ds if s=="O" and d<=2.6]
    else:
        nN=nO=0; ds=[]; qm_=0.0; Ecplx=0.0; dN=dO=[]
    ch=at.get_charges(); pos=at.get_positions()
    rows.append(dict(chelator=chel,metal=metal,n_atoms=len(at),charge=q,uhf=u,source=src,
        n_candidates=len(cands),coord_ok=ok,
        E_eV=float(E),E_complexation_eV=float(Ecplx),q_metal=qm_,
        dipole_D=float(np.linalg.norm(at.get_dipole_moment())*4.8032),
        CN_N=nN,CN_O=nO,CN_total=nN+nO,
        d_MN_mean=float(np.mean(dN)) if dN else np.nan,
        d_MO_mean=float(np.mean(dO)) if dO else np.nan,
        q_abs_sum=float(np.abs(ch).sum()),q_min=float(ch.min()),q_max=float(ch.max()),
        rgyr=float(np.sqrt(((pos-pos.mean(0))**2).sum(1).mean()))))
    write(f"{OUT}/{chel}_{metal}.xyz",at)
    pd.DataFrame(rows).to_csv("qm_descriptors_solv.csv",index=False,encoding="utf-8-sig")
    kc=Ecplx*23.0605
    print(f"{metal:>4}-{chel:<7} cand={len(cands)} CN=N{nN}O{nO} {'OK ' if ok else 'SUS'} "
          f"Ecplx={Ecplx:7.2f} eV = {kc:8.1f} kcal/mol  ({time.time()-t0:.0f}s)",flush=True)
json.dump({"ion_E":ION,"lig_E":LIGE,"solvation":list(SOLV)},open("qm_solv_refs.json","w"),indent=1)
print(f"\n완료 {len(rows)}/{len(combos)} · {time.time()-t0:.0f}s")
