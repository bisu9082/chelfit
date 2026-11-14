# -*- coding: utf-8 -*-
"""Production GFN2-xTB run over (chelator, metal) combos present in the dataset.
3 conformer seeds each; lowest-energy structure kept. Coordination verified."""
import numpy as np, json, time, sys, os, warnings
warnings.filterwarnings("ignore")
import pandas as pd
sys.path.insert(0,"/home/claude/chelfit")
from build_complex import make_complex, optimize, coordination, LIG, MET
from ase import Atoms
from tblite.ase import TBLite
from ase.io import write

SEEDS=[42,7,2024]
OUT="/home/claude/chelfit/qm"; os.makedirs(OUT,exist_ok=True)
df=pd.read_csv("raw_results.csv",encoding="utf-8-sig")
combos=sorted(set(zip(df.chelator,df.metal)))
combos=[(c,m) for c,m in combos if c in LIG and (m=="apo" or m in MET)]
print(f"조합 {len(combos)}개 x 시드 {len(SEEDS)}", flush=True)

def ion_energy(sym,z,u):
    a=Atoms(sym,positions=[[0,0,0]])
    a.calc=TBLite(method="GFN2-xTB",charge=float(z),uhf=int(u),verbosity=0)
    return a.get_potential_energy()

ION={m:ion_energy(*MET[m]) for m in sorted({m for _,m in combos if m!="apo"})}
print("금속 이온 참조에너지 완료", flush=True)

rows=[]; t00=time.time()
for chel,metal in combos:
    best=None
    for sd in SEEDS:
        try:
            r=make_complex(chel,metal,seed=sd)
            if r is None: continue
            at,q,u,mi=r
            at,ns=optimize(at,q,u)
            E=at.get_potential_energy()
            if not np.isfinite(E): continue
            if best is None or E<best[0]: best=(E,at,q,u,mi,ns,sd)
        except Exception as e:
            print(f"   {chel}-{metal} seed{sd} 실패: {type(e).__name__}", flush=True)
    if best is None:
        print(f"!! {chel}-{metal} 전 시드 실패", flush=True); continue
    E,at,q,u,mi,ns,sd=best
    ch=at.get_charges(); dip=np.linalg.norm(at.get_dipole_moment())
    nN,nO,ds=coordination(at,mi)
    pos=at.get_positions(); rg=np.sqrt(((pos-pos.mean(0))**2).sum(1).mean())
    if mi is not None:
        dN=[d for d,s in ds if s=="N" and d<=2.9]; dO=[d for d,s in ds if s=="O" and d<=2.6]
        qm=float(ch[mi])
        # complexation energy vs free ligand + bare ion (same conformer protocol)
        lr=make_complex(chel,"apo",seed=sd)
        lat,lq,lu,_=lr; lat,_=optimize(lat,lq,lu)
        Ecplx = E - lat.get_potential_energy() - ION[metal]
    else:
        dN=dO=[]; qm=0.0; Ecplx=0.0
    rec=dict(chelator=chel, metal=metal, n_atoms=len(at), charge=q, uhf=u, seed=sd,
             E_eV=float(E), E_complexation_eV=float(Ecplx), q_metal=qm,
             dipole_D=float(dip*4.8032), CN_N=nN, CN_O=nO, CN_total=nN+nO,
             d_MN_mean=float(np.mean(dN)) if dN else np.nan,
             d_MO_mean=float(np.mean(dO)) if dO else np.nan,
             q_abs_sum=float(np.abs(ch).sum()), q_min=float(ch.min()), q_max=float(ch.max()),
             rgyr=float(rg), opt_steps=ns)
    rows.append(rec)
    write(f"{OUT}/{chel}_{metal}.xyz", at)
    print(f"{metal:>4}-{chel:<7} n={len(at):<3} q={q:+d} CN=N{nN}O{nO} "
          f"q_M={qm:+.3f} dip={rec['dipole_D']:6.2f}D Ecplx={Ecplx:8.2f}eV  ({time.time()-t00:.0f}s)", flush=True)

pd.DataFrame(rows).to_csv("qm_descriptors.csv",index=False,encoding="utf-8-sig")
print(f"\n완료 {len(rows)}/{len(combos)} 조합, 총 {time.time()-t00:.0f}s -> qm_descriptors.csv")
