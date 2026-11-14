# -*- coding: utf-8 -*-
"""Metal-chelate complexes, conjugated form (one arm = N-methylamide to the vector).
Metal is threaded into the cavity via explicit DATIVE bonds before embedding,
then relaxed with GFN2-xTB. Coordination is verified, not assumed."""
import numpy as np, os, json, time
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
RDLogger.DisableLog("rdApp.*")

LIG = {
 "DOTA":   ("[O-]C(=O)CN1CCN(CC(=O)[O-])CCN(CC(=O)[O-])CCN(CC(=O)NC)CC1", -3, 4),
 "DOTAGA": ("[O-]C(=O)CN1CCN(CC(=O)[O-])CCN(CC(=O)[O-])CCN(C(CCC(=O)NC)C(=O)[O-])CC1", -4, 4),
 "NOTA":   ("[O-]C(=O)CN1CCN(CC(=O)[O-])CCN(CC(=O)NC)CC1", -2, 3),
 "NODAGA": ("[O-]C(=O)CN1CCN(CC(=O)[O-])CCN(C(CCC(=O)NC)C(=O)[O-])CC1", -3, 3),
 "DTPA":   ("[O-]C(=O)CN(CC(=O)[O-])CCN(CC(=O)NC)CCN(CC(=O)[O-])CC(=O)[O-]", -4, 3),
}
MET = {"Ga":("Ga",3,0),"In":("In",3,0),"Sc":("Sc",3,0),"Y":("Y",3,0),"Lu":("Lu",3,0),
       "Bi":("Bi",3,0),"Zr":("Zr",4,0),"Cu":("Cu",2,1),"Co":("Co",2,3)}

def donors(m):
    N=[a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=="N" and a.IsInRing()]
    if not N: N=[a.GetIdx() for a in m.GetAtoms()
                 if a.GetSymbol()=="N" and a.GetTotalNumHs()==0 and a.GetDegree()==3]
    O=[a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=="O" and a.GetFormalCharge()==-1]
    return N,O

def make_complex(chel, metal, seed=42):
    smi, lq, nN = LIG[chel]
    base = Chem.MolFromSmiles(smi)
    N,O = donors(base)
    if metal=="apo":
        m=Chem.AddHs(base); mi=None; q=lq; u=0
    else:
        msym,mz,u = MET[metal]
        rw=Chem.RWMol(base)
        ma=Chem.Atom(msym); ma.SetFormalCharge(mz); ma.SetNoImplicit(True)
        mi=rw.AddAtom(ma)
        # dative bonds: all ring N + up to 4 carboxylate O (one O per COO-)
        seenC=set(); Osel=[]
        for o in O:
            c=[nb.GetIdx() for nb in rw.GetAtomWithIdx(o).GetNeighbors()][0]
            if c in seenC: continue
            seenC.add(c); Osel.append(o)
        for d in N+Osel[:4]:
            rw.AddBond(d, mi, Chem.BondType.DATIVE)
        m=rw.GetMol(); Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_ALL ^
                                           Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        m=Chem.AddHs(m); q=lq+mz
    ps=AllChem.ETKDGv3(); ps.randomSeed=seed; ps.useRandomCoords=True
    if AllChem.EmbedMolecule(m, ps)!=0: return None
    try: AllChem.MMFFOptimizeMolecule(m, maxIters=3000)
    except Exception:
        try: AllChem.UFFOptimizeMolecule(m, maxIters=3000)
        except Exception: pass
    conf=m.GetConformer()
    pos=np.array([list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    sym=[a.GetSymbol() for a in m.GetAtoms()]
    mi2 = sym.index(MET[metal][0]) if metal!="apo" else None
    return Atoms(symbols=sym, positions=pos), q, u, mi2

def optimize(at,q,u,fmax=0.08,steps=500):
    at.calc=TBLite(method="GFN2-xTB",charge=float(q),uhf=int(u),accuracy=1.0,
                   verbosity=0,max_iterations=350)
    o=BFGS(at,logfile=None); o.run(fmax=fmax,steps=steps)
    return at,o.get_number_of_steps()

def coordination(at, mi, rmax_N=2.90, rmax_O=2.60):
    if mi is None: return 0,0,[]
    ds=[(at.get_distance(mi,i),at[i].symbol) for i in range(len(at)) if i!=mi]
    ds.sort()
    nN=sum(1 for d,s in ds if s=="N" and d<=rmax_N)
    nO=sum(1 for d,s in ds if s=="O" and d<=rmax_O)
    return nN,nO,ds[:10]

if __name__=="__main__":
    for chel,metal,exp in [("DOTA","Ga","N4O2-3"),("NOTA","Ga","N3O2-3"),("DOTA","Cu","N4O1-2")]:
        t0=time.time(); r=make_complex(chel,metal)
        if r is None: print(f"{chel}-{metal}: embed 실패"); continue
        at,q,u,mi=r; at,ns=optimize(at,q,u)
        nN,nO,ds=coordination(at,mi)
        print(f"{metal}-{chel:<7} n={len(at):<3} q={q:+d} u={u} {ns:>3}steps {time.time()-t0:>4.0f}s "
              f"| 배위 N{nN} O{nO} (기대 {exp}) | " + ", ".join(f"{s}{d:.2f}" for d,s in ds[:6]))
