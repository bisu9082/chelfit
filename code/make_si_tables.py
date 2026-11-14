# -*- coding: utf-8 -*-
"""ESI 표 3종을 CSV에서 직접 생성. 손으로 옮겨적지 않음."""
import pandas as pd,numpy as np,re,json
def esc(s):
    s="" if pd.isna(s) else str(s)
    for a,b in [("&",r"\&"),("%",r"\%"),("_",r"\_"),("#",r"\#"),("~",r"\textasciitilde{}"),
                ("–","--"),("α",r"$\alpha$"),("β",r"$\beta$"),("μ",r"$\mu$")]:
        s=s.replace(a,b)
    return s
def ce(m):   # 금속 라벨 → mhchem
    return "---" if m=="apo" else (r"[\ce{AlF}]$^{2+}$" if m=="AlF" else f"\\ce{{{m}}}")
def shortdoi(d): return esc(d.replace("https://doi.org/",""))

# ---------- Table S1 : 친화도 원데이터 106건 ----------
a=pd.read_csv("affinity_raw.csv",encoding="utf-8-sig")
a=a.sort_values(["vector_family","vector","doi","chelator","metal"]).reset_index(drop=True)
L=[r"\small", r"\begin{longtable}{@{}p{1.8cm}p{3.0cm}p{1.7cm}p{0.9cm}p{1.2cm}r@{\hskip 10pt}r@{\hskip 12pt}p{3.3cm}@{}}",
   r"\caption{Complete matched affinity dataset. Every value was transcribed by hand from the "
   r"open-access full text of the cited study; no value was parsed automatically, interpolated or "
   r"estimated. Metal `---' denotes the metal-free conjugate.}\label{tab:S1}\\",
   r"\toprule",
   r"Family & Vector & Chelator & Metal & Endpoint & Value (nM) & SD & Source DOI \\",
   r"\midrule \endfirsthead",
   r"\multicolumn{8}{@{}l}{\small\itshape Table S1 continued}\\ \toprule",
   r"Family & Vector & Chelator & Metal & Endpoint & Value (nM) & SD & Source DOI \\",
   r"\midrule \endhead",
   r"\bottomrule \endlastfoot"]
for r_ in a.itertuples():
    sd="" if pd.isna(r_.sd_nM) else f"{r_.sd_nM:g}"
    L.append(f"{esc(r_.vector_family)} & {esc(r_.vector)} & {esc(r_.chelator)} & {ce(r_.metal)} & "
             f"{esc(r_.endpoint)} & {r_.value_nM:g} & {sd} & \\tiny{{{shortdoi(r_.doi)}}} \\\\")
L.append(r"\bottomrule")
L.append(r"\end{longtable}")
L.append(r"\normalsize")
open("tabS1.tex","w").write("\n".join(L))
print(f"Table S1: {len(a)}행")

# ---------- Table S2 : 양자화학 복합체 32종 ----------
q=pd.read_csv("qm_descriptors_solv.csv",encoding="utf-8-sig")
e=pd.read_csv("exchange_energy.csv",encoding="utf-8-sig")[["chelator","metal","dE_exch_kcal"]]
q=q.merge(e,on=["chelator","metal"],how="left").sort_values(["chelator","metal"])
L=[r"\small", r"\begin{longtable}{@{}p{1.6cm}p{1.2cm}r@{\hskip 10pt}r@{\hskip 10pt}l@{\hskip 10pt}c@{\hskip 10pt}r@{\hskip 10pt}r@{\hskip 10pt}r@{}}",
   r"\caption{Optimised metal--chelate complexes (GFN2-xTB with ALPB implicit solvation in water). "
   r"`Coord.\ verified' records whether the optimised coordination mode matched the mode known for "
   r"that chelator; all 27 metal-bearing complexes passed. The ligand-exchange energy is referenced "
   r"to DOTA at the same metal, so it is zero by definition for DOTA and undefined for the "
   r"metal-free ligands.}"
   r"\label{tab:S2}\\",
   r"\toprule",
   r"Chelator & Metal & Charge & Atoms & CN & Coord.\ verified & "
   r"$\Delta\Delta E$ (kcal\,mol$^{-1}$) & $\mu$ (D) & $q_\mathrm{M}$ \\",
   r"\midrule \endfirsthead",
   r"\multicolumn{9}{@{}l}{\small\itshape Table S2 continued}\\ \toprule",
   r"Chelator & Metal & Charge & Atoms & CN & Coord.\ verified & "
   r"$\Delta\Delta E$ (kcal\,mol$^{-1}$) & $\mu$ (D) & $q_\mathrm{M}$ \\",
   r"\midrule \endhead", r"\bottomrule \endlastfoot"]
for r_ in q.itertuples():
    apo=(r_.metal=="apo")
    cn="---" if apo else f"{int(r_.CN_total)} (N$_{{{int(r_.CN_N)}}}$O$_{{{int(r_.CN_O)}}}$)"
    ok="---" if apo else (r"\checkmark" if r_.coord_ok else "no")
    dE=("---" if pd.isna(r_.dE_exch_kcal) else
        ("0 (ref.)" if r_.chelator=="DOTA" else f"{r_.dE_exch_kcal:+.1f}"))
    qm="---" if apo else f"{r_.q_metal:+.3f}"
    L.append(f"{esc(r_.chelator)} & {ce(r_.metal)} & {int(r_.charge):+d} & {int(r_.n_atoms)} & "
             f"{cn} & {ok} & {dE} & {r_.dipole_D:.1f} & {qm} \\\\")
L += [r"\bottomrule", r"\end{longtable}", r"\normalsize"]
open("tabS2.tex","w").write("\n".join(L))
print(f"Table S2: {len(q)}행 (배위검증 통과 {int(q[q.metal!='apo'].coord_ok.sum())}/{int((q.metal!='apo').sum())})")

# ---------- Table S3 : 친유성 ----------
l=pd.read_csv("lipophilicity.csv",encoding="utf-8-sig").sort_values(["doi","vector","chelator"])
L=[r"\small", r"\begin{longtable}{@{}p{3.0cm}p{1.7cm}p{1.0cm}r@{\hskip 14pt}p{3.4cm}p{3.2cm}@{}}",
   r"\caption{Reported lipophilicity values transcribed alongside the affinity data, used to test "
   r"lipophilicity as a confounder.}\label{tab:S3}\\",
   r"\toprule Vector & Chelator & Metal & log$D$ / log$P$ & Note & Source DOI \\",
   r"\midrule \endfirsthead",
   r"\multicolumn{6}{@{}l}{\small\itshape Table S3 continued}\\ \toprule",
   r"Vector & Chelator & Metal & log$D$ / log$P$ & Note & Source DOI \\",
   r"\midrule \endhead", r"\bottomrule \endlastfoot"]
for r_ in l.itertuples():
    L.append(f"{esc(r_.vector)} & {esc(r_.chelator)} & {ce(r_.metal)} & {r_.logD:g} & "
             f"\\tiny{{{esc(r_.note)}}} & \\tiny{{{shortdoi(r_.doi)}}} \\\\")
L += [r"\bottomrule", r"\end{longtable}", r"\normalsize"]
open("tabS3.tex","w").write("\n".join(L))
print(f"Table S3: {len(l)}행")
json.dump(dict(n_affinity=int(len(a)),n_complexes=int(len(q)),
               n_verified=int(q[q.metal!='apo'].coord_ok.sum()),
               n_metal_complexes=int((q.metal!='apo').sum()),n_logd=int(len(l))),
          open("si_counts.json","w"),indent=1)
