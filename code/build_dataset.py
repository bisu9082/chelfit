# -*- coding: utf-8 -*-
"""chelfit dataset. EVERY value hand-transcribed from an open-access full-text table.
[C-3] No value is estimated, interpolated or invented. Provenance per row."""
import csv, math

# (doi, vector, vector_family, chelator, metal, endpoint, value_nM, sd_nM, note)
D = []
def add(doi, vec, fam, chel, metal, ep, val, sd, note=""):
    D.append(dict(doi=doi, vector=vec, vector_family=fam, chelator=chel, metal=metal,
                  endpoint=ep, value_nM=val, sd_nM=sd, note=note))

# ---------------------------------------------------------------- CXCR4 pentixafor
# 10.1186/s13550-016-0193-8  Table: "IC50 values of different metal-chelate conjugates
# consisting of cyclo(-d-Tyr-N-Me-d-Orn(spacer-[M3+]chelator)-l-Arg-l-2-Nal-Gly-)"
# Jurkat cells, [125I]CPCR4.3 radioligand. ONE assay, ONE scaffold, 8 chelators.
S = "10.1186/s13550-016-0193-8"; V = "pentixafor_scaffold"; F = "CXCR4"
for chel, vals in {
  "DOTA":    [("apo",102,17),("Ga",24.8,2.5),("Bi",22.1,7.0),("Cu",131,11),("Lu",40.9,12),("Y",40.8,27)],
  "DOTAGA":  [("apo",654,263),("Ga",380,102),("Bi",29.4,7.6),("Cu",1165,220),("Lu",188,2.8),("Y",126,89)],
  "NOTA":    [("apo",253,49),("Ga",17.8,7.7),("Cu",46.1,26),("AlF",220,57)],
  "NODAGA":  [("apo",275,106),("Ga",342,72),("Cu",343,10)],
  "NODA_Bn": [("apo",115,24),("AlF",329,49)],
  "DTPA":    [("apo",53,7.9),("Y",156,6.7),("Lu",111,77),("In",90,42),("Cu",175,66)],
  "DTPA_Bn": [("apo",200,73),("In",74.1,2.4)],
  "DFO_Bn":  [("apo",105,35),("Zr",148,27)],
}.items():
    for m,v,s in vals: add(S,V,F,chel,m,"IC50",v,s)

# 10.1186/s41181-016-0020-6 — pentixather scaffold (NEW values only; its pentixafor
# rows carry footnote 'a' = reproduced from the paper above -> EXCLUDED as duplicates)
S2="10.1186/s41181-016-0020-6"; V2="pentixather_scaffold"
add(S2,V2,F,"NOTA","apo","IC50",27.5,2.1); add(S2,V2,F,"NOTA","Cu","IC50",14.9,2.1)
add(S2,V2,F,"NODAGA","apo","IC50",64.9,6.5); add(S2,V2,F,"NODAGA","Cu","IC50",47.5,1.8)

# ---------------------------------------------------------------- GRPR
# 10.1038/s41598-019-52914-y  natCo-X-PEG2-RM26, PC-3 cells, 125I-Tyr4-BBN
S3="10.1038/s41598-019-52914-y"; V3="PEG2-RM26"; F3="GRPR"
for chel,ic,ics,kd,kds in [("NOTA",2.8,0.2,0.056,0.003),("NODAGA",4.7,0.7,0.066,0.025),
                           ("DOTA",2.0,0.2,0.033,0.012),("DOTAGA",4.4,0.4,0.167,0.022)]:
    add(S3,V3,F3,chel,"Co","IC50",ic,ics)
    add(S3,V3,F3,chel,"Co","KD",kd,kds,"pM in source, converted to nM")

# 10.1186/s13550-025-01241-7  RM2 derivatives, PC-3, [3-125I]MJ9.  chelator x peptide x metal
S4="10.1186/s13550-025-01241-7"; F4="GRPR"
RM2 = {
 ("DOTA","RM2"):          [("Ga",2.1,0.3),("Cu",3.2,0.6),("Lu",3.5,0.2)],
 ("DOTA","RM2_Hse7"):     [("Ga",7.7,0.7),("Cu",33.0,2.7),("Lu",19.7,1.6)],
 ("DOTA","RM2_Bta8"):     [("Ga",2.3,0.1),("Cu",4.2,0.6),("Lu",4.6,0.2)],
 ("DOTA","RM2_aMeTrp8"):  [("Ga",1.5,0.2),("Cu",4.0,0.5),("Lu",3.0,0.0)],
 ("NOTA","RM2"):          [("Ga",1.2,0.1),("Cu",1.7,0.3)],
 ("NOTA","RM2_Hse7"):     [("Ga",8.4,1.0),("Cu",14.2,2.7)],
 ("NOTA","RM2_Bta8"):     [("Ga",1.8,0.0),("Cu",2.8,0.2)],
 ("NOTA","RM2_aMeTrp8"):  [("Ga",5.0,0.7),("Cu",5.1,0.3)],
 ("NODAGA","RM2"):        [("Ga",1.4,0.2),("Cu",3.1,0.4)],
 ("NODAGA","RM2_Hse7"):   [("Ga",7.2,0.7),("Cu",11.8,2.1)],
 ("NODAGA","RM2_Bta8"):   [("Ga",2.9,0.6),("Cu",8.9,0.3)],
 ("NODAGA","RM2_aMeTrp8"):[("Ga",4.9,0.7),("Cu",12.6,0.4)],
}
for (chel,vec),vals in RM2.items():
    for m,v,s in vals: add(S4,vec,F4,chel,m,"IC50",v,s)

# ---------------------------------------------------------------- PSMA
# 10.1186/s13550-014-0063-1  KuE-based PSMA ligands
S5="10.1186/s13550-014-0063-1"; F5="PSMA"
add(S5,"Ahx-KuE",F5,"HBED_CC","apo","IC50",5.7,0.5); add(S5,"Ahx-KuE",F5,"HBED_CC","Ga","IC50",6.1,0.8)
for chel,vec,vals in [("DOTA","FFK_SubKuE",[("apo",13.1,2.3),("Ga",29.5,6.6),("Lu",54.7,6.1)]),
                      ("DOTAGA","FFK_SubKuE",[("apo",10.2,1.5),("Ga",12.1,3.9),("Lu",15.1,1.5)]),
                      ("DOTAGA","ffk_SubKuE",[("apo",13.9,0.4),("Ga",15.9,0.5),("Lu",13.1,2.2)])]:
    for m,v,s in vals: add(S5,vec,F5,chel,m,"IC50",v,s)

# 10.1186/s41181-020-00107-8  PSMA-617 scaffold, Ki.  rows citing prior work are flagged.
S6="10.1186/s41181-020-00107-8"
add(S6,"PSMA-617",F5,"DOTA","Ga","Ki",6.4,1.0,"cited: Benesova 2015")
add(S6,"PSMA-617",F5,"AAZTA5","Ga","Ki",8.7,0.9)
add(S6,"PSMA-617",F5,"DOTA","Sc","Ki",4.7,0.8,"cited: Eppard 2017")
add(S6,"PSMA-617",F5,"AAZTA5","Sc","Ki",30.6,11.5)
add(S6,"PSMA-617",F5,"DOTA","Lu","Ki",6.9,1.3,"cited: Benesova 2015")
add(S6,"PSMA-617",F5,"AAZTA5","Lu","Ki",26.6,11.1)

# 10.2967/jnumed.123.267126  PSMA-I&T with Cu-61
S7="10.2967/jnumed.123.267126"
add(S7,"PSMA-I&T",F5,"DOTAGA","Cu","IC50",11.2,2.3); add(S7,"PSMA-I&T",F5,"NODAGA","Cu","IC50",9.3,1.8)

# ---------------------------------------------------------------- SSTR
# 10.3389/fnume.2024.1481343  TOC, HEK-SST2 membranes
S8="10.3389/fnume.2024.1481343"; F8="SSTR"
add(S8,"TOC",F8,"DOTA","Cu","IC50",0.23,None,"95%CI 0.21-0.26")
add(S8,"TOC",F8,"NODAGA","Cu","IC50",0.34,None,"95%CI 0.30-0.38")
add(S8,"TOC",F8,"DOTA","Ga","IC50",0.18,None,"95%CI 0.16-0.20")
add(S8,"TOC",F8,"DOTA","Cu","KD",0.105,0.011); add(S8,"TOC",F8,"NODAGA","Cu","KD",0.078,0.007)

# ---------------------------------------------------------------- affibody / HER3
# 10.1038/s41598-018-36827-w  Z08698 affibody, SPR
S9="10.1038/s41598-018-36827-w"; V9="Z08698"; F9="affibody"
for chel,apo,ind in [("DOTA",0.012,0.080),("DOTAGA",0.015,0.008),("NOTA",0.040,0.160),("NODAGA",0.011,0.021)]:
    add(S9,V9,F9,chel,"apo","KD",apo,None,"pM in source -> nM; SPR 1500 RU")
    add(S9,V9,F9,chel,"In","KD",ind,None,"pM in source -> nM; LigandTracer")
# 10.3390/ijms21061972  same affibody family, Co-labeled, BxPC-3 LigandTracer
S10="10.3390/ijms21061972"
for chel,kd,sd in [("NOTA",0.1,0.1),("NODAGA",0.09,0.07),("DOTA",0.2,0.1),("DOTAGA",0.2,0.1)]:
    add(S10,"Z_HER3_affibody",F9,chel,"Co","KD",kd,sd,"wide CI; low discrimination")

with open("/home/claude/chelfit/affinity_raw.csv","w",newline="",encoding="utf-8-sig") as f:
    w=csv.DictWriter(f,fieldnames=list(D[0].keys())); w.writeheader(); w.writerows(D)

from collections import Counter,defaultdict
print(f"총 수치 행: {len(D)}")
print(f"논문: {len(set(r['doi'] for r in D))}  벡터: {len(set(r['vector'] for r in D))}  "
      f"킬레이터: {len(set(r['chelator'] for r in D))}  금속: {len(set(r['metal'] for r in D))}")
print("\n계열별:", dict(Counter(r['vector_family'] for r in D)))
print("엔드포인트:", dict(Counter(r['endpoint'] for r in D)))
print("킬레이터:", dict(Counter(r['chelator'] for r in D)))
print("금속:", dict(Counter(r['metal'] for r in D)))
# matched sets = same doi+vector+metal+endpoint, >=2 chelators
sets=defaultdict(set)
for r in D: sets[(r['doi'],r['vector'],r['metal'],r['endpoint'])].add(r['chelator'])
usable={k:v for k,v in sets.items() if len(v)>=2}
print(f"\n매치드 세트(킬레이터 2종 이상): {len(usable)}")
print(f"독립 대비 합계 sum(k-1): {sum(len(v)-1 for v in usable.values())}")
print(f"전체 쌍 합계 sum C(k,2):  {sum(len(v)*(len(v)-1)//2 for v in usable.values())}")
