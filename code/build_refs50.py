# -*- coding: utf-8 -*-
"""Build the expanded reference list. Every DOI verified against CrossRef."""
import json,time,urllib.request,urllib.parse,sys
UA={"User-Agent":"chelfit-refs/2.0 (mailto:noreply@example.org)"}
REFS=[
# --- core framing (existing)
("PriceOrvig2014","10.1039/c3cs60304k"),("Fani2011","10.2967/jnumed.111.087999"),
("Varasteh2015","10.1016/j.nucmedbio.2014.12.009"),("Poschenrieder2016a","10.1186/s13550-016-0193-8"),
("Mitran2016","10.3892/ijo.2016.3429"),("Mitran2019","10.1038/s41598-019-52914-y"),
("Poschenrieder2016b","10.1186/s13550-016-0227-2"),("Poschenrieder2017","10.1186/s41181-016-0020-6"),
("Fischer2025","10.1186/s13550-025-01241-7"),("Weineisen2014","10.1186/s13550-014-0063-1"),
("Sinnes2020","10.1186/s41181-020-00107-8"),("Basaco2024a","10.2967/jnumed.123.267126"),
("Basaco2024b","10.3389/fnume.2024.1481343"),("Rinne2019","10.1038/s41598-018-36827-w"),
("Rinne2020","10.3390/ijms21061972"),
# --- our own direct predecessors
("Kim2026ChelatorDFT","10.1021/acs.molpharmaceut.5c01955"),
("Kim2026PbSSTR2","10.1002/advs.78419"),
("Kim2026TAS2R","10.1021/acs.jcim.6c02473"),
("Kang2026Docking","10.1021/acsomega.6c08355"),
("Kang2026GasBench","10.3390/molecules31111884"),
("Kang2025QSAR","10.1016/j.comtox.2025.100363"),
# --- further matched chelator comparisons
("Poschenrieder2020","10.7150/thno.45537"),("Abouzayed2019","10.1002/ijc.32401"),
("Guenther2022","10.2967/jnumed.121.263323"),("Abouzayed2023","10.1186/s13550-023-01005-1"),
("Bodei2025","10.1186/s41181-025-00402-2"),("Stallons2022","10.2967/jnumed.122.264707"),
("Nicolas2018","10.1371/journal.pone.0195802"),("Fani2021","10.3390/ph14121265"),
("Maina2018","10.1186/s13550-018-0428-y"),("Rinne2019b","10.1016/j.ejpb.2019.05.008"),
("Oroujeni2025","10.1186/s41181-025-00347-6"),("Krasniqi2024","10.1186/s41181-024-00243-5"),
("Osl2020","10.1021/acs.molpharmaceut.9b00808"),("Greifenstein2021","10.1021/acs.bioconjchem.1c00305"),
# --- methods
("Ehlert2021ALPB","10.1021/acs.jctc.1c00471"),("Pracht2024CREST","10.1063/5.0197592"),
("Bannwarth2019","10.1021/acs.jctc.8b01176"),("Bannwarth2020","10.1002/wcms.1493"),
# --- evaluation / leakage / reproducibility
("Schuh2026Audit","10.1039/d6sc01799a"),("Tropsha2025Eval","10.1021/acs.jcim.5c00475"),
("Boldini2025UMAP","10.1186/s13321-025-01039-8"),("McGreivy2024Weak","10.1038/s42256-024-00897-5"),
("Kapoor2023REFORMS","10.48550/arXiv.2308.07832"),
# --- background: chelator & radiometal chemistry
("Spang2023GaChelators","10.3390/molecules28010203"),("Sinnes2020b","10.1186/s41181-020-00115-8"),
("Deberle2025Th227","10.1021/acs.jmedchem.4c02423"),("Hennrich2026TAT","10.1016/j.tips.2026.01.001"),
("Rudd2025ChemSci","10.1039/d4sc06969b"),
# --- background: targets
("Zhu2024SSTRanta","10.1007/s00259-024-06639-4"),("Mcneil2023Pb","10.1007/s00259-023-06494-9"),
("Zhang2024JR11","10.1007/s00259-024-06978-2"),
]
def get(doi):
    u=f"https://api.crossref.org/works/{urllib.parse.quote(doi)}"
    for _ in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=45) as r:
                return json.load(r)["message"]
        except Exception: time.sleep(2.5)
    return None
ok=[];bad=[]
for k,d in REFS:
    m=get(d); time.sleep(0.45)
    if not m: bad.append((k,d)); print(f"XX {k:<22} {d}"); continue
    jr=(m.get("short-container-title") or m.get("container-title") or [""])[0]
    yr=(m.get("issued",{}).get("date-parts") or [[None]])[0][0]
    ok.append(dict(key=k,doi=d,journal=jr,year=yr,volume=m.get("volume"),
                   page=m.get("page") or m.get("article-number"),
                   title=(m.get("title") or [""])[0],
                   authors=[{"family":a.get("family") or a.get("name",""),"given":a.get("given","")}
                            for a in (m.get("author") or [])]))
    print(f"OK {k:<22} {jr[:34]:<34} {yr}")
json.dump({"ok":ok,"bad":bad},open("refs50.json","w"),ensure_ascii=False,indent=1)
print(f"\n검증 {len(ok)}/{len(REFS)}  실패 {len(bad)}")
if bad: print("실패 목록:",bad)
