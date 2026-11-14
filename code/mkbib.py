# -*- coding: utf-8 -*-
import json,re
d=json.load(open("refs50.json"))["ok"]
def esc(s):
    s=s or ""
    for a,b in [("&",r"\&"),("%",r"\%"),("_",r"\_"),("#",r"\#"),("–","--"),("—","---"),
                ("’","'"),("‘","'"),("“","``"),("”","''"),("α",r"$\alpha$"),("β",r"$\beta$"),
                ("γ",r"$\gamma$"),("μ",r"$\mu$"),("×",r"$\times$"),("≥",r"$\ge$"),("≤",r"$\le$")]:
        s=s.replace(a,b)
    s=re.sub(r"<[^>]+>","",s)
    return s.strip()
def fixcase(n):
    if n and n.isupper() and len(n)>2: return n.capitalize()
    return n
def auth(al):
    out=[]
    for a in al:
        g=fixcase(a["given"]); f=fixcase(a["family"])
        out.append(f"{g} {f}".strip() if g else f)
    if len(out)>12: out=out[:12]+["others"]
    return " and ".join(esc(x) for x in out)
JMAP={"The Journal of Chemical Physics":"J. Chem. Phys.",
 "Journal of Chemical Theory and Computation":"J. Chem. Theory Comput.",
 "Trends in Pharmacological Sciences":"Trends Pharmacol. Sci.",
 "Eur J Nucl Med Mol Imaging":"Eur. J. Nucl. Med. Mol. Imaging",
 "EJNMMI radiopharm. chem.":"EJNMMI Radiopharm. Chem.",
 "J Cheminform":"J. Cheminform.","Nat Mach Intell":"Nat. Mach. Intell.",
 "WIREs Comput Mol Sci":"WIREs Comput. Mol. Sci.",
 "European Journal of Pharmaceutics and Biopharmaceutics":"Eur. J. Pharm. Biopharm.",
 "Computational Toxicology":"Comput. Toxicol.","Adv Sci":"Adv. Sci.",
 "Advanced Science":"Adv. Sci.","ACS Omega":"ACS Omega","Molecules":"Molecules",
 "Theranostics":"Theranostics","Int J Cancer":"Int. J. Cancer",
 "International Journal of Cancer":"Int. J. Cancer","PLoS ONE":"PLoS One",
 "Pharmaceuticals":"Pharmaceuticals","Sci Rep":"Sci. Rep.",
 "Int. J. Mol. Sci.":"Int. J. Mol. Sci.","EJNMMI Res":"EJNMMI Res.",
 "Front. Nucl. Med.":"Front. Nucl. Med.","Chem. Soc. Rev.":"Chem. Soc. Rev.",
 "Journal of Nuclear Medicine":"J. Nucl. Med.","J Nucl Med":"J. Nucl. Med.",
 "Int. J. Oncol.":"Int. J. Oncol.","Nucl. Med. Biol.":"Nucl. Med. Biol.",
 "Nuclear Medicine and Biology":"Nucl. Med. Biol."}
out=[]
for r in d:
    jr=JMAP.get(r["journal"],r["journal"])
    e=[f"@article{{{r['key']},"]
    e.append(f"  author = {{{auth(r['authors'])}}},")
    e.append(f"  title = {{{esc(r['title'])}}},")
    e.append(f"  journal = {{{esc(jr)}}},")
    e.append(f"  year = {{{r['year']}}},")
    if r.get("volume"): e.append(f"  volume = {{{r['volume']}}},")
    if r.get("page"):   e.append(f"  pages = {{{str(r['page']).replace('-','--')}}},")
    e.append(f"  doi = {{{r['doi']}}},")
    e.append("}")
    out.append("\n".join(e))
open("references_new.bib","w").write("\n\n".join(out)+"\n")
print(f"{len(out)} entries -> references_new.bib")
print("\n".join(r["key"] for r in d))
