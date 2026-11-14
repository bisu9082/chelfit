# -*- coding: utf-8 -*-
"""Metal-chelate complex descriptors. Formal values from coordination chemistry;
charges assume full deprotonation of free carboxylates at pH 7.4 (approximation - flagged)."""

# conjugated form: one arm consumed by the amide bond to the vector
# ring_N = macrocycle nitrogen count; ring_size = macrocycle atoms; free_COO = carboxylates left
CHEL = {
 # name      ring_size ring_N free_COO denticity macrocyclic
 "DOTA":     dict(ring_size=12, ring_N=4, free_COO=3, denticity=7, macrocyclic=1),
 "DOTAGA":   dict(ring_size=12, ring_N=4, free_COO=4, denticity=8, macrocyclic=1),
 "NOTA":     dict(ring_size=9,  ring_N=3, free_COO=2, denticity=5, macrocyclic=1),
 "NODAGA":   dict(ring_size=9,  ring_N=3, free_COO=3, denticity=6, macrocyclic=1),
 "NODA_Bn":  dict(ring_size=9,  ring_N=3, free_COO=2, denticity=5, macrocyclic=1),
 "DTPA":     dict(ring_size=0,  ring_N=3, free_COO=4, denticity=7, macrocyclic=0),
 "DTPA_Bn":  dict(ring_size=0,  ring_N=3, free_COO=4, denticity=7, macrocyclic=0),
 "DFO_Bn":   dict(ring_size=0,  ring_N=0, free_COO=0, denticity=6, macrocyclic=0),
 "HBED_CC":  dict(ring_size=0,  ring_N=2, free_COO=2, denticity=6, macrocyclic=0),
 "AAZTA5":   dict(ring_size=7,  ring_N=2, free_COO=4, denticity=6, macrocyclic=1),
}
# metal: formal charge, Shannon ionic radius (CN6, pm) - standard reference values
METAL = {
 "apo": dict(z=0,  r=None),
 "Ga":  dict(z=3,  r=62.0),   "Lu": dict(z=3, r=86.1),  "Y": dict(z=3, r=90.0),
 "In":  dict(z=3,  r=80.0),   "Sc": dict(z=3, r=74.5),  "Bi": dict(z=3, r=103.0),
 "Cu":  dict(z=2,  r=73.0),   "Co": dict(z=2, r=74.5),  "Zr": dict(z=4, r=72.0),
 "AlF": dict(z=2,  r=53.5),   # [AlF]2+ unit; Al3+ r=53.5
}
def complex_charge(chel, metal):
    """net formal charge of the metal-chelate unit on the conjugate"""
    return METAL[metal]["z"] - CHEL[chel]["free_COO"]

def featurize(chel, metal):
    c, m = CHEL[chel], METAL[metal]
    f = dict(
        q_complex   = complex_charge(chel, metal),
        free_COO    = c["free_COO"],
        ring_N      = c["ring_N"],
        ring_size   = c["ring_size"],
        denticity   = c["denticity"],
        macrocyclic = c["macrocyclic"],
        metal_z     = m["z"],
        metal_r     = m["r"] if m["r"] is not None else 0.0,
        has_metal   = 0 if metal == "apo" else 1,
    )
    # cavity-fit proxy: small ring prefers small ion
    f["fit_mismatch"] = 0.0 if m["r"] is None else abs(m["r"] - (62.0 if c["ring_size"]==9 else 90.0))/10.0
    return f

FEATURES = ["q_complex","free_COO","ring_N","ring_size","denticity","macrocyclic",
            "metal_z","metal_r","has_metal","fit_mismatch"]
