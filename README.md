# chelfit

**A matched benchmark of chelate–vector coupling in targeted conjugates.**

Data, quantum-chemical geometries and analysis code for the article
*A metal chelate's effect on molecular recognition is not a property of the chelate.*

---

## What this is

When a metal chelate is appended to a molecule that binds a biological target, it changes how
tightly that molecule binds. Whether that change is a property of the chelate — and therefore
transferable from one vector to the next — had never been tested. We tested it.

| | |
|---|---|
| Matched affinity measurements | **106**, from 10 open-access studies |
| Matched sets | **32** (5 target families, 13 vectors, 9 chelators, 10 metal states) |
| Effect size within a set | median **2.6×**, max **25.3×** |
| Measurement noise | **11.1 %** of within-set variance (ceiling R² = 0.889) |
| Chelator main effect | R² = **0.133** |
| Chelator × target family | R² = **0.498** |
| Best transfer model (LOPO) | **all worse than predicting zero** |
| Permutation null | p = 0.88 (solvated QM), p = 0.76 (classical) |
| Rank recovery | 8–9 of 18 sets — chance |

**The effect is real and large, and it is not encoded in the chelate.**

## Layout

```
data/          affinity values, descriptor tables, result summaries
code/          the full pipeline (see "Reproduce")
geometries/    optimised Cartesian coordinates
  gas_phase/     GFN2-xTB, 32 structures
  solvated/      GFN2-xTB + ALPB(water), 32 structures
figures/       Fig. 1–4 at 400 dpi
manuscript/    LaTeX source (RSC Chemical Science format), bibliography, cover letter
```

### Key data files

| File | Contents |
|---|---|
| `data/affinity_raw.csv` | 106 affinity values, **one source DOI per row** |
| `data/raw_results.csv` | 95 rows in 32 matched sets, within-set centred |
| `data/qm_descriptors_solv.csv` | solvated descriptors, 32 complexes |
| `data/exchange_energy.csv` | ligand-exchange energies ΔΔE(L;M) referenced to DOTA |
| `data/lipophilicity.csv` | 38 logD₇.₄ / logP values transcribed from the same tables |
| `data/candidate_pool.csv` | the 42 screened candidate studies |

## Reproduce

```bash
pip install -r code/requirements.txt
cd code
python build_dataset.py        # 106 values -> affinity_raw.csv
python analysis_main.py        # matched-set analysis, sanity checks, LOPO
python run_qm.py               # GFN2-xTB, gas phase        (~20 min)
python rerun_suspect.py        # wider seed sweep for failures
python run_qm_solv.py          # GFN2-xTB + ALPB(water)     (~27 min)
python analysis_qm_solv.py     # the transfer audit
python perm_test.py            # within-set permutation null
python power.py                # injection power analysis
python rank_metric.py          # scale-free rank recovery
python lipophilicity.py        # lipophilicity confounder test
python make_figs_final.py && python make_fig4.py
```

Running this from the deposited files reproduces every number and figure in the article.

## Data integrity

Every affinity value was **transcribed by hand** from an open-access full-text table and carries its
source DOI. No value was parsed automatically, interpolated, estimated or generated. Values that a
source reproduced from earlier work were identified and excluded as duplicates — three pentixafor
studies share one underlying dataset. Coordination geometry of every complex was verified against
the known coordination mode rather than assumed; the five initial conformational-search failures are
documented in the code and were resolved by a wider seed sweep, not by relaxing the criterion.

## Citation

See `CITATION.cff`. Please cite both this deposit and the article, and credit the original
measurement authors listed in `data/affinity_raw.csv`.

## Licence

CC BY 4.0 — see `LICENSE`.
