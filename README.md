# Oral fluoroquinolone step-down for non-HACEK Gram-negative infective endocarditis — simulation code, version 1.1

**No new clinical data were generated.** All exposure inputs come from published
pharmacokinetic sources (FDA product labels and published population PK models);
all MIC inputs come from the public EUCAST MIC distribution database. The output
is a probability-of-target-attainment (PTA) curve and a PK/PD eligibility verdict
for each regimen, MIC value and PK/PD target.

Corresponding author: Huiwen Chen, Department of Cardiac Surgery, Wuhan Asia Heart
Hospital. ORCID 0009-0006-4948-0174.

---

## What changed in v1.1

1. **STEP 4 correction** (already in v1.0.1): the SOME-scenario column of the
   claim-range calculation used `min(ps) >= 0.90` by mistake. It now uses
   `max(ps) >= 0.90`, while EVERY uses `min(ps) >= 0.90`. Five cells of Table 2
   changed.
2. **STEP 2 correction (new in v1.1)**: a scenario in which no MIC on the grid
   attains 90% was previously dropped from `min(got)`. It is now retained as
   `null`, so `BREAKPOINT_PTA90` agrees exactly with `CLAIM_MIC_RANGE` in all
   16 regimen × target cells. No exposure model, parameter, target, MIC grid or
   seed was changed.
3. **Three new executed scripts** close the gaps identified in the v1.0 release:
   external validation (Table 1), dose-specific renal stratification (Table 6
   and Figure 3), and sensitivity analyses (Table 3 and Section 3.6).
4. **Table 6 and Figure 3 recomputed.** The superseded table applied the label's
   500 mg renal AUC values to the 750 mg rows. v1.1 uses a dose-specific
   baseline: the 500 mg renal table is scaled by `90.7 / 47.9 = 1.8935` for the
   750 mg dose (the label gives no separate 750 mg renal table), with the
   label's dose-specific coefficient of variation and Monte Carlo PTA >= 90%.
   `Supplemental_Results_renal.json` contains a `legacy_reconstruction` block
   that reproduces the superseded values for audit.
5. **Section 3.6 rewritten.** Two statements in the earlier manuscript were not
   supported by the code and have been replaced with the measured results:
   - label-AUC CV 0.15/0.25/0.40 leaves the ciprofloxacin 500 mg fAUC >= 72
     breakpoint unchanged (EVERY 0.0625, SOME 0.125 mg/L);
   - the 80% and 95% PTA criteria shift 8 of 16 breakpoint cells by one
     doubling dilution and change 11 of 58 decision-cell eligibility classes;
     the six principal calls are unchanged at 95%, and five of six at 80%.

---

## Files

| File | Contents |
|---|---|
| `Supplemental_Code_gnb_pkpd_v3_1.py` | Primary Monte Carlo simulation (200,000 subjects/scenario). STEP 1 PTA vs MIC, STEP 2 breakpoints, STEP 3 boundary cells, STEP 4 claim ranges. |
| `Supplemental_Results_gnb_pkpd_v3_1.json` | Output of the above. |
| `Supplemental_Code_validate.py` | Seven external-validation benchmarks (Table 1). |
| `Supplemental_Results_validation.json` | Output of the above; 6 of 7 benchmarks met. |
| `Supplemental_Code_renal.py` | Dose-specific renal stratification (Table 6, Figure 3). `--plot` also regenerates the figure. |
| `Supplemental_Results_renal.json` | Output of the above, with `legacy_reconstruction`. |
| `Supplemental_Code_sensitivity.py` | S1 unbound fraction, S2 label CV, S3 PTA threshold, S4 candidate AUC variants (Table 3 and Section 3.6). |
| `Supplemental_Results_sensitivity.json` | Output of the above. |
| `Figure_3_v2.png` | Regenerated Figure 3, plotted from `Supplemental_Results_renal.json`. |
| `gnb_pkpd_v3.py`, `gnb_pkpd_v3_results.json`, `gnlog_v3.txt` | Archived v1.0 release, retained for comparison. Do not use for the current manuscript. |

## Reproducing

Requires Python 3.13 and NumPy. `matplotlib` is needed only for `--plot`.
All scripts use the fixed seed `20261008` and 200,000 virtual subjects per
scenario, so the JSON outputs are byte-for-byte reproducible.

```bash
pip install numpy matplotlib
python Supplemental_Code_gnb_pkpd_v3_1.py   # -> Supplemental_Results_gnb_pkpd_v3_1.json
python Supplemental_Code_validate.py        # -> Supplemental_Results_validation.json
python Supplemental_Code_renal.py --plot    # -> Supplemental_Results_renal.json, Figure_3_v2.png
python Supplemental_Code_sensitivity.py     # -> Supplemental_Results_sensitivity.json
```

Output map:

| Output block | Manuscript location |
|---|---|
| `PTA_VS_MIC` | Tables 4, 5; Figures 1, 2, 4 |
| `BREAKPOINT_PTA90` | Cross-check for Table 2 |
| `BOUNDARY_CELLS` | Tables 4, 5 |
| `CLAIM_MIC_RANGE` | Table 2 |
| `AUC_BENCHMARKS`, `IV_CIPRO_PTA`, `PUBLISHED_BREAKPOINTS`, `INTERNAL_CONSISTENCY` | Table 1 |
| `TABLE6`, `LEGACY_RECONSTRUCTION` | Table 6, Figure 3 |
| `S1_PROTEIN_BINDING`, `S2_LABEL_CV`, `S3_THRESHOLD_*`, `S4_CANDIDATE_AUC` | Table 3, Section 3.6 |

Note on the two PTA tables: `PTA_VS_MIC` and `BOUNDARY_CELLS` are computed from
independently drawn exposure vectors with the same fixed seed, so individual
percentages can differ by 0.1–0.3 percentage points. Tables 4 and 5 use
`BOUNDARY_CELLS`; Figures 1, 2 and 4 use `PTA_VS_MIC`.

## Inputs

| Input | Value | Source |
|---|---|---|
| Ciprofloxacin label AUC, 500 / 750 mg q12h | 27.4 / 31.6 µg·h/mL | FDA label |
| Levofloxacin label AUC, 500 / 750 mg q24h | 47.5 ± 6.7 / 90.7 ± 17.6 µg·h/mL | FDA label |
| Ciprofloxacin popPK | CL 25.4 L/h, log-scale IIV 0.678, F 0.70 | published ICU model |
| Levofloxacin popPK | CL/F 8.32 L/h, log-scale IIV 0.226, F 0.99 | published model |
| Unbound fraction, cipro / levo | 0.60/0.70/0.80, 0.62/0.69/0.76 | label protein binding |
| Levofloxacin renal AUC, 500 mg (CL_CR >80 / 50–80 / 20–49 / <20) | 47.9 / 95.6 / 182.1 / 263.5 µg·h/mL | FDA label |
| Targets | fAUC24/MIC >= 72; AUC24/MIC >= 125, 250, 500 | Delphi consensus and literature |
| MIC grid | 2^-5 … 2^3 mg/L (0.03125 – 8) | two-fold dilution series |

ECOFF values (EUCAST MIC distribution database, accessed 12 October 2026):
*E. coli* ciprofloxacin 0.06, levofloxacin 0.125; *K. pneumoniae* ciprofloxacin
0.125, levofloxacin 0.25; *P. aeruginosa* ciprofloxacin 0.5, levofloxacin
2.0 mg/L. Clinical breakpoints are from EUCAST v16.1 and CLSI M100 Ed36; those
documents are not redistributed here.

## Known limitations of the code

- No pharmacodynamic / outcome model is included; the output is eligibility,
  not efficacy.
- The 750 mg renal AUC values are derived by dose-proportional scaling of the
  500 mg label renal table because no 750 mg-specific renal table is published.
- Figures are generated from the JSON outputs; the archived `make_figures.py`
  from v1.0 still draws the superseded Figure 3. Use
  `Supplemental_Code_renal.py --plot` for the current Figure 3.

## Citation

> Chen H, Li Y, Hua B, Sun Y. *Oral fluoroquinolone step-down for non-HACEK
> Gram-negative infective endocarditis: a Monte Carlo PK/PD eligibility
> analysis.* Simulation code and result files (version 1.1). Science Data
> Bank; 2026. https://doi.org/10.57760/sciencedb.015m5

The DOI `10.57760/sciencedb.015m5` is the dataset DOI (not a preprint DOI) and
is cited as reference [25] in the manuscript. CSTR:
`31253.11.sciencedb.015m5`.

## License

MIT — see `LICENSE`.
