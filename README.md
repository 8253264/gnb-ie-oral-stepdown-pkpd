# Oral fluoroquinolone step-down for non-HACEK Gram-negative infective endocarditis: a Monte Carlo PK/PD eligibility analysis

Simulation code and result files for the manuscript.

**No new clinical data were generated.** All exposure inputs come from published
pharmacokinetic sources (FDA product labels and published population PK models);
all MIC inputs come from the EUCAST MIC distribution database. The output is a
probability-of-target-attainment (PTA) curve and a PK/PD eligibility verdict per
regimen, MIC value and PK/PD target.

## What the analysis asks

The 2026 American Heart Association statement on infective endocarditis lists an
oral continuation regimen for HACEK organisms and for *Pseudomonas aeruginosa*,
but none for non-HACEK Enterobacterales. This analysis asks a narrower,
pharmacological question: **which oral fluoroquinolone regimens can actually
attain a PK/PD target at the MIC values that laboratories report as susceptible?**

Four regimens were evaluated:

| Regimen | Dosing |
|---|---|
| Ciprofloxacin | 500 mg twice daily |
| Ciprofloxacin | 750 mg twice daily |
| Levofloxacin | 500 mg once daily |
| Levofloxacin | 750 mg once daily |

against four PK/PD targets:

| Target | Basis |
|---|---|
| fAUC₂₄/MIC ≥ 72 | free drug; modified Delphi consensus |
| AUC₂₄/MIC ≥ 125 | conventional |
| AUC₂₄/MIC ≥ 250 | clinically validated in Gram-negative bloodstream infection |
| AUC₂₄/MIC ≥ 500 | extreme upper bound / high-exposure proxy |

MIC is an **explicit output axis**, not an input assumption: the PTA-versus-MIC
curve is reported over the whole scenario grid, so a verdict can be read off the
curve at whatever MIC a reader accepts (ECOFF, clinical breakpoint, or a local
laboratory value).

## Files

| File | Contents |
|---|---|
| `gnb_pkpd_v3.py` | Monte Carlo simulation. 200,000 subjects per scenario; two parallel exposure models (label steady-state AUC converted to log-normal distributions, and published population PK models). Writes `gnb_pkpd_v3_results.json`. |
| `gnb_pkpd_v3_results.json` | Full result set: PTA grids, PK/PD breakpoints, P. aeruginosa boundary cells. |
| `gnlog_v3.txt` | Plain-text log of the simulation run. |
| `make_figures.py` | Generates Figures 1-4 from the result JSON (vector PDF + 600 dpi PNG). |
| `fetch_eucast_ecoff.py` | Retrieves the six ECOFF values, with observation and source counts, from the EUCAST MIC distribution database. |
| `Figure_1.png` … `Figure_4.png` | 600 dpi previews of the manuscript figures. |

## Reproducing

Requires Python 3 and NumPy.

```bash
pip install numpy
python gnb_pkpd_v3.py     # writes gnb_pkpd_v3_results.json  (a few minutes)
python make_figures.py    # writes Figure_1..4 .pdf and .png
python fetch_eucast_ecoff.py   # re-derives the ECOFF values used
```

The random seed is fixed (`default_rng(20261008)`), so the results are exactly
reproducible.

## MIC inputs (epidemiological cut-offs)

All six values were read directly from the EUCAST MIC distribution database
(https://mic.eucast.org/search/, Method = MIC) on 2026-10-12, together with the
number of observations and contributing sources. `fetch_eucast_ecoff.py`
reproduces the retrieval.

| Species | Agent | ECOFF (mg/L) | Observations | Sources | 95% CI |
|---|---|---|---|---|---|
| *Escherichia coli* | Ciprofloxacin | 0.06 | 15,667 | 53 | 0.03 – 0.06 |
| *Escherichia coli* | Levofloxacin | 0.125 | 2,778 | 5 | 0.03 – 0.25 |
| *Klebsiella pneumoniae* | Ciprofloxacin | 0.125 | 3,788 | 40 | 0.06 – 0.25 |
| *Klebsiella pneumoniae* | Levofloxacin | 0.25 | 3,674 | 6 | 0.06 – 0.25 |
| *Pseudomonas aeruginosa* | Ciprofloxacin | 0.5 | 26,996 | 52 | 0.5 – 1 |
| *Pseudomonas aeruginosa* | Levofloxacin | 2.0 | 14,871 | 11 | 1 – 4 |

**Note on species specificity.** The ECOFF is species-specific within the
Enterobacterales. The 0.125 / 0.25 mg/L pair that is sometimes applied to the
order as a whole is in fact the *K. pneumoniae* pair; *E. coli* sits one doubling
dilution below it. Both species are therefore modelled and reported separately.

## Notes

- The `fetch_eucast_ecoff.py` output depends on the current state of the EUCAST
  database; values are subject to revision by EUCAST.
- Clinical breakpoints cited in the manuscript are from EUCAST v16.1 and from
  **CLSI M100 Ed36** (read from Tables 2A-1 and 2B-1 of the standard itself).
  Those documents are not redistributed here — please consult EUCAST and CLSI
  directly.

## Citation

The manuscript is under review. If you use this code before it appears, please
cite this archived release:

> Chen H, Li Y, Hua B, Sun Y. *Oral fluoroquinolone step-down for non-HACEK
> Gram-negative infective endocarditis: a Monte Carlo PK/PD eligibility
> analysis.* Simulation code and result files (version 1.0.0). Science Data
> Bank; 2026. https://doi.org/10.57760/sciencedb.015m5

The archived, DOI-bearing copy is deposited in the **Science Data Bank**
(ScienceDB, https://www.scidb.cn/, operated by the Computer Network Information
Center of the Chinese Academy of Sciences), which issued the DOI and a CSTR on
publication (19 October 2026):

- DOI: `10.57760/sciencedb.015m5` — https://doi.org/10.57760/sciencedb.015m5
- CSTR: `31253.11.sciencedb.015m5` — https://cstr.cn/31253.11.sciencedb.015m5
- Landing page: https://www.scidb.cn/detail?dataSetId=35c0b02b9fc84a32a5586a206e6220db

The same DOI appears in the manuscript's Data availability statement, where the
archived release is cited as reference [25].

Corresponding author: Huiwen Chen (陈慧文), Department of Cardiac Surgery,
Wuhan Asia Heart Hospital. ORCID 0009-0006-4948-0174.

## License

MIT — see `LICENSE`.
