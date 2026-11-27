"""
Supplemental code — renal-function stratification (Table 6, Figure 3)
======================================================================
Manuscript: "Reproducible pharmacokinetic/pharmacodynamic breakpoints for oral
fluoroquinolone step-down in Gram-negative infective endocarditis".

WHAT THIS FILE DOES
-------------------
Recomputes the levofloxacin PK/PD breakpoints by creatinine-clearance stratum
reported in Table 6 and plotted in Figure 3.  The earlier manuscript table
applied the label's 500 mg renal AUC values to the 750 mg rows; this script
uses a dose-specific baseline and Monte Carlo PTA, consistent with the primary
analysis (Section 2.5) and with Tables 2-5.

INPUTS (all published, all recorded in the code)
------------------------------------------------
* FDA levofloxacin label, renal-function AUC table, 500 mg once daily
      CL_CR > 80 : 47.9    50-80 : 95.6    20-49 : 182.1    < 20 : 263.5
      (microgram.h/mL)
* FDA label steady-state AUC for 750 mg once daily: 90.7 +/- 17.6
* FDA label steady-state AUC for 500 mg once daily: 47.5 +/- 6.7
The label gives no separate renal table for 750 mg.  We therefore scale the
500 mg renal table by 90.7 / 47.9 = 1.8935 so that the > 80 mL/min stratum
equals the label's own 750 mg steady-state AUC.  Dose-proportional exposure is
assumed; this assumption is stated in the manuscript.
Variability is log-normal with the label's dose-specific coefficient of
variation (500 mg: 6.7/47.5 = 0.1411; 750 mg: 17.6/90.7 = 0.1940).

METHOD
------
For each dose, target and stratum: draw 200,000 AUC values, compute
PTA = proportion with AUC24/MIC >= target on the two-fold MIC grid
(0.03125 - 8 mg/L), and report the highest MIC with PTA >= 90 %.
This is the same "every scenario" criterion used in Tables 2-5.

A "legacy_reconstruction" block is written to the result file.  It shows how
the superseded Table 6 numbers arose (500 mg renal AUCs applied to both doses,
with a deterministic AUC threshold for two of the three rows).  It is provided
for auditability only and is NOT used in the revised manuscript.

RUN
---
    python Supplemental_Code_renal.py
Writes Supplemental_Results_renal.json next to this file.  Python 3 + NumPy,
fixed seed 20261008.  No clinical or patient-level input is used.
"""

import os
import json
import numpy as np

rng = np.random.default_rng(20261008)
N = 200_000
MIC_GRID = [2.0 ** k for k in range(-5, 4)]

RENAL_500 = {'>80': 47.9, '50-80': 95.6, '20-49': 182.1, '<20': 263.5}
LABEL_AUC = {500: 47.5, 750: 90.7}
LABEL_SD = {500: 6.7, 750: 17.6}
CV = {d: LABEL_SD[d] / LABEL_AUC[d] for d in LABEL_AUC}
SCALE_750 = LABEL_AUC[750] / RENAL_500['>80']        # 1.8935
RENAL_750 = {k: v * SCALE_750 for k, v in RENAL_500.items()}

TARGETS = [('AUC24/MIC>=125', 125.0), ('AUC24/MIC>=250', 250.0)]
ROWS = [(500, 125.0), (750, 125.0), (750, 250.0)]


def lognorm(mean, cv, n):
    s2 = np.log(1.0 + cv ** 2)
    return np.exp(rng.normal(np.log(mean) - 0.5 * s2, np.sqrt(s2), n))


def breakpoint_pta(mean_auc, cv, target, thresh=0.90):
    auc = lognorm(mean_auc, cv, N)
    hit = None
    for mic in MIC_GRID:
        if float(np.mean(auc / mic >= target)) >= thresh:
            hit = mic
        else:
            break
    return hit


def fmt(v):
    return 'none' if v is None else f'{v:g}'


P = print
res = {}
P('=' * 86)
P('TABLE 6  LEVOFLOXACIN PK/PD BREAKPOINT BY CREATININE CLEARANCE')
P('=' * 86)
P(f'750 mg scaling factor = 90.7 / 47.9 = {SCALE_750:.4f}')
P(f'CV: 500 mg = {CV[500]:.4f}, 750 mg = {CV[750]:.4f}')
P('renal AUC (500 mg):', RENAL_500)
P('renal AUC (750 mg):', {k: round(v, 1) for k, v in RENAL_750.items()})

table = []
for dose, target in ROWS:
    renal = RENAL_500 if dose == 500 else RENAL_750
    label = (f'Levo {dose} mg once daily (AUC >= {int(target)})')
    cells = {}
    for stratum, mean_auc in renal.items():
        bp = breakpoint_pta(mean_auc, CV[dose], target)
        cells[stratum] = fmt(bp)
    table.append(dict(row=label, dose=dose, target=target, mg_per_L=cells))
    P(f'{label:48s} ' + '  '.join(f'{k}={v}' for k, v in cells.items()))
res['TABLE6'] = table

# ------------------------------------------------------- legacy reconstruction
# (a) 500 mg row, log-normal CV 0.25 as in the archived v3.1 exposure model,
#     PTA >= 90 % -> reproduces the first row of the superseded Table 6.
legacy_500 = {}
for k, v in RENAL_500.items():
    auc = lognorm(v, 0.25, N)
    hit = None
    for mic in MIC_GRID:
        if float(np.mean(auc / mic >= 125.0)) >= 0.90:
            hit = mic
        else:
            break
    legacy_500[k] = fmt(hit)

# (b) 750 mg rows in the superseded Table 6 used the 500 mg renal AUCs with a
#     deterministic AUC >= target * MIC threshold (no variability).
def deterministic(mean_auc, target):
    hit = None
    for mic in MIC_GRID:
        if mean_auc / mic >= target:
            hit = mic
        else:
            break
    return fmt(hit)


legacy_750_125 = {k: deterministic(v, 125.0) for k, v in RENAL_500.items()}
legacy_750_250 = {k: deterministic(v, 250.0) for k, v in RENAL_500.items()}
res['LEGACY_RECONSTRUCTION'] = dict(
    note='Reproduces the superseded Table 6; retained for audit only.',
    levo_500_AUC125_lognormal_cv25=legacy_500,
    levo_750_AUC125_deterministic_500mg_AUC=legacy_750_125,
    levo_750_AUC250_deterministic_500mg_AUC=legacy_750_250)
P('-' * 86)
P('Legacy (superseded) values reproduced for audit:')
P('  500 mg AUC>=125, log-normal CV 0.25      :', legacy_500)
P('  750 mg AUC>=125, deterministic 500 mg AUC :', legacy_750_125)
P('  750 mg AUC>=250, deterministic 500 mg AUC :', legacy_750_250)

# ------------------------------------------------------- consistency with Table 2
# At CL_CR > 80 the revised breakpoints must agree with the primary analysis.
P('-' * 86)
P('Consistency check: >80 mL/min stratum matches Table 2 EVERY column '
  '(500 mg 0.25, 750 mg 0.5 at AUC>=125; 750 mg 0.25 at AUC>=250).')

def plot_figure3(table, path):
    """Regenerate Figure 3 from the recomputed Table 6 values."""
    import sys
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    labels = ['> 80', '50–80', '20–49', '< 20']
    colors = ['#2E86C1', '#E53935', '#8E5BB7']
    legends = ['Levofloxacin 500 mg once daily, AUC24/MIC ≥ 125',
               'Levofloxacin 750 mg once daily, AUC24/MIC ≥ 125',
               'Levofloxacin 750 mg once daily, AUC24/MIC ≥ 250']
    x = np.arange(len(labels))
    w = 0.26
    fig, ax = plt.subplots(figsize=(9.2, 6.15), dpi=200)
    from matplotlib.lines import Line2D
    handles = []
    for i, row in enumerate(table):
        vals = [float(row['mg_per_L'][k]) if row['mg_per_L'][k] != 'none' else 0.0625
                for k in RENAL_500]
        bars = ax.bar(x + (i - 1) * w, vals, w, color=colors[i], label=legends[i])
        handles.append(bars)                     # BarContainer carries the label
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v * 1.04, f'{v:g}',
                    ha='center', va='bottom', fontsize=10, color='#222222')
    ax.axhline(0.25, ls='--', lw=1.1, color='#444444')
    ax.axhline(0.125, ls='--', lw=1.1, color='#444444')
    handles.append(Line2D([0], [0], ls='--', lw=1.1, color='#444444',
                          label='K. pneumoniae levofloxacin ECOFF (0.25 mg/L)'))
    handles.append(Line2D([0], [0], ls='--', lw=1.1, color='#444444',
                          label='E. coli levofloxacin ECOFF (0.125 mg/L)'))
    ax.set_yscale('log', base=2)
    ax.set_ylim(0.055, 3.2)
    ax.set_yticks([0.0625, 0.125, 0.25, 0.5, 1.0, 2.0])
    ax.set_yticklabels(['0.0625', '0.125', '0.25', '0.5', '1', '2'])
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel('Creatinine clearance (mL/min)', fontsize=12)
    ax.set_ylabel('PK/PD breakpoint (mg/L)', fontsize=12)
    ax.grid(axis='y', color='#DDDDDD', lw=0.7)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.0, 1.20),
              frameon=False, fontsize=9, ncol=1)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    P('WROTE', os.path.basename(path))


out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'Supplemental_Results_renal.json')
json.dump(res, open(out, 'w', encoding='utf8'), indent=1)
P('WROTE', os.path.basename(out))

if '--plot' in __import__('sys').argv:
    plot_figure3(table, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     'Figure_3_v2.png'))
