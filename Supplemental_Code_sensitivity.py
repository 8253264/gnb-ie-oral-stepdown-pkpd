"""
Supplemental code — sensitivity analyses (Section 3.6)
======================================================
Manuscript: "Reproducible pharmacokinetic/pharmacodynamic breakpoints for oral
fluoroquinolone step-down in Gram-negative infective endocarditis".

WHAT THIS FILE DOES
-------------------
Runs the three sensitivity analyses named in Section 3.6, using the SAME
exposure models as the primary simulation.  Self-contained: the parameter
block below is copied verbatim from Supplemental_Code_gnb_pkpd_v3_1.py.

  S1  Unbound fraction (fu) varied one value at a time.
  S2  Ciprofloxacin label AUC variability (CV 0.15 / 0.25 / 0.40).
  S3  PTA criterion varied over 80 / 90 / 95 %.
      For S3 the eligibility classes of the organism-anchored decision cells
      that make up Tables 4 and 5 are recomputed at each threshold, so that
      the statement about threshold dependence in Section 3.6 is quantitative
      rather than rhetorical.

DEFINITIONS
-----------
* "EVERY breakpoint" = highest MIC on the two-fold grid at which PTA >= thr
  in every scenario of the grid (exposure models x unbound fractions).
* Eligibility class at a given MIC: E if the scenario minimum >= thr,
  N if the scenario maximum < thr, otherwise C (conditional).

RUN
---
    python Supplemental_Code_sensitivity.py
Writes Supplemental_Results_sensitivity.json next to this file.
Python 3 + NumPy, fixed seed 20261008, 200,000 subjects per scenario.
"""

import os
import json
import numpy as np

rng = np.random.default_rng(20261008)
N = 200_000
MIC_GRID = [2.0 ** k for k in range(-5, 4)]
TARGETS = [
    ('fAUC24/MIC>=72',  'free',   72.0),
    ('AUC24/MIC>=125',  'total', 125.0),
    ('AUC24/MIC>=250',  'total', 250.0),
    ('AUC24/MIC>=500',  'total', 500.0),
]
THRESHOLDS = (0.80, 0.90, 0.95)
DRUGS = {
    'cipro': dict(fu=(0.60, 0.70, 0.80),
                  label_auc={500: 27.4, 750: 31.6},
                  label_sd={500: None, 750: None},
                  poppk=dict(mode='iv', CL=25.4, iiv=0.678, F=0.70), fe=0.45),
    'levo': dict(fu=(0.62, 0.69, 0.76),
                 label_auc={500: 47.5, 750: 90.7},
                 label_sd={500: 6.7, 750: 17.6},
                 poppk=dict(mode='oral', CL=8.32, iiv=0.226, F=0.99)),
}
REGIMENS = [('cipro', 500, 2), ('cipro', 750, 2), ('levo', 500, 1), ('levo', 750, 1)]


def lognorm(mean, sd, n):
    cv2 = (sd / mean) ** 2
    s2 = np.log(1.0 + cv2)
    return np.exp(rng.normal(np.log(mean) - 0.5 * s2, np.sqrt(s2), n))


def auc_label(drug, dose, n, cv_override=None):
    d = DRUGS[drug]
    m = d['label_auc'][dose]
    sd = d['label_sd'][dose]
    if sd is None:
        sd = m * (cv_override if cv_override else 0.25)
    return lognorm(m, sd, n)


def auc_poppk(drug, dose, nday, n):
    d = DRUGS[drug]
    p = d['poppk']
    if drug == 'cipro':
        cl = np.exp(rng.normal(np.log(p['CL']), p['iiv'], n))
        return p['F'] * dose * nday / cl
    clf = np.exp(rng.normal(np.log(p['CL']), p['iiv'], n))
    return p['F'] * dose * nday / clf


def exposure_models(drug, dose, nday):
    m = {'label': auc_label(drug, dose, N),
         'poppk': auc_poppk(drug, dose, nday, N)}
    if drug == 'cipro':
        m['label_cv15'] = auc_label(drug, dose, N, cv_override=0.15)
        m['label_cv40'] = auc_label(drug, dose, N, cv_override=0.40)
        for name, val in (('cand_ss_table', 31.6), ('cand_iv_equiv', 32.9),
                          ('cand_single_x2', 40.4), ('cand_poppk_icu', 41.3)):
            if dose == 750:
                m[name] = lognorm(val, val * 0.25, N)
    return m


def pta_array(aucs, fus, mic, basis, target):
    return [float(np.mean(((a * f if basis == 'free' else a) / mic >= target)))
            for a in aucs for f in fus]


def every_bp(aucs, fus, basis, target, thr):
    hit = None
    for mic in MIC_GRID:
        if min(pta_array(aucs, fus, mic, basis, target)) >= thr:
            hit = mic
        else:
            break
    return hit


def fmt(v):
    return 'none' if v is None else f'{v:g}'


def dilutions(a, b):
    if a is None or b is None:
        return None
    return int(abs(np.log2(a / b)))


P = print
res = {}
P('=' * 86)
P('SECTION 3.6  SENSITIVITY ANALYSES  (seed 20261008, N=200,000)')
P('=' * 86)

# ------------------------------------------------------------------ S1  fu
P('\nS1  Unbound fraction varied one value at a time (EVERY breakpoint)')
s1 = []
changed = 0
for drug, dose, nday in REGIMENS:
    models = exposure_models(drug, dose, nday)
    aucs = list(models.values())
    for tname, basis, t in TARGETS:
        vals = {str(fu): fmt(every_bp(aucs, (fu,), basis, t, 0.90))
                for fu in DRUGS[drug]['fu']}
        bps = [every_bp(aucs, (fu,), basis, t, 0.90) for fu in DRUGS[drug]['fu']]
        spread = 0 if any(b is None for b in bps) else int(np.log2(max(bps) / min(bps)))
        cells_have_effect = len(set(vals.values())) > 1
        changed += int(cells_have_effect)
        s1.append(dict(regimen=f'{drug} {dose} x{nday}', target=tname,
                       by_fu=vals, spread_dilutions=spread,
                       breakpoint_changed=cells_have_effect))
        P(f'  {drug:6s}{dose:4d} | {tname:15s} ' +
          '  '.join(f'fu={k}:{v}' for k, v in vals.items()) +
          f' | spread {spread} dilution(s)')
res['S1_PROTEIN_BINDING'] = dict(cells=len(s1), cells_changed=changed, detail=s1)
P(f'  -> {changed} of {len(s1)} regimen x target cells change; '
  f'maximum spread {max(x["spread_dilutions"] for x in s1)} dilution(s).')

# ------------------------------------------------------------------ S2  CV
P('\nS2  Ciprofloxacin 500 mg label AUC variability, fAUC24/MIC >= 72')
s2 = []
for cv in (0.15, 0.25, 0.40):
    lab = auc_label('cipro', 500, N, cv_override=cv)
    pop = auc_poppk('cipro', 500, 2, N)
    aucs = [lab, pop]
    fus = DRUGS['cipro']['fu']
    e = every_bp(aucs, fus, 'free', 72.0, 0.90)
    some = None
    for mic in MIC_GRID:
        if max(pta_array(aucs, fus, mic, 'free', 72.0)) >= 0.90:
            some = mic
        else:
            break
    s2.append(dict(cv=cv, every=fmt(e), some=fmt(some)))
    P(f'  CV {cv:.2f} -> EVERY {fmt(e)} ; SOME {fmt(some)}')
res['S2_LABEL_CV'] = s2

# ------------------------------------------------------------------ S3  threshold
P('\nS3  PTA criterion 80 / 90 / 95 %: EVERY breakpoints')
s3_bp = {}
for thr in THRESHOLDS:
    for drug, dose, nday in REGIMENS:
        models = exposure_models(drug, dose, nday)
        aucs = list(models.values())
        key = f'{drug} {dose} x{nday}'
        s3_bp.setdefault(key, {})[str(thr)] = {
            tname: fmt(every_bp(aucs, DRUGS[drug]['fu'], basis, t, thr))
            for tname, basis, t in TARGETS}
for k, v in s3_bp.items():
    P(f'  {k:16s} ' + ' | '.join(f'{thr}: ' +
      ','.join(v[str(thr)][t[0]] for t in TARGETS) for thr in THRESHOLDS))

n_cells = sum(len(v[str(0.90)]) for v in s3_bp.values())
changed_vs90 = sum(1 for v in s3_bp.values() for t in v[str(0.90)]
                   if v[str(0.80)][t] != v[str(0.90)][t]
                   or v[str(0.95)][t] != v[str(0.90)][t])
max_shift = 0
for v in s3_bp.values():
    for t in v[str(0.90)]:
        base = v[str(0.90)][t]
        for thr in (0.80, 0.95):
            x, y = v[str(thr)][t], base
            if x != 'none' and y != 'none':
                max_shift = max(max_shift, abs(int(round(np.log2(float(x) / float(y))))))

res['S3_THRESHOLD_BREAKPOINTS'] = dict(
    breakpoints=s3_bp, cells=n_cells,
    cells_changed_at_80_or_95=changed_vs90,
    max_shift_dilutions=max_shift)

# Decision cells of Tables 4 and 5 (regimen, MIC, targets)
T4 = [('cipro', 500, 0.0625), ('cipro', 500, 0.125),
      ('cipro', 750, 0.0625), ('cipro', 750, 0.125), ('cipro', 750, 0.25),
      ('levo', 500, 0.125), ('levo', 500, 0.25),
      ('levo', 750, 0.125), ('levo', 750, 0.25), ('levo', 750, 0.5)]
T5 = [('cipro', 500, 0.5), ('cipro', 750, 0.5), ('levo', 500, 2.0),
      ('levo', 750, 2.0), ('levo', 750, 1.0), ('levo', 750, 0.25)]


def classify(drug, dose, mic, basis, target, thr):
    nday = 2 if drug == 'cipro' else 1
    models = exposure_models(drug, dose, nday)
    ps = pta_array(list(models.values()), DRUGS[drug]['fu'], mic, basis, target)
    if min(ps) >= thr:
        return 'E'
    if max(ps) < thr:
        return 'N'
    return 'C'


classes = {}
for name, cells, tset in (('Table 4', T4, TARGETS), ('Table 5', T5, TARGETS[:3])):
    for drug, dose, mic in cells:
        for tname, basis, t in tset:
            key = f'{name}|{drug} {dose}|MIC {mic:g}|{tname}'
            classes[key] = {str(thr): classify(drug, dose, mic, basis, t, thr)
                            for thr in THRESHOLDS}
n_dec = len(classes)
changed_dec = sum(1 for v in classes.values()
                  if v['0.8'] != v['0.9'] or v['0.95'] != v['0.9'])
# Principal calls (the conclusions of the paper)
principal = [
    ('Table 4|levo 750|MIC 0.125|fAUC24/MIC>=72', 'E'),
    ('Table 4|levo 750|MIC 0.125|AUC24/MIC>=250', 'E'),
    ('Table 4|levo 750|MIC 0.25|AUC24/MIC>=125', 'E'),
    ('Table 4|cipro 750|MIC 0.25|AUC24/MIC>=125', 'N'),
    ('Table 5|levo 750|MIC 2|AUC24/MIC>=125', 'N'),
    ('Table 5|levo 750|MIC 1|AUC24/MIC>=125', 'N'),
]
principal_ok = all(classes[k][str(thr)] == want
                   for k, want in principal for thr in THRESHOLDS)
res['S3_THRESHOLD_CLASSES'] = dict(
    decision_cells=n_dec, cells_changed_at_80_or_95=changed_dec,
    principal_calls_unchanged=bool(principal_ok),
    principal_calls=[dict(cell=k, expected=w,
                          at_80=classes[k]['0.8'], at_90=classes[k]['0.9'],
                          at_95=classes[k]['0.95']) for k, w in principal],
    detail=classes)
P(f'\n  Breakpoint cells: {n_cells}; changed at 80 % or 95 %: {changed_vs90}; '
  f'maximum shift {max_shift} dilution(s).')
P(f'  Decision cells in Tables 4-5: {n_dec}; class changed at 80 % or 95 %: '
  f'{changed_dec}.')
P(f'  All six principal eligibility calls unchanged at all three thresholds: '
  f'{principal_ok}.')
for x in res['S3_THRESHOLD_CLASSES']['principal_calls']:
    P(f'    {x["cell"]:48s} 80%={x["at_80"]} 90%={x["at_90"]} 95%={x["at_95"]}')

# ------------------------------------------------------------------ S4  candidate AUC
# Table 3: the four ciprofloxacin 750 mg candidate AUC values carried through
# the grid are tabulated one at a time, at the central unbound fraction
# (fu = 0.70), so that the reader can see which target is assumption-sensitive.
P('\nS4  Ciprofloxacin 750 mg candidate AUC variants (Table 3), fu = 0.70')
s4 = []
for mean_auc, source in ((31.6, 'label steady-state table'),
                         (32.9, 'label IV-equivalence statement'),
                         (40.4, '2 x single-dose AUC0-8'),
                         (41.3, 'ICU population PK model')):
    auc = lognorm(mean_auc, mean_auc * 0.25, N)
    med = float(np.median(auc))
    bp = {tname: fmt(every_bp([auc], (0.70,), basis, t, 0.90))
          for tname, basis, t in TARGETS[:3]}
    s4.append(dict(input_mean=mean_auc, median=round(med, 1), source=source,
                   breakpoints=bp))
    P(f'  mean {mean_auc:5.1f} (median {med:5.1f}) | ' +
      '  '.join(f'{k}={v}' for k, v in bp.items()))
res['S4_CANDIDATE_AUC'] = s4

out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'Supplemental_Results_sensitivity.json')
json.dump(res, open(out, 'w', encoding='utf8'), indent=1)
P('WROTE', os.path.basename(out))
