import os
"""
GNB-IE oral step-down -- v3.0 : MIC-EXPLICIT build
================================================================
Purpose of this build
----------------------------------------------------------------
v2.0 produced verdicts that DEPEND on a single MIC input (the ECOFF).
Two of those inputs were flagged as unverified:
    * P. aeruginosa / levofloxacin  ECOFF 2.0 mg/L  (from a 2014 EUCAST doc)
    * Enterobacterales / levofloxacin ECOFF 0.25 mg/L
If either value is wrong, the verdict flips.  That is a fatal dependency for
a paper whose central claim is a verdict.

v3.0 removes the dependency by making MIC an EXPLICIT axis of the result
instead of an input assumption.  We report, for every regimen and every
PK/PD target, the full PTA-vs-MIC curve over the entire scenario grid.
Any reader (and any reviewer) can then read the verdict off the curve for
whatever MIC value they accept -- ECOFF, clinical breakpoint, or their own
laboratory's MIC.  The paper no longer asserts a verdict that depends on a
number we cannot verify; it asserts a CURVE, from which verdicts follow
conditionally.

Consequences that must be stated in the manuscript
----------------------------------------------------------------
1. For P. aeruginosa, the claim "0 % target attainment" is TRUE ONLY AT
   MIC = 2.0 mg/L.  At MIC = 0.5 mg/L the strictest and the conventional
   targets fail but the most permissive target (fAUC24/MIC >= 72) is
   borderline.  The v2.0 wording must therefore be narrowed.
2. The Enterobacterales claim is robust: it holds across the whole
   MIC range that any plausible ECOFF or breakpoint could take.

Added axes relative to v2.0
----------------------------------------------------------------
* MIC grid          2^-5 ... 2^3 mg/L (0.03125 - 8), explicit output axis
* target            AUC24/MIC >= 500 added (conservative proxy for the
                    high-inoculum / vegetation-penetration effect, which
                    cannot be modelled explicitly and was previously
                    handled only by AUC >= 250)
* boundary cells    P. aeruginosa at MIC 0.25/0.5/1.0/2.0 reported
                    individually with scenario-min / scenario-max

Numerical provenance is unchanged from v2.0 (FDA label primary exposure
model; popPK secondary).  No new clinical data, no outcome model.

Verification status of the MIC inputs
----------------------------------------------------------------
Updated 2026-10-12.  The two flagged ECOFF values, and in fact all six
ECOFF values used to read verdicts off the curves, were read directly off
the EUCAST MIC distribution database (https://mic.eucast.org/search/,
Method = MIC) and are reproduced by fetch_eucast_ecoff.py in this
repository:

    E. coli            ciprofloxacin 0.06  mg/L  (15,667 obs / 53 sources)
    E. coli            levofloxacin  0.125 mg/L  ( 2,778 obs /  5 sources)
    K. pneumoniae      ciprofloxacin 0.125 mg/L  ( 3,788 obs / 40 sources)
    K. pneumoniae      levofloxacin  0.25  mg/L  ( 3,674 obs /  6 sources)
    P. aeruginosa      ciprofloxacin 0.5   mg/L  (26,996 obs / 52 sources)
    P. aeruginosa      levofloxacin  2.0   mg/L  (14,871 obs / 11 sources)

Note that the ECOFF is species-specific within the Enterobacterales: the
0.125 / 0.25 mg/L pair that earlier drafts applied to the order as a whole
is in fact the K. pneumoniae pair, one doubling dilution above E. coli.
The MIC axis is explicit in this build, so such a correction changes which
cell of the grid a verdict is read from but never the curve itself.
"""

import numpy as np, json

rng = np.random.default_rng(20261008)
N = 200_000

MIC_GRID = [2.0 ** k for k in range(-5, 4)]          # 0.03125 ... 8

TARGETS = [
    ('fAUC24/MIC>=72',  'free',   72.0),
    ('AUC24/MIC>=125',  'total', 125.0),
    ('AUC24/MIC>=250',  'total', 250.0),
    ('AUC24/MIC>=500',  'total', 500.0),   # NEW conservative proxy
]
THRESH = (0.80, 0.90, 0.95)

DRUGS = {
    'cipro': dict(fu=(0.60, 0.70, 0.80),
                  label_auc={500: 27.4, 750: 31.6},
                  label_sd={500: None, 750: None},
                  poppk=dict(mode='iv', CL=25.4, iiv=0.678, F=0.70), fe=0.45),
    'levo': dict(fu=(0.62, 0.69, 0.76),
                 label_auc={500: 47.5, 750: 90.7},
                 label_sd={500: 6.7, 750: 17.6},
                 poppk=dict(mode='oral', CL=8.32, iiv=0.226, F=0.99),
                 renal={'>80': 47.9, '50-80': 95.6, '20-49': 182.1, '<20': 263.5}),
}
REGIMENS = [('cipro', 500, 2), ('cipro', 750, 2), ('levo', 500, 1), ('levo', 750, 1)]


def lognorm(mean, sd, n):
    cv2 = (sd / mean) ** 2
    s2 = np.log(1.0 + cv2)
    return np.exp(rng.normal(np.log(mean) - 0.5 * s2, np.sqrt(s2), n))


def auc_label(drug, dose, n, cv_override=None):
    d = DRUGS[drug]
    m = d['label_auc'][dose]; sd = d['label_sd'][dose]
    if sd is None:
        sd = m * (cv_override if cv_override else 0.25)
    return lognorm(m, sd, n)


def auc_poppk(drug, dose, nday, n, egfr=100):
    d = DRUGS[drug]; p = d['poppk']
    if drug == 'cipro':
        cl = np.exp(rng.normal(np.log(p['CL']), p['iiv'], n))
        cl = cl * (1 - d['fe'] + d['fe'] * egfr / 100.0)
        return p['F'] * dose * nday / cl
    clf = np.exp(rng.normal(np.log(p['CL']), p['iiv'], n))
    clf = clf * (1 - 0.66 + 0.66 * egfr / 100.0)
    return p['F'] * dose * nday / clf


def exposure_models(drug, dose, nday):
    """All exposure scenarios carried forward from v2.0."""
    m = {'label': auc_label(drug, dose, N),
         'poppk': auc_poppk(drug, dose, nday, N)}
    if drug == 'cipro':
        m['label_cv15'] = auc_label(drug, dose, N, cv_override=0.15)
        m['label_cv40'] = auc_label(drug, dose, N, cv_override=0.40)
        for name, val in (('cand_ss_table', 31.6), ('cand_iv_equiv', 32.9),
                          ('cand_single_x2', 40.4), ('cand_poppk_icu', 41.3)):
            if dose == 750:
                m[name] = lognorm(val, val * 0.25, N)
    else:
        m['label'] = auc_label(drug, dose, N)
    return m


def pta(auc, fu, mic, target, basis):
    a = auc * fu if basis == 'free' else auc
    return float(np.mean(a / mic >= target))


def mc_halfwidth(p, n=N):
    return round(1.96 * np.sqrt(max(p, 1e-9) * (1 - min(p, 1 - 1e-9)) / n) * 100, 3)


res = {}
P = print

# ============================================================ STEP 1  CURVES
P('=' * 86)
P('STEP 1   PTA vs MIC  --  the MIC value is now an OUTPUT AXIS, not an assumption')
P('=' * 86)
curves = {}
for drug, dose, nday in REGIMENS:
    models = exposure_models(drug, dose, nday)
    per_t = {}
    for tname, basis, t in TARGETS:
        row = {}
        for mic in MIC_GRID:
            ps = [pta(a, fu, mic, t, basis)
                  for a in models.values() for fu in DRUGS[drug]['fu']]
            row[f'{mic:g}'] = [round(min(ps) * 100, 1), round(max(ps) * 100, 1)]
        per_t[tname] = row
        s = '  '.join(f'{m:g}:{v[0]:.0f}-{v[1]:.0f}' for m, v in
                      zip(MIC_GRID, row.values()))
        P(f'{drug:6s}{dose:4d} x{nday} | {tname:15s} {s}')
    curves[f'{drug} {dose} x{nday}'] = per_t
res['PTA_VS_MIC'] = curves

# ============================================================ STEP 2  BREAKPOINTS
P('\n' + '=' * 86)
P('STEP 2   Highest MIC attaining PTA >= 90 %  (scenario min - max)')
P('=' * 86)
P('v3.1b: a scenario in which NO MIC on the grid attains 90 % is kept as None.')
P('       EVERY therefore becomes None if any scenario fails altogether;')
P('       SOME is the maximum over the scenarios that do attain it.')
bps = {}
for drug, dose, nday in REGIMENS:
    models = exposure_models(drug, dose, nday)
    row = {}
    for tname, basis, t in TARGETS:
        got = []
        for a in models.values():
            for fu in DRUGS[drug]['fu']:
                b = None
                for mic in MIC_GRID:
                    if pta(a, fu, mic, t, basis) >= 0.90:
                        b = mic
                    else:
                        break
                got.append(b)                      # None = this scenario never attains 90 %
        vals = [v for v in got if v is not None]
        every = min(got) if all(v is not None for v in got) else None
        some = max(vals) if vals else None
        row[tname] = [every, some]
        fmt = lambda v: 'none' if v is None else f'{v:g}'
        P(f'{drug:6s}{dose:4d} x{nday} | {tname:15s} -> '
          f'{fmt(every)} - {fmt(some)} mg/L')
    bps[f'{drug} {dose} x{nday}'] = row
res['BREAKPOINT_PTA90'] = bps

# ============================================================ STEP 3  P. AERUGINOSA BOUNDARY
P('\n' + '=' * 86)
P('STEP 3   P. aeruginosa boundary -- the cell whose verdict depends on the MIC input')
P('=' * 86)
BOUND = {
    # MIC anchors now follow the EUCAST aggregated ECOFF values confirmed
    # against the EUCAST ECOFF change log (2026-10-12):
    #   E. coli        cipro 0.0625 (reported as 0.06) , levo 0.125
    #   K. pneumoniae  cipro 0.125                    , levo 0.25
    'E. coli+K. pneumoniae / cipro': ('cipro', (0.0625, 0.125, 0.25, 0.5)),
    'E. coli+K. pneumoniae / levo':  ('levo',  (0.125, 0.25, 0.5, 1.0)),
    'P. aeruginosa / cipro':    ('cipro', (0.25, 0.5, 1.0, 2.0)),
    'P. aeruginosa / levo':     ('levo',  (0.25, 0.5, 1.0, 2.0)),
}
bound = {}
for name, (drug, mics) in BOUND.items():
    nday = 2 if drug == 'cipro' else 1
    for dose in (500, 750):
        models = exposure_models(drug, dose, nday)
        key = f'{name} | {drug} {dose} mg'
        row = {}
        for tname, basis, t in TARGETS:
            cells = {}
            for mic in mics:
                ps = [pta(a, fu, mic, t, basis)
                      for a in models.values() for fu in DRUGS[drug]['fu']]
                lo, hi = min(ps), max(ps)
                if lo >= 0.90:
                    v = 'ELIGIBLE in every scenario'
                elif hi < 0.90:
                    v = 'NOT ELIGIBLE in any scenario'
                else:
                    v = 'CONDITIONAL'
                cells[f'{mic:g}'] = dict(pta_min=round(lo * 100, 1),
                                         pta_max=round(hi * 100, 1), verdict=v)
            row[tname] = cells
            P(f'{key:38s} {tname:15s} ' +
              '  '.join(f'MIC{m}={c["pta_min"]:.0f}-{c["pta_max"]:.0f}%'
                        f'[{c["verdict"][0]}]' for m, c in
                        zip(mics, cells.values())))
        bound[key] = row
res['BOUNDARY_CELLS'] = bound

# ============================================================ STEP 4  INVARIANCE OF THE CENTRAL CLAIM
P('\n' + '=' * 86)
P('STEP 4   Over what MIC range does each central claim hold?')
P('=' * 86)
claim = {}
for drug, dose, nday in REGIMENS:
    models = exposure_models(drug, dose, nday)
    row = {}
    for tname, basis, t in TARGETS:
        # v3.1: separate boolean arrays. SOME must use max(ps), not min(ps).
        ok_every, ok_some = [], []
        for mic in MIC_GRID:
            ps = [pta(a, fu, mic, t, basis)
                  for a in models.values() for fu in DRUGS[drug]['fu']]
            ok_every.append(min(ps) >= 0.90)
            ok_some.append(max(ps) >= 0.90)
        hi_every = None
        for mic, o in zip(MIC_GRID, ok_every):
            if o:
                hi_every = mic
            else:
                break
        hi_some = None
        for mic, o in zip(MIC_GRID, ok_some):
            if o:
                hi_some = mic
            else:
                break
        row[tname] = dict(highest_MIC_eligible_in_EVERY_scenario=hi_every,
                          highest_MIC_eligible_in_SOME_scenario=hi_some)
        P(f'{drug:6s}{dose:4d} x{nday} | {tname:15s} '
          f'EVERY <= {row[tname]["highest_MIC_eligible_in_EVERY_scenario"]} , '
          f'SOME <= {row[tname]["highest_MIC_eligible_in_SOME_scenario"]} mg/L')
    claim[f'{drug} {dose} x{nday}'] = row
res['CLAIM_MIC_RANGE'] = claim

_here = os.path.dirname(os.path.abspath(__file__))
_base = os.path.basename(__file__)
if _base.startswith('Supplemental_Code_'):
    # This copy is supplied as a journal supplementary file; its output must
    # carry the matching Supplemental_Results_* name.
    _outname = 'Supplemental_Results_' + _base[len('Supplemental_Code_'):-3] + '.json'
else:
    _outname = _base[:-3] + '_results.json'
json.dump(res, open(os.path.join(_here, _outname), 'w'), indent=1)
P('\nWROTE ' + _outname)
