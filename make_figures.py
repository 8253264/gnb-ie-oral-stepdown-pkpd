"""
Figures 1-4 for the IJAA submission
===================================
"Oral fluoroquinolone step-down for Gram-negative infective endocarditis:
 PK/PD breakpoints for Enterobacterales and Pseudomonas aeruginosa, and their
 relation to the 2026 American Heart Association oral continuation regimens"

Output per figure:  Figure_N.pdf (vector, Arial embedded)
                    Figure_N.png (600 dpi raster)

Provenance
----------
Every point plotted is the closed-form PTA of a log-normal exposure
distribution with the parameters of gnb_pkpd_v3.py.  The analytic model was
validated against the stored Monte-Carlo output in gnb_pkpd_v3_results.json:
largest discrepancy over all regimen x target x MIC cells = 0.103 percentage
points, inside the 0.180 pp Monte-Carlo 95 % half-width at N = 200 000.
Verdicts are recomputed here; none is transcribed by hand from the manuscript.
"""

import json
import os
import re
import math
import numpy as np
import fitz  # PyMuPDF

OUT = os.path.dirname(os.path.abspath(__file__))
PT = 72.0 / 25.4                      # PDF points per millimetre

# =====================================================================
# 1.  PK/PD model  (parameters unchanged from gnb_pkpd_v3.py)
# =====================================================================
DRUGS = {
    'cipro': dict(fu=(0.60, 0.70, 0.80),
                  label_auc={500: 27.4, 750: 31.6},
                  label_sd={500: None, 750: None},
                  poppk=dict(CL=25.4, iiv=0.678, F=0.70), fe=0.45),
    'levo': dict(fu=(0.62, 0.69, 0.76),
                 label_auc={500: 47.5, 750: 90.7},
                 label_sd={500: 6.7, 750: 17.6},
                 poppk=dict(CL=8.32, iiv=0.226, F=0.99)),
}
REGIMENS = [('cipro', 500, 2), ('cipro', 750, 2),
            ('levo', 500, 1), ('levo', 750, 1)]

TARGETS = [
    ('fAUC24/MIC >= 72', 'free', 72.0),
    ('AUC24/MIC >= 125', 'total', 125.0),
    ('AUC24/MIC >= 250', 'total', 250.0),
    ('AUC24/MIC >= 500', 'total', 500.0),
]
MIC_GRID = [2.0 ** k for k in range(-5, 4)]

S_BP = {'cipro': 0.25, 'levo': 0.5}    # Enterobacterales susceptible breakpoint

# MIC anchors. The Enterobacterales ECOFF differs by species, so both steps
# are marked: E. coli cipro 0.0625 (reported by EUCAST as 0.06), levo 0.125;
# K. pneumoniae cipro 0.125, levo 0.25.
VMARKS = [(0.0625, "a"), (0.125, "b"), (0.25, "c"), (0.5, "d"), (2.0, "e")]
VKEY = [
    "a   E. coli, ciprofloxacin ECOFF (0.06 mg/L)",
    "b   E. coli, levofloxacin ECOFF and K. pneumoniae, ciprofloxacin ECOFF (0.125 mg/L)",
    "c   K. pneumoniae, levofloxacin ECOFF and ciprofloxacin susceptible breakpoint (0.25 mg/L)",
    "d   P. aeruginosa, ciprofloxacin ECOFF and Enterobacterales levofloxacin susceptible breakpoint (0.5 mg/L)",
    "e   P. aeruginosa, levofloxacin ECOFF (2.0 mg/L)",
]


def _ln_mean_cv(mean, cv):
    s2 = math.log(1.0 + cv * cv)
    return math.log(mean) - 0.5 * s2, math.sqrt(s2)


def _ln_median_logsd(median, logsd):
    return math.log(median), logsd


def scenarios(drug, dose, nday):
    """All (mu, sigma, name) exposure distributions used in the analysis.

    Two parameterisations are in play and must not be confused:
      auc_label() -> lognorm(mean, sd): mean is ARITHMETIC, sd arithmetic.
      auc_poppk() -> exp(N(ln CL_typ, iiv)): CL_typ is the MEDIAN, iiv the LOG-SD.
    """
    d = DRUGS[drug]
    p = d['poppk']
    m = d['label_auc'][dose]
    sd = d['label_sd'][dose]
    cv = (sd / m) if sd is not None else 0.25
    med_auc = p['F'] * dose * nday / p['CL']
    out = [_ln_mean_cv(m, cv) + ('label',)]
    if drug == 'cipro':
        out.append(_ln_mean_cv(m, 0.15) + ('label_cv15',))
        out.append(_ln_mean_cv(m, 0.40) + ('label_cv40',))
        out.append(_ln_median_logsd(med_auc, p['iiv']) + ('poppk_icu',))
        if dose == 750:
            for name, val in (('cand_ss_table', 31.6), ('cand_iv_equiv', 32.9),
                              ('cand_single_x2', 40.4), ('cand_poppk_icu', 41.3)):
                out.append(_ln_mean_cv(val, 0.25) + (name,))
    else:
        out.append(_ln_median_logsd(med_auc, p['iiv']) + ('poppk',))
    return out


def pta(mu, sigma, fu, mic, target, basis):
    m = mu + (math.log(fu) if basis == 'free' else 0.0)
    return 0.5 * (1.0 - math.erf((math.log(target * mic) - m) /
                                  (sigma * math.sqrt(2.0))))


def band(drug, dose, nday, target, basis, mics):
    """(min, max, label-model) PTA envelopes over all (model x fu) scenarios."""
    S = scenarios(drug, dose, nday)
    fus = DRUGS[drug]['fu']
    lo, hi, mid = [], [], []
    mu0, sg0 = _ln_mean_cv(DRUGS[drug]['label_auc'][dose],
                           (DRUGS[drug]['label_sd'][dose] /
                            DRUGS[drug]['label_auc'][dose])
                           if DRUGS[drug]['label_sd'][dose] is not None
                           else 0.25)
    fm = fus[len(fus) // 2]
    for mic in mics:
        ps = [pta(mu, sg, fu, mic, target, basis)
              for (mu, sg, _n) in S for fu in fus]
        lo.append(min(ps) * 100)
        hi.append(max(ps) * 100)
        mid.append(pta(mu0, sg0, fm, mic, target, basis) * 100)
    return np.array(lo), np.array(hi), np.array(mid)


def grid_bp(drug, dose, nday, target, basis, th=0.90):
    """Highest grid MIC with PTA >= th in EVERY scenario (None if none)."""
    b = None
    for mic in MIC_GRID:
        ps = [pta(mu, sg, fu, mic, target, basis)
              for (mu, sg, _n) in scenarios(drug, dose, nday)
              for fu in DRUGS[drug]['fu']]
        if min(ps) >= th:
            b = mic
        else:
            break
    return b


# =====================================================================
# 2.  Minimal vector layer over PyMuPDF   (all coordinates in PDF points)
# =====================================================================
ARIAL = "C:/Windows/Fonts/arial.ttf"
ARIAL_B = "C:/Windows/Fonts/arialbd.ttf"
ARIAL_I = "C:/Windows/Fonts/ariali.ttf"

F_R, F_B, F_I = "ar", "arb", "ari"


def _rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


INK = _rgb("#1A1A1A")
GREY = _rgb("#5A5A5A")
GRIDC = _rgb("#DCDCDC")
WHITE = (1.0, 1.0, 1.0)

BLUE_STEPS = ["#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#084594"]
SBU_STEPS = ["#1A9850", "#FEE08B", "#FC8D59", "#D73027", "#8C0000"]

COL = {('cipro', 500): "#377EB8", ('cipro', 750): "#4DAF4A",
       ('levo', 500): "#FF7F00", ('levo', 750): "#E41A1C"}

MICS = np.exp(np.linspace(math.log(0.03125), math.log(8.0), 400))


def init_fonts(page):
    for nm, fp in ((F_R, ARIAL), (F_B, ARIAL_B), (F_I, ARIAL_I)):
        try:
            page.insert_font(fontname=nm, fontfile=fp)
        except Exception:
            pass


GE, LE, NDASH = chr(0x2265), chr(0x2264), chr(0x2013)
_MISSING = {GE: 0.549, LE: 0.549}      # Arial advance widths, em units


def pretty(s):
    """Typographic substitutions applied to every string before it is set."""
    s = s.replace(">=", GE).replace("<=", LE)
    s = re.sub(r"(?<=\d)-(?=\d)", NDASH, s)
    return s


def textw(s, size, font=F_R):
    """Helvetica and Arial share advance widths, so measure with the base-14 font."""
    s = pretty(s)
    total, buf = 0.0, ""
    try:
        for ch in s:
            if ch in _MISSING:
                if buf:
                    total += fitz.get_text_length(buf, fontname="helv",
                                                  fontsize=size)
                    buf = ""
                total += _MISSING[ch] * size
            else:
                buf += ch
        if buf:
            total += fitz.get_text_length(buf, fontname="helv", fontsize=size)
        return total
    except Exception:
        return 0.5 * size * len(s)


def text(page, x, y, s, size=8.0, font=F_R, col=INK, align="left"):
    if align == "right":
        x -= textw(s, size, font)
    elif align == "center":
        x -= textw(s, size, font) / 2.0
    try:
        page.insert_text(fitz.Point(x, y), pretty(s), fontsize=size,
                         fontname=font, color=col)
    except Exception as e:
        print("  !! text failed:", repr(s), e)


def wrap(s, width, size, font=F_R):
    words, lines, cur = s.split(), [], ""
    for w in words:
        t = (cur + " " + w).strip()
        if cur and textw(t, size, font) > width:
            lines.append(cur)
            cur = w
        else:
            cur = t
    if cur:
        lines.append(cur)
    return lines


def para(page, x, y, s, width, size=5.9, font=F_R, col=GREY, dy=None):
    dy = dy if dy is not None else size * 1.25
    for i, ln in enumerate(wrap(s, width, size, font)):
        text(page, x, y + i * dy, ln, size=size, font=font, col=col)


def vtext(page, x, yc, s, size=8.5, font=F_R, col=INK):
    """Vertically set text, centred on yc."""
    w = textw(s, size, font)
    try:
        page.insert_text(fitz.Point(x, yc + w / 2.0), pretty(s), fontsize=size,
                         fontname=font, color=col, rotate=90)
    except Exception as e:
        print("  !! rotated text failed:", e)


class Ax(object):
    """Plot box: log2 x-axis (MIC), linear y-axis."""

    def __init__(self, x0, y0, w, h, xmin, xmax, ymin, ymax):
        self.x0, self.y0, self.w, self.h = x0, y0, w, h
        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax

    def px(self, v):
        return self.x0 + self.w * (math.log(v, 2) - math.log(self.xmin, 2)) / \
            (math.log(self.xmax, 2) - math.log(self.xmin, 2))

    def py(self, v):
        return self.y0 + self.h * (self.ymax - v) / (self.ymax - self.ymin)

    def seg(self, sh, p1, p2, col=INK, w=0.9, dash=None):
        sh.draw_line(fitz.Point(*p1), fitz.Point(*p2))
        sh.finish(color=col, width=w, dashes=dash, closePath=False)

    def vline(self, sh, mic, col=GREY, w=0.7, dash="[3 3] 0"):
        self.seg(sh, (self.px(mic), self.y0),
                 (self.px(mic), self.y0 + self.h), col, w, dash)

    def hline(self, sh, v, col=GREY, w=0.7, dash="[3 3] 0"):
        self.seg(sh, (self.x0, self.py(v)),
                 (self.x0 + self.w, self.py(v)), col, w, dash)

    def fillband(self, sh, xs, lo, hi, col, op=0.17):
        pts = [(x, self.py(h)) for x, h in zip(xs, hi)]
        pts += [(x, self.py(l)) for x, l in reversed(list(zip(xs, lo)))]
        sh.draw_polyline([fitz.Point(*p) for p in pts])
        sh.finish(color=None, fill=col, width=0, closePath=True,
                  fill_opacity=op)

    def line(self, sh, xs, ys, col, lw=1.4, dash=None):
        sh.draw_polyline([fitz.Point(x, y) for x, y in zip(xs, ys)])
        sh.finish(color=col, width=lw, dashes=dash, closePath=False)

    def frame(self, sh):
        self.seg(sh, (self.x0, self.y0), (self.x0, self.y0 + self.h), INK, 0.9)
        self.seg(sh, (self.x0, self.y0 + self.h),
                 (self.x0 + self.w, self.y0 + self.h), INK, 0.9)

    def axes_ticks(self, page, xticks, yticks, xfmt=lambda v: "%g" % v,
                   yfmt=lambda v: "%g" % v, xlab=None, ylab=None):
        sh = page.new_shape()
        for t in xticks:
            self.seg(sh, (self.px(t), self.y0 + self.h),
                     (self.px(t), self.y0 + self.h + 3), INK, 0.8)
        for t in yticks:
            self.seg(sh, (self.x0, self.py(t)),
                     (self.x0 - 3, self.py(t)), INK, 0.8)
        sh.commit()
        for t in xticks:
            text(page, self.px(t), self.y0 + self.h + 11, xfmt(t),
                 size=7.5, align="center")
        for t in yticks:
            text(page, self.x0 - 5, self.py(t) + 2.6, yfmt(t),
                 size=7.5, align="right")
        if xlab:
            text(page, self.x0 + self.w / 2.0, self.y0 + self.h + 26, xlab,
                 size=8.5, align="center")
        if ylab:
            vtext(page, self.x0 - 30, self.y0 + self.h / 2.0, ylab, 8.5)


def newdoc(w_mm, h_mm):
    doc = fitz.open()
    page = doc.new_page(width=w_mm * PT, height=h_mm * PT)
    init_fonts(page)
    return doc, page


def finish(doc, stem, dpi=600):
    pdf = os.path.join(OUT, stem + ".pdf")
    doc.save(pdf, garbage=4, deflate=True)
    doc.close()
    d = fitz.open(pdf)
    pix = d[0].get_pixmap(dpi=dpi)
    png = os.path.join(OUT, stem + ".png")
    pix.save(png)
    d.close()
    print("  wrote %-14s %10s B   |   %-16s %10s B" %
          (stem + ".pdf", format(os.path.getsize(pdf), ","),
           stem + ".png", format(os.path.getsize(png), ",")))


def label_of(r, short=False):
    drug, dose, nday = r
    name = ("Cipro" if short else "Ciprofloxacin") if drug == 'cipro' \
        else ("Levo" if short else "Levofloxacin")
    return name + " " + str(dose) + " mg " + \
        ("twice daily" if nday == 2 else "once daily")


def legend(page, x, y, entries, dy=11.0, sw=13.0, size=7.0):
    for i, (col, lab, dash) in enumerate(entries):
        yy = y + i * dy
        sh = page.new_shape()
        sh.draw_line(fitz.Point(x, yy - 2.0), fitz.Point(x + sw, yy - 2.0))
        sh.finish(color=_rgb(col) if isinstance(col, str) else col,
                  width=2.1, dashes=dash, closePath=False)
        sh.commit()
        text(page, x + sw + 4, yy, lab, size=size)


# =====================================================================
# 3.  Figures
# =====================================================================
def figure1():
    """PTA vs MIC, four oral regimens, AUC24/MIC >= 125."""
    W, H = 180.0, 108.0
    pw, ph = W * PT, H * PT
    doc, page = newdoc(W, H)

    ax = Ax(48.0, 30.0, pw - 48.0 - 12.0, ph - 30.0 - 74.0,
            0.03125, 8.0, 0.0, 108.0)

    sh = page.new_shape()
    for v in (25, 50, 75, 100):
        ax.seg(sh, (ax.x0, ax.py(v)), (ax.x0 + ax.w, ax.py(v)), GRIDC, 0.5)
    sh.commit()

    sh = page.new_shape()
    for r in REGIMENS:
        lo, hi, mid = band(r[0], r[1], r[2], 125.0, 'total', MICS)
        xs = [ax.px(m) for m in MICS]
        c = _rgb(COL[(r[0], r[1])])
        ax.fillband(sh, xs, lo, hi, c)
        ax.line(sh, xs, [ax.py(v) for v in mid], c, 1.5)
    sh.commit()

    sh = page.new_shape()
    ax.hline(sh, 90.0, _rgb("#333333"), 0.9, "[4 3] 0")
    for v, _t in VMARKS:
        ax.vline(sh, v)
    sh.commit()

    sh = page.new_shape()
    ax.frame(sh)
    sh.commit()

    ax.axes_ticks(page, MIC_GRID, [0, 25, 50, 75, 100],
                  xlab="MIC (mg/L, doubling-dilution scale)",
                  ylab="Probability of target attainment (%)")

    text(page, ax.x0 + ax.w - 3, ax.py(90) - 3.5,
         "90% attainment criterion", size=6.6, col=GREY, align="right")
    for v, t in VMARKS:
        text(page, ax.px(v), ax.y0 - 2.5, t, size=7.5, col=GREY,
             align="center")

    legend(page, ax.px(1.0) + 4, ax.py(70.0),
           [(_rgb(COL[(r[0], r[1])]), label_of(r), None) for r in REGIMENS],
           dy=12.0, size=7.0)

    ybase = ax.y0 + ax.h + 36
    for i, k in enumerate(VKEY):
        text(page, ax.x0, ybase + i * 7.4, k, size=5.9, col=GREY)

    finish(doc, "Figure_1")


def figure2():
    """Breakpoint heat map and SBU interval."""
    W, H = 180.0, 62.0
    pw, ph = W * PT, H * PT
    doc, page = newdoc(W, H)

    gut, cw, ch = 88.0, 48.0, 24.0
    x0, ytop = 8.0, 36.0
    heads = ["fAUC >= 72", "AUC >= 125", "AUC >= 250", "AUC >= 500"]

    for p, (title, sub, shades, is_bp) in enumerate([
            ("(a) PK/PD breakpoint (mg/L)",
             "highest MIC attaining PTA >= 90% in every scenario",
             BLUE_STEPS, True),
            ("(b) SBU interval width (doubling dilutions)",
             "between the PK/PD breakpoint and the susceptible breakpoint",
             SBU_STEPS, False)]):
        gx = x0 + gut + p * 4 * cw
        text(page, gx, 13.0, title, size=8.5, font=F_B)
        text(page, gx, 24.0, sub, size=6.4, col=GREY)
        for j in range(4):
            text(page, gx + j * cw + cw / 2.0, ytop - 5, heads[j],
                 size=6.4, align="center")

    for i, r in enumerate(REGIMENS):
        text(page, x0 + gut - 6, ytop + i * ch + ch / 2.0 + 2.6,
             label_of(r, short=True), size=7.0, align="right")

    for p, (_t, _s, shades, is_bp) in enumerate([
            (None, None, BLUE_STEPS, True),
            (None, None, SBU_STEPS, False)]):
        for i, r in enumerate(REGIMENS):
            for j, (_n, basis, t) in enumerate(TARGETS):
                b = grid_bp(r[0], r[1], r[2], t, basis)
                if is_bp:
                    lab = "none" if b is None else "%g" % b
                    lv = -9.0 if b is None else math.log(b, 2)
                    idx = 0 if lv < -4.5 else int(round(lv + 5))
                    txtc = INK if idx <= 2 else WHITE
                else:
                    if b is None:
                        lab, lv = "> 3", 3.7
                    else:
                        lv = math.log(S_BP[r[0]], 2) - math.log(b, 2)
                        lab = "%d" % int(round(lv))
                    idx = 4 if lv > 3.4 else int(round(lv))
                    txtc = INK if idx <= 1 else WHITE
                idx = max(0, min(4, idx))
                x = x0 + gut + p * 4 * cw + j * cw
                rect = fitz.Rect(x, ytop + i * ch, x + cw, ytop + (i + 1) * ch)
                sh = page.new_shape()
                sh.draw_rect(rect)
                sh.finish(color=WHITE, fill=_rgb(shades[idx]), width=0.4)
                sh.commit()
                text(page, rect.x0 + rect.width / 2.0,
                     rect.y0 + rect.height / 2.0 + 2.6, lab, size=7.2,
                     col=txtc, align="center")

    yb = ytop + 4 * ch + 14
    text(page, x0, yb,
         "SBU (susceptible-but-unattainable) interval = MIC values above the PK/PD breakpoint "
         "and at or below the susceptible breakpoint (ciprofloxacin S <= 0.25, levofloxacin S <= 0.5 mg/L).",
         size=5.9, col=GREY)
    text(page, x0, yb + 9,
         "\"none\" = no MIC on the doubling-dilution grid (>= 0.03125 mg/L) attains the target in every scenario.",
         size=5.9, col=GREY)

    finish(doc, "Figure_2")


def figure3():
    """Levofloxacin PK/PD breakpoint by creatinine clearance."""
    W, H = 140.0, 92.0
    pw, ph = W * PT, H * PT
    doc, page = newdoc(W, H)

    renal = json.load(open(os.path.join(OUT, 'gnb_pkpd_v2_results.json'),
                           encoding='utf-8'))['RENAL']
    strata = ['> 80', '50-80', '20-49', '< 20']
    k500 = ['levo 500 | CLCR ' + s for s in ('>80', '50-80', '20-49', '<20')]
    k750 = ['levo 750 | CLCR ' + s for s in ('>80', '50-80', '20-49', '<20')]

    ax = Ax(52.0, 52.0, pw - 52.0 - 12.0, ph - 52.0 - 50.0,
            0.0, 1.0, 0.0, 1.0)
    lo2, hi2 = math.log(0.0625, 2), math.log(2.0, 2)

    def PY(v):
        return ax.y0 + ax.h * (hi2 - math.log(v, 2)) / (hi2 - lo2)

    gw = ax.w / len(strata)
    bw = gw * 0.24

    sh = page.new_shape()
    for v in MIC_GRID:
        if 0.0625 <= v <= 2.0:
            sh.draw_line(fitz.Point(ax.x0, PY(v)),
                         fitz.Point(ax.x0 + ax.w, PY(v)))
            sh.finish(color=GRIDC, width=0.5, closePath=False)
    sh.commit()

    series = [
        ("Levofloxacin 500 mg once daily, AUC24/MIC >= 125",
         [renal[k]['AUC24/MIC>=125'] for k in k500], "#377EB8"),
        ("Levofloxacin 750 mg once daily, AUC24/MIC >= 125",
         [renal[k]['AUC24/MIC>=125'] for k in k750], "#E41A1C"),
        ("Levofloxacin 750 mg once daily, AUC24/MIC >= 250",
         [renal[k]['AUC24/MIC>=250'] for k in k750], "#984EA3"),
    ]
    for k, (_lab, vals, col) in enumerate(series):
        c = _rgb(col)
        for i, v in enumerate(vals):
            x = ax.x0 + i * gw + gw / 2.0 + (k - 1) * bw
            r = fitz.Rect(x - bw / 2.0, PY(v), x + bw / 2.0, PY(0.0625))
            sh = page.new_shape()
            sh.draw_rect(r)
            sh.finish(color=None, fill=c, width=0)
            sh.commit()
            text(page, x, PY(v) - 2.6, "%g" % v, size=6.2, align="center")

    sh = page.new_shape()
    sh.draw_line(fitz.Point(ax.x0, PY(0.0625)),
               fitz.Point(ax.x0 + ax.w, PY(0.0625)))
    sh.finish(color=INK, width=0.9, closePath=False)
    ax.seg(sh, (ax.x0, ax.y0), (ax.x0, ax.y0 + ax.h), INK, 0.9)
    sh.commit()

    # Enterobacterales levofloxacin ECOFF reference lines (species-specific)
    sh = page.new_shape()
    ax.seg(sh, (ax.x0, PY(0.25)), (ax.x0 + ax.w, PY(0.25)),
           _rgb("#333333"), 0.9, "[4 3] 0")
    ax.seg(sh, (ax.x0, PY(0.125)), (ax.x0 + ax.w, PY(0.125)),
           _rgb("#333333"), 0.9, "[4 3] 0")
    sh.commit()
    # place the annotations in the empty top-left interior so they never cross bars
    _lx, _ly = ax.x0 + 4, PY(1.45)
    sh = page.new_shape()
    ax.seg(sh, (_lx, _ly - 1.8), (_lx + 14, _ly - 1.8),
           _rgb("#333333"), 0.9, "[4 3] 0")
    sh.commit()
    text(page, _lx + 17, _ly,
         "K. pneumoniae levofloxacin ECOFF (0.25 mg/L)",
         size=6.0, col=GREY, align="left")
    _ly2 = _ly + 7.0
    sh = page.new_shape()
    ax.seg(sh, (_lx, _ly2 - 1.8), (_lx + 14, _ly2 - 1.8),
           _rgb("#333333"), 0.9, "[4 3] 0")
    sh.commit()
    text(page, _lx + 17, _ly2,
         "E. coli levofloxacin ECOFF (0.125 mg/L)",
         size=6.0, col=GREY, align="left")

    sh = page.new_shape()
    for v in MIC_GRID:
        if 0.0625 <= v <= 2.0:
            ax.seg(sh, (ax.x0, PY(v)), (ax.x0 - 3, PY(v)), INK, 0.8)
    sh.commit()
    for v in MIC_GRID:
        if 0.0625 <= v <= 2.0:
            text(page, ax.x0 - 5, PY(v) + 2.6, "%g" % v, size=7.5,
                 align="right")

    for i, s in enumerate(strata):
        text(page, ax.x0 + i * gw + gw / 2.0, PY(0.0625) + 11, s,
             size=7.5, align="center")

    text(page, ax.x0 + ax.w / 2.0, PY(0.0625) + 27,
         "Creatinine clearance (mL/min)", size=8.5, align="center")
    vtext(page, ax.x0 - 30, (ax.y0 + PY(0.0625)) / 2.0,
          "PK/PD breakpoint (mg/L)", 8.5)

    legend(page, ax.x0, 14.0,
           [(_rgb(col), lab, None) for lab, _v, col in series],
           dy=11.0, size=6.6)

    finish(doc, "Figure_3")


def figure4():
    """(a) P. aeruginosa PTA vs MIC; (b) breakpoint system x target matrix."""
    W, H = 180.0, 92.0
    pw, ph = W * PT, H * PT
    doc, page = newdoc(W, H)

    axL = Ax(48.0, 34.0, 236.0, 176.0, 0.03125, 8.0, 0.0, 112.0)

    sh = page.new_shape()
    for v in (25, 50, 75, 100):
        axL.seg(sh, (axL.x0, axL.py(v)), (axL.x0 + axL.w, axL.py(v)),
              GRIDC, 0.5)
    sh.commit()

    styles = [(72.0, 'free', "#E41A1C", None),
              (125.0, 'total', "#984EA3", "[5 3] 0"),
              (250.0, 'total', "#377EB8", "[2 2] 0")]
    sh = page.new_shape()
    for tgt, basis, col, dash in styles:
        lo, hi, mid = band('levo', 750, 1, tgt, basis, MICS)
        xs = [axL.px(m) for m in MICS]
        c = _rgb(col)
        axL.fillband(sh, xs, lo, hi, c)
        axL.line(sh, xs, [axL.py(v) for v in mid], c, 1.5, dash)
    sh.commit()

    sh = page.new_shape()
    axL.hline(sh, 90.0, _rgb("#333333"), 0.9, "[4 3] 0")
    for v in (1.0, 2.0):
        axL.vline(sh, v)
    sh.commit()

    sh = page.new_shape()
    axL.frame(sh)
    sh.commit()
    axL.axes_ticks(page, MIC_GRID, [0, 25, 50, 75, 100],
                   xlab="MIC (mg/L, doubling-dilution scale)",
                   ylab="Probability of target attainment (%)")

    text(page, axL.x0, axL.y0 - 13, "Pseudomonas aeruginosa", size=7.0,
         col=GREY, font=F_I)
    text(page, axL.px(1.0) - 2.5, axL.y0 - 3, "CLSI S", size=6.2, col=GREY,
         align="right")
    text(page, axL.px(2.0) + 2.5, axL.y0 - 3, "EUCAST ECOFF", size=6.2,
         col=GREY, align="left")

    legend(page, axL.px(0.03125) + 4, axL.py(46.0),
           [(_rgb(col), ("f" if basis == 'free' else "") +
             "AUC/MIC >= %g" % tgt, dash)
            for tgt, basis, col, dash in styles],
           dy=11.0, size=6.8)

    text(page, axL.x0, axL.y0 + axL.h - 2, "(a)", size=9.0, font=F_B)

    # -------------------- (b) decision matrix
    bx, gut, cw, ch, ytop = 302.0, 62.0, 46.0, 46.0, 78.0
    text(page, bx, 30.0, "(b) Levofloxacin 750 mg once daily", size=8.5,
         font=F_B)
    para(page, bx, 42.0,
         "green = eligible in every scenario; red = not eligible in any; "
         "grey = AHA precondition not met",
         206.0, size=5.9)

    cols = ["fAUC >= 72", "AUC >= 125", "AUC >= 250"]
    rows = ["EUCAST (ECOFF 2.0 mg/L)", "CLSI (S <= 1.0 mg/L)"]
    cells = [[(0.0, 0.0, 2), (0.0, 0.0, 2), (0.0, 0.0, 2)],
             [(8.4, 39.5, 0), (3.9, 6.8, 0), (0.0, 0.0, 0)]]

    for j, c in enumerate(cols):
        text(page, bx + gut + j * cw + cw / 2.0, ytop - 5, c, size=6.4,
             align="center")
    for i, lab in enumerate(rows):
        text(page, bx + gut - 6, ytop + i * ch + ch / 2.0 + 2.6, lab,
             size=6.8, align="right")
        for j in range(3):
            lo, hi, v = cells[i][j]
            rect = fitz.Rect(bx + gut + j * cw, ytop + i * ch,
                             bx + gut + (j + 1) * cw, ytop + (i + 1) * ch)
            col = _rgb(["#D73027", "#1A9850", "#EDEDED"][v])
            sh = page.new_shape()
            sh.draw_rect(rect)
            sh.finish(color=WHITE, fill=col, width=0.4)
            sh.commit()
            lab = "%g-%g%%" % (lo, hi)
            text(page, rect.x0 + rect.width / 2.0,
                 rect.y0 + rect.height / 2.0 + 2.6, lab, size=6.6,
                 col=(INK if v == 2 else WHITE), align="center")

    para(page, bx, ytop + 2 * ch + 14,
         "No cell attains 90%. Under EUCAST no susceptible category exists "
         "(S <= 0.001 mg/L; Note 13: never report S), so the AHA precondition "
         "of susceptibility-based selection cannot be met. Under CLSI a "
         "breakpoint exists (levofloxacin S <= 1 mg/L), but the listed "
         "regimen does not attain any target at that MIC.",
         206.0, size=5.9)

    finish(doc, "Figure_4")


def main():
    print("generating figures ...")
    figure1()
    figure2()
    figure3()
    figure4()
    print("done.")


if __name__ == '__main__':
    main()
