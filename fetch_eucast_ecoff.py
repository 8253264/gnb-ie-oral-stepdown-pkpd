#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fetch EUCAST epidemiological cut-off (ECOFF) values together with the number of
observations and the number of contributing distributions (sources), directly
from the EUCAST MIC distribution database.

    https://mic.eucast.org/search/

Why this script exists
----------------------
The MIC distribution search page renders its result table server-side, but the
form accepts only ONE filter at a time: changing the Antimicrobial dropdown
resets Species to -1, and changing Species resets Antimicrobial to -1. Sending
both parameters at once returns an empty result container. The browser
therefore shows, for a chosen species, the full table of every agent, and the
ECOFF of interest is read off the corresponding row.

This script reproduces that request and then parses the row, so the numbers in
the manuscript can be re-derived at any time without a browser.

Usage
-----
    python fetch_eucast_ecoff.py                  # fetch all six combinations
    python fetch_eucast_ecoff.py --json out.json  # also write machine-readable output
    python fetch_eucast_ecoff.py --species 261    # restrict to one species id

Only the Python standard library is required.
"""

import argparse
import http.cookiejar as cookiejar
import io
import json
import re
import ssl
import sys
import time
import urllib.parse
import urllib.request

ORIGIN = "https://mic.eucast.org"
SEARCH = ORIGIN + "/search/"
UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/124.0 Safari/537.36")

# Species ids are the option values of the `search[species]` select on the page.
SPECIES = {
    "Escherichia coli": 261,
    "Klebsiella pneumoniae": 288,
    "Pseudomonas aeruginosa": 411,
}

# Antibiotic ids are the option values of the `search[antibiotic]` select.
AGENTS = {
    "Ciprofloxacin": 60,
    "Levofloxacin": 117,
}

TARGETS = ("Ciprofloxacin", "Levofloxacin")

# Column layout of the result table:
#   [0] agent name | [1..19] counts per doubling dilution | [20] Distributions
#   | [21] Observations | [22] (T)ECOFF | [23] Confidence interval
N_COL = 24
COL_DIST, COL_OBS, COL_ECOFF, COL_CI = 20, 21, 22, 23


def _clean(cell_html):
    t = re.sub(r"<[^>]+>", " ", cell_html)
    for a, b in (("&nbsp;", " "), ("&#x20;", " "), ("&amp;", "&"),
                 ("&#8804;", "<="), ("&le;", "<="), ("&#x2264;", "<="),
                 ("&#8805;", ">=")):
        t = t.replace(a, b)
    return re.sub(r"\s+", " ", t).strip()


def build_opener():
    ctx = ssl.create_default_context()
    jar = cookiejar.CookieJar()
    op = urllib.request.build_opener(
        urllib.request.HTTPCookieProcessor(jar),
        urllib.request.HTTPSHandler(context=ctx))
    op.addheaders = [("User-Agent", UA), ("Accept", "text/html,*/*"),
                     ("Referer", SEARCH)]
    return op


def page_url(species_id, page=None, limit=50):
    kv = [("search[method]", "mic"),
          ("search[species]", str(species_id)),
          ("search[limit]", str(limit))]
    if page and page > 1:
        kv.append(("page", str(page)))
    return SEARCH + "?" + urllib.parse.urlencode(kv)


def parse_rows(html):
    i = html.find('id="table-obs"')
    if i < 0 or "<table" not in html[i:]:
        return []
    seg = html[i:]
    rows = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", seg, re.S):
        cells = [_clean(c) for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", tr, re.S)]
        if len(cells) >= N_COL and cells[0] and re.match(r"^[A-Z]", cells[0]):
            rows.append(cells)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", dest="json_out", help="write results to this JSON file")
    ap.add_argument("--species", type=int, help="restrict to one species id")
    ap.add_argument("--html-dir", help="also archive raw HTML of the first page per species")
    args = ap.parse_args()

    op = build_opener()
    op.open(SEARCH, timeout=120).read()      # establish the PHP session first

    wanted = {k: v for k, v in SPECIES.items()
              if args.species is None or v == args.species}
    out = []
    for sp, sid in wanted.items():
        page, found = 1, {}
        while page <= 12:
            url = page_url(sid, page)
            html = op.open(url, timeout=150).read().decode("utf-8", "replace")
            if page == 1 and args.html_dir:
                io.open("%s/EUCAST_MICdistribution_%s_full.html"
                        % (args.html_dir, sp.replace(" ", "_")),
                        "w", encoding="utf-8").write(html)
            rows = parse_rows(html)
            if not rows:
                break
            for c in rows:
                if c[0] in TARGETS and c[0] not in found:
                    found[c[0]] = c
            if len(found) == len(TARGETS) or len(rows) < 50:
                break
            page += 1
            time.sleep(0.5)

        for t in TARGETS:
            if t not in found:
                print("  !! %s / %s : row not found" % (sp, t), file=sys.stderr)
                continue
            c = found[t]
            rec = {
                "species": sp,
                "agent": c[0],
                "distributions": c[COL_DIST],
                "observations": c[COL_OBS],
                "ecoff_mg_L": c[COL_ECOFF],
                "confidence_interval": c[COL_CI],
            }
            out.append(rec)
            print("%-24s %-14s ECOFF WT <= %-6s  obs=%-7s dist=%-4s CI=%s"
                  % (sp, t, rec["ecoff_mg_L"], rec["observations"],
                     rec["distributions"], rec["confidence_interval"]))

    if args.json_out:
        io.open(args.json_out, "w", encoding="utf-8").write(
            json.dumps(out, ensure_ascii=False, indent=2))
        print("\nwrote %s" % args.json_out)


if __name__ == "__main__":
    main()
