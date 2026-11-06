# -*- coding: utf-8 -*-
"""Fetches title/authors/journal/volume/pages/year for candidate DOIs from Crossref to build refs.json (only verified entries are cited)."""
import json, sys, urllib.request
DOIS = """10.1007/BF02710419 10.1103/PhysRevD.35.1747 10.1016/0370-2693(87)90786-4 10.1016/0550-3213(90)90357-Y 10.1103/PhysRevD.41.2638 10.1103/PhysRevD.21.3305
10.1119/1.15620 10.1103/PhysRevD.56.4745 10.1142/S0218271899000286 10.1103/PhysRevD.49.6467 10.1103/PhysRev.136.B571 10.1103/PhysRevD.51.5558 10.12942/lrr-2007-5
10.1007/BF01205930 10.1103/PhysRevLett.14.57 10.1103/PhysRevD.7.2333 10.1103/PhysRevD.15.2738 10.1126/science.1070462 10.1088/0264-9381/9/1/016
10.1016/0370-2693(89)91114-3 10.1051/0004-6361/201833910 10.1103/PhysRevLett.104.231301 10.1016/0550-3213(85)90304-0 10.1103/PhysRevD.36.2919 10.1103/RevModPhys.74.825
10.1103/PhysRevLett.70.9 10.1016/0370-2693(87)91336-9 10.1016/0370-2693(87)90786-4""".split()
out = {}
for d in dict.fromkeys(DOIS):
    try:
        r = json.load(urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + d, headers={"User-Agent": "refcheck (mailto:leon.sandler@alumni.swinburne.edu)"}), timeout=30))["message"]
        au = ", ".join("%s %s" % ((a.get("given") or "")[:1] + ".", a.get("family", "")) for a in r.get("author", [])[:6])
        out[d] = dict(title=r["title"][0], authors=au, journal=(r.get("container-title") or [""])[0], volume=r.get("volume"), page=r.get("page") or r.get("article-number"), year=(r.get("issued", {}).get("date-parts") or [[None]])[0][0])
        print(d, "|", au, "|", out[d]["title"][:80], "|", out[d]["journal"], out[d]["volume"], out[d]["page"], out[d]["year"])
    except Exception as e:
        print(d, "FAILED", e)
json.dump(out, open("refs_fetched.json", "w"), indent=1)
