# -*- coding: utf-8 -*-
"""Builds manuscript/Nontraversable_ER_Bridge_Cyclic_Region.docx from phase2e_paper_content.py.   python phase2e_paper_build.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import docx_engine as E
import paper_build as PB
import phase2e_paper_content as PC

OUT = os.path.join(PB.OUT, "Nontraversable_ER_Bridge_Cyclic_Region.docx")
if __name__ == "__main__":
    os.makedirs(PB.OUT, exist_ok=True)
    cit = E.Citations(PC.LIB)
    for b in PC.B:
        for x in b[1:]:
            if isinstance(x, str):
                cit.scan(x)
            elif isinstance(x, list):
                for row in x:
                    for c in (row if isinstance(row, list) else [row]):
                        if isinstance(c, str):
                            cit.scan(c)
    PB.build(OUT, PC.B, cit, journal=True)
    abstract = [b[1] for b in PC.B if b[0] == "p0"][0]
    print("built %s; %d references; abstract %d words" % (OUT, len(PC.LIB), len(abstract.split())))
    assert len(abstract.split()) <= 250, "abstract too long (%d)" % len(abstract.split())
