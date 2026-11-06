# -*- coding: utf-8 -*-
"""Shared helpers for the manuscript builders: results (read from the JSON files the scripts write), number formatting, block lists."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
FIG = os.path.join(ROOT, "figures")
OUT = os.path.join(ROOT, "manuscript")
TW = json.load(open(os.path.join(HERE, "thinwall.json"), encoding="utf-8"))
VAL = json.load(open(os.path.join(HERE, "validation.json"), encoding="utf-8"))
SCAN = json.load(open(os.path.join(HERE, "collapse_scan.json"), encoding="utf-8"))
SCANQ = json.load(open(os.path.join(HERE, "collapse_scan_q.json"), encoding="utf-8"))
VE = json.load(open(os.path.join(HERE, "validation_extra.json"), encoding="utf-8"))
DOIS = json.load(open(os.path.join(HERE, "dois.json"), encoding="utf-8"))
GITHUB = "https://github.com/sandlerleon/vacuum-transition-throats"
VERSION = "1.1.0"
TITLE = "Wormhole Formation in Matter-Coupled Vacuum Transitions: Gravitational Collapse, Horizon Constraints, and the Limits of Cosmological Rebirth"


def n(x, d=0):
    s = "{:,.{d}f}".format(x, d=d)
    if set(s) <= set("-0.,"):
        return s.lstrip("-")
    return s.replace("-", "−")


def sci(x, d=1):
    """scientific notation as '2.7 × 10^{-30}' for the inline-math engine"""
    m, e = ("%.*e" % (d, x)).split("e")
    return "%s × 10$^{%d}$" % (m, int(e))


def pc(x, d=0):
    return n(100.0 * x, d) + "%"


B = []      # manuscript blocks
BS = []     # supplementary blocks


def add(b):
    B.append(b)


def adds(b):
    BS.append(b)
