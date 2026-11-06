# -*- coding: utf-8 -*-
"""Cover letter for the submission to Classical and Quantum Gravity (IOP Publishing); numbers from the JSON result files, DOIs from dois.json.

    python build_cover_letter.py   -> manuscript/CQG_Cover_Letter.docx
"""
import datetime
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from docx import Document
from docx.shared import Pt, Inches

from pcommon import *

M1 = [r for r in TW["mass_table"] if r["q"] == 1.0][0]
SC = VAL["scan"]
d = Document()
st = d.styles["Normal"]; st.font.name = "Calibri"; st.font.size = Pt(11); st.paragraph_format.space_after = Pt(8)
for s in d.sections:
    s.left_margin = s.right_margin = Inches(1.0); s.top_margin = s.bottom_margin = Inches(0.9)


def para(t, bold=False, after=None):
    p = d.add_paragraph(); r = p.add_run(t); r.bold = bold
    if after is not None:
        p.paragraph_format.space_after = Pt(after)
    return p


para("Leon Sandler", True, 0)
for t in ("Independent Researcher", "Northbrook, Illinois, USA", "leon.sandler@alumni.swinburne.edu", "ORCID 0009-0007-4584-808X"):
    para(t, after=0)
para("")
para(datetime.date.today().strftime("%d %B %Y").lstrip("0"))
para("The Editors", after=0); para("Classical and Quantum Gravity", after=0); para("IOP Publishing")
para("Dear Editors,")
para("I submit the manuscript “%s” as a Paper for consideration in Classical and Quantum Gravity." % TITLE)
para("Proposals for cyclic and rebirth cosmologies often rely on a spacetime throat that forms near the end of the universe and carries matter into a new expanding region. The manuscript asks whether a "
     "specific, density-triggered vacuum transition can produce such a throat, and answers by calculation. The model is a canonical scalar field coupled to dark matter, whose false vacuum carries the "
     "observed dark-energy density and whose true vacuum lies below it, so that the critical sizes and masses are numbers fixed by the observed dark-energy scale and not free parameters.")
resolved = sum(1 for r in VE["inconclusive"] if any(v["first_trapped"] is not None for k, v in r.items() if k.startswith("refined")))
CRE = VE["creation"]
para("The manuscript makes five contributions. (i) A short proof that the whole model satisfies the null energy condition, so that no null ray can have a minimum of the areal radius and no traversable "
     "throat can form; the question is narrowed to a non-traversable bridge. (ii) A finite-tension thin-wall criterion from the Israel junction conditions: a critical exterior mass m_c(sigma), "
     "equal to (1+q) R_dS up to corrections of order sigma^2 (q the depth of the true vacuum, R_dS the de Sitter radius), below which every regular false-vacuum configuration is confined and collapses, "
     "and the observation, stated as a geometric coincidence of the small-tension limit and not as a thermodynamic criterion, that at that limit the exterior horizon and the interior de Sitter horizon have "
     "equal area. (iii) A spherically symmetric Einstein–scalar code in double-null coordinates with a regular centre, validated against exact de Sitter and anti-de Sitter solutions, whose constraints and "
     "null expansions converge on the dynamical solutions, and which forms a trapped surface in %d of %d runs. (iv) A separation of persistence from creation: an expanding core coexists with the trapped "
     "region only when the initial data already contain a false-vacuum core larger than its own de Sitter radius, and in none of the %d runs that start with a smaller core does one appear; the simulations "
     "therefore do not show the dynamical creation of a daughter universe. (v) An accounting at the observed dark-energy scale: the required region has a radius of %s Gpc and a "
     "vacuum-energy deficit of %s solar masses, %s times the matter inside a sphere of that radius today. The conclusions are also tested for a Minkowski and a lower de Sitter exterior."
     % (SC["n_trapped"] + resolved, SC["n_runs"], CRE["n_initially_normal_core"], n(TW["R_dS_Gpc"], 1), "%s × 10^23" % n(M1["M_c_Msun"] / 1e23, 1),
        n(M1["deficit_over_matter_in_RdS_sphere"], 1)))
para("The main finding is a conditional negative, stated as such. Localized collapse does not produce a daughter region, matter is not transferred to it, and the expanding region that does exist in the model is a "
     "pre-existing, horizon-sized false-vacuum core that the collapse of its surroundings isolates; the mechanism is therefore not supported as a cosmological reset, and even a created daughter would be "
     "cosmological reproduction, not a cycle. The manuscript is explicit about what the numerical work is and is not: the walls it can resolve are far thicker and heavier than the physical wall "
     "(by thirty and sixty orders of magnitude), so the experiments test the structure of the thin-wall solution and not the physical values; the evolution is spherical and classical; the dark matter that "
     "motivates the transition is not simulated; no spacelike slice through the bridge is constructed numerically; and the code does not follow the solution beyond the singularity. "
     "It lists six criteria by which each of its statements can be falsified.")
para("The geometry of false-vacuum bubbles, baby universes and the obstacle to creating them from regular data is classical, and the manuscript credits it (Blau, Guendelman and Guth; Farhi and Guth; "
     "Farhi, Guth and Guven; Fischler, Morgan and Polchinski). The contribution is the application to a model with a fixed vacuum scale, the exact critical mass and its entropy interpretation, "
     "and the numerical test with a regular centre.")
para("The code and results are openly available at %s and archived at https://doi.org/%s (code, all versions); the manuscript and its Supplementary Information are archived at https://doi.org/%s "
     "(all versions). Every number in the text is read from the result files that the scripts write, fifteen automated checks run on the stored results, and one command per script regenerates each figure and table."
     % (GITHUB, DOIS["code_concept"], DOIS["paper_concept"]))
para("Author affiliation and email address. I am an independent researcher with no current institutional affiliation, and I am the sole author. My contact and corresponding-author address is the alumni "
     "email address provided by my former university (Swinburne University of Technology, leon.sandler@alumni.swinburne.edu); it identifies me as a graduate and does not indicate a current "
     "affiliation or employer.")
para("A generative AI tool (Claude, Anthropic) assisted with code development, numerical experiments, drafting and editing; I reviewed and checked all of it and take full responsibility for the content, as "
     "disclosed in the manuscript. The manuscript is original, has not been published elsewhere and is not under consideration at any other journal. I am the sole author, the work received no external funding "
     "and I have no competing interests.")
para("Thank you for considering this submission.")
para("Sincerely,")
para("Leon Sandler", True)
out = os.path.join(OUT, "CQG_Cover_Letter.docx")
d.save(out)
print("saved", out)
