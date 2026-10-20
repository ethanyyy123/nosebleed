# nosebleed -- computational study of fibrin "knob A" peptide variants at fibrinogen hole "a"

Open, CPU-only pipeline and manuscript for an in silico study motivated by
fibrin-targeted topical hemostats (e.g. for nosebleeds): can an inexpensive
ensemble MM/GBSA protocol rank single-residue variants of the fibrin knob-A
model peptide **Gly-Pro-Arg-Pro (GPRP)** for binding to hole "a" of the
fibrinogen gamma-module, and is the protocol trustworthy enough to say so?

**Short answer: not fully, and the paper is about why.** The pipeline
correctly and strongly reproduces the known importance of the pocket's
Asp364, but misses a statistical-power threshold on one pre-registered check
and fails a specificity control on another -- so the 76-variant saturation
scan is reported as **unvalidated**: it turns out to be dominated by a known
MM/GBSA artifact (predicted affinity correlates with substituted side-chain
size at rho = -0.79, p = 3x10^-17), not real knob-hole chemistry. See
**[`manuscript/draft.md`](manuscript/draft.md)** /
**[`manuscript/submission/manuscript.pdf`](manuscript/submission/manuscript.pdf)**
for the full paper.

> **Evidentiary boundary.** Everything here is *computed*. No peptide was made
> or tested. Predicted energies are model outputs with unquantified systematic
> error; they are hypotheses for experiment (SPR/ITC, clot-turbidity), not
> evidence of improved binding or hemostatic activity. No candidate from the
> saturation scan is a validated or confirmed tighter binder.

## Start here
- **[`manuscript/draft.md`](manuscript/draft.md)** -- the manuscript, with a
  status header listing exactly what an author still needs to add (name,
  affiliation, mentor) before submission.
- **[`manuscript/submission/manuscript.pdf`](manuscript/submission/manuscript.pdf)**
  -- a compiled, readable PDF of the same content (a convenience rendering;
  Reinvention's actual submission format is Word/RTF/OpenOffice, Methods/docs/01).
- **[`docs/02_study_design.md`](docs/02_study_design.md)** -- the hypotheses
  and numeric decision rules, frozen in writing *before* production runs, plus
  a dated deviation log (every departure from the frozen plan, and why).

## Why this idea (of the three considered)
See `docs/01_journal_venue_requirements.md` and the manuscript Discussion. In
short: the botanical-docking idea mis-targets the plasmin active site
(tranexamic acid acts at kringle lysine-binding sites, millimolar affinity);
the RADA16 pH/ion idea needs assembly timescales and coarse-grained tooling
beyond four CPU cores, and simulations of RADA16 variants in salt already
exist in the literature; the knob-mimic idea had real crystal structures,
published binding data to check against, and a gap in the literature (no
published free-energy-type ranking of designed knob mimics was found) -- which
is also what made it possible to discover the pipeline's own limitations.

## Layout
| Path | Content |
|---|---|
| `docs/01_journal_venue_requirements.md` | target venue (Reinvention), submission rules, provenance warnings |
| `docs/02_study_design.md` | hypotheses and decision rules, frozen before production runs; deviation log |
| `docs/03_reference_verification.md` | what was and was not verified for each citation |
| `scripts/` | pipeline (`knobhole_lib.py`, `run_jobs.py`, `make_jobs.py`, `aggregate.py`, `06_analyze.py`, figure scripts `08`-`10`) |
| `data/structures/` | PDB entries used (fetched from the wwPDB AWS mirror); `structure_survey.tsv` lists 29 knob-hole complexes |
| `data/jobs/` | the frozen Phase A (validation) and Phase B (saturation scan) job lists |
| `data/results/` | raw per-replicate JSON (`raw/`), trajectories (`../trajectories/`), summary CSVs, H1/scan/diagnostic outputs |
| `data/literature_values.csv` | external numbers used, with provenance flags |
| `figures/` | the three publication figures (PNG + PDF) |
| `manuscript/` | manuscript source (`draft.md`), LaTeX preamble, compiled `submission/manuscript.pdf` |

## Reproduce
```
pip install -r requirements.txt          # + apt: autodock-vina openbabel texlive-* pandoc pymol
python scripts/00_fetch_structures.py 2FIB 2ffd 3bvh 1fzc 2hpc 3fib   # wwPDB snapshot bucket
python scripts/make_jobs.py
python scripts/run_jobs.py data/jobs/phaseA.json --workers 4   # validation panel (H1), 29 replicates
python scripts/run_jobs.py data/jobs/phaseB.json --workers 4   # saturation scan (H2), 152 replicates
python scripts/06_analyze.py h1      # validation panel table + pass/fail
python scripts/06_analyze.py scan    # scan table, size-bias diagnostic, stage-1 candidates
python scripts/09_fig2_h1panel.py && python scripts/10_fig3_scan.py && python scripts/08b_fig1.py
```
Throughput on 4 CPU cores: ~7-13 min per 200-ps job depending on pocket-shell
size; the full campaign (181 replicates) took on the order of 10 CPU-hours.
