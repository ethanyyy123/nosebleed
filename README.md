# nosebleed -- computational study of fibrin "knob A" peptide variants at fibrinogen hole "a"

Open, CPU-only pipeline and manuscript for an in silico study motivated by
fibrin-targeted topical hemostats (e.g. for nosebleeds): can an inexpensive
ensemble MM/GBSA protocol rank single-residue variants of the fibrin knob-A
model peptide **Gly-Pro-Arg-Pro (GPRP)** for binding to hole "a" of the
fibrinogen gamma-module, and is the protocol trustworthy enough to say so?

> **Evidentiary boundary.** Everything here is *computed*. No peptide was made
> or tested. Predicted energies are model outputs with unquantified systematic
> error; they are hypotheses for experiment (SPR/ITC, clot-turbidity), not
> evidence of improved binding or hemostatic activity.

## Why this idea (of the three considered)
See `docs/01_journal_venue_requirements.md` and the manuscript. In short:
the botanical-docking idea mis-targets the plasmin active site (tranexamic acid
acts at kringle lysine-binding sites, millimolar affinity); the RADA16
pH/ion idea needs assembly timescales (and coarse-grained tooling) beyond four
CPU cores, and simulations of RADA16 variants in salt already exist; the knob-mimic idea has real crystal structures,
published binding data to check against, and a gap in the literature (no
published free-energy-type ranking of designed knob mimics was found).

## Layout
| Path | Content |
|---|---|
| `docs/01_journal_venue_requirements.md` | target venue, submission rules, provenance warnings |
| `docs/02_study_design.md` | hypotheses and decision rules, frozen before production runs |
| `docs/03_reference_verification.md` | what was and was not verified for each citation |
| `scripts/` | pipeline (`knobhole_lib.py`, `run_jobs.py`, `make_jobs.py`, `aggregate.py`, analysis and figure scripts) |
| `data/structures/` | PDB entries used (fetched from the wwPDB AWS mirror); `structure_survey.tsv` lists 29 knob-hole complexes |
| `data/results/` | raw per-job JSON, pilot data, summary tables |
| `data/literature_values.csv` | external numbers used, with provenance flags |
| `figures/` | publication figures (PNG/PDF) |
| `manuscript/` | manuscript source (Markdown), DOCX in JEI order, PDF |

## Reproduce
```
pip install -r requirements.txt          # + apt: texlive, pandoc
python scripts/00_fetch_structures.py 2FIB   # from the wwPDB snapshot bucket
python scripts/make_jobs.py
python scripts/run_jobs.py data/jobs/phaseA.json --workers 4
python scripts/run_jobs.py data/jobs/phaseB.json --workers 4
```
Throughput on 4 CPU cores: ~7.8 min per 214-ps job (6 A shell), 4 jobs in parallel.
