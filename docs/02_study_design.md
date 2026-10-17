# Study design and analysis plan (frozen before production runs)

**Status:** written after the pilot (`data/results/pilot*`) and before any
production job was started. Pilot data were used only to estimate run time,
replicate noise and equilibration time of the *reference peptide GPRP* (plus
three arbitrary non-benchmark variants, GPKP/GPWP/GADP, used to check that the
build procedure works for substitutions). **No result for any benchmark
comparison (GPR, GPRV, GHRP, Asp364 mutants, control mutants) had been
examined when these rules were fixed.** Any later change is listed in the
deviation log at the end of this file.

## 1. Question

Fibrin assembly depends on "knob A" (the new N-terminus Gly-Pro-Arg-... that
thrombin exposes) binding "hole a" in the fibrinogen gamma-module. Can a cheap,
CPU-only ensemble MM/GBSA protocol (a) reproduce documented determinants of
that interaction and, if so, (b) suggest single-residue changes to the model
knob GPRP that are predicted to bind hole a more tightly?

## 2. Model and protocol (frozen)

| Item | Setting |
|---|---|
| Template | PDB 2FIB (human gamma-module 143-411 + GPRP, 2.01 A, wwPDB snapshot 2024-01-01) |
| Pocket model | all residues with a heavy atom within **6 A** of the template GPRP (+ whole Ca2+ site when any of its ligands is inside), segments capped with ACE/NME built from the real neighbouring atoms |
| Peptide | N-term free (NH3+), C-term N-methylamide; positions P1-P4 = G-P-R-P |
| Force field / solvent | AMBER ff14SB; Ca2+ from the AMBER/OpenMM ion set (Li-Merz); OBC-II generalized Born (OpenMM `GBSAOBCForce`), eps_in = 1, eps_out = 78.5, no salt term |
| Restraints | all receptor C-alpha atoms and Ca2+ harmonic 1 kcal/mol/A^2; flat-bottom peptide "leash" (6 A, backbone centroid) |
| Dynamics | 4-fs Langevin (300 K, 1/ps), H mass 3 amu, H-bond constraints; minimise -> 7 x 2 ps anneal (150-250-350-450-400-350-300 K) -> **100 ps equilibration (discarded)** -> **100 ps production, 20 frames** |
| Score | single-trajectory MM/GBSA, dG = G(complex) - G(receptor) - G(peptide) on the same frames, no entropy term; reported in kcal/mol, more negative = tighter |
| Unit of analysis | one independent MD replicate (mean over its 20 frames); frames are not treated as independent |
| Replicates | reference GPRP n = 6; benchmark mutants n = 3; scan variants n = 2 (n = 6 for candidates, see Phase C) |

Why **6 A**: a full-protein model (3,932 atoms) ran at 3.7 ns/day on the 4
available CPU cores versus ~39 ns/day/core for the 6 A pocket (619 atoms).
Shell-size sensitivity is checked (Phase C, 8 A).

## 3. Hypotheses and decision rules

Notation: dG = replicate-level MM/GBSA value; ddG(X) = mean dG(X) - mean dG(reference).
Negative ddG = predicted tighter binding. Welch t-test on replicate means.
A practical-significance threshold of **3.0 kcal/mol** is fixed in advance
(~4-8x the replicate SD seen for GPRP in the pilot, 0.4-0.8 kcal/mol).

**H1 -- the protocol is consistent with published determinants of hole-a binding.**
* H1a. GPR (tripeptide) is predicted to bind less tightly than GPRP
  (Laudano & Doolittle 1978): ddG(GPR) > 0, one-sided Welch p < 0.05 (n = 3 vs 6).
* H1b. Replacing hole residue Asp364 by Ala, His or Val weakens GPRP binding
  (Bowley et al. 2008; Ala partial, His/Val no detectable binding): for each
  mutant ddG >= +3.0 kcal/mol and one-sided Welch p < 0.05 vs the same-shell
  wild type (n = 3 vs 3; 7 A shell).
* H1c. Specificity control: Ala substitution of four **neutral** second-shell
  residues whose side chains lie 6-9 A from the peptide (Phe293, Ser300,
  Met336, Asn337; chosen by this rule, not by outcome) changes dG by less than
  3.0 kcal/mol (n = 2 each).
* "Consistent" = H1a and H1b (all three mutants) and H1c (all four controls)
  satisfied. Anything less is reported as a failure of the corresponding
  check, and scan results are then labelled *not validated*.
* Secondary, descriptive only: ddG(D364A) < ddG(D364H), ddG(D364V) (the assay
  reports partial vs no protection); order of GPRP, GPRV, GHRP (a patent
  summary ranks GPRP > GPRV > GHRP at hole a; weak provenance).

**H2 -- saturation scan identifies candidate tighter binders.** Scan all 19
substitutions at P4, P2, P3 and P1 (76 variants, n = 2).
* Stage-1 candidate: mean ddG <= -3.0 kcal/mol, *both* replicates <= -2.0, and
  mean backbone RMSD to the crystal pose < 3.0 A.
* Confirmed candidate (Phase C: seeds 3-6 at 6 A -> n = 6; 3 seeds at 8 A):
  n = 6 mean ddG <= -3.0 with Benjamini-Hochberg-adjusted Welch p < 0.05
  (across all tested variants), 8 A-shell ddG <= -2.0, and the sign of ddG
  preserved on rescoring with eps_in = 2 and 4 and with GBn2.
* H2 is **not supported** if no candidate is confirmed, *or* if H1 is not
  consistent (in which case confirmed candidates are reported as unvalidated
  model outputs).

**Diagnostics declared in advance (reported whatever they show):**
(i) rank correlation of ddG with side-chain heavy-atom count at each position
(size bias of MM/GBSA); (ii) per-position tolerance profile; (iii) rescoring
sensitivity (eps_in 1/2/4, GBn2) on the existing trajectories; (iv)
gas-phase vs solvation split of every ddG.

## 4. What these results can and cannot show
* Replicate SD measures sampling noise only. Systematic error (force field,
  GB model, truncation, no entropy, one crystal template, fixed protonation,
  restrained backbone) is **not** captured and is likely larger.
* The 3.0 kcal/mol threshold is a convention, not a calibrated error bar.
* Absolute dG values (tens of kcal/mol) are not comparable to experimental
  affinities; only differences between peptides evaluated with the same model
  are interpreted, and even those only as hypotheses for experiment.
* No claim is made about hemostatic efficacy.

## 5. Pilot summary (used to set the protocol)
| set | n | production | mean dG | between-replicate SD |
|---|---|---|---|---|
| GPRP (long) | 4 | 300 ps | -55.5 | 0.38 |
| GPRP (protocol-length, 30 ps equil.) | 2 | 100 ps | -55.4 | 0.79 |
| GPKP / GPWP / GADP (build check only) | 2 each | 100 ps | -45.4 / -30.7 / -22.8 | 0.3-2.6 |

Window means of the four long GPRP runs: 0-100 ps = **-53.9** (SD 1.5);
100-200 ps = -56.3 (SD 0.8); 200-300 ps = -56.4 (SD 0.5). The first 100 ps after
annealing still relaxed, so the equilibration was lengthened from 30 to
100 ps (deviation D1 below). Pilot wall time ~315 s per 144 ps job with 4
parallel jobs.

## 6. Deviation log
* **D1 (pre-production, pilot-driven):** equilibration 30 -> 100 ps because the
  first 100 ps window of the long GPRP runs was biased by ~2.4 kcal/mol.
  Applied to all production jobs alike.
