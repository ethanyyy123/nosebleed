---
STATUS: DRAFT IN PROGRESS. Results and Discussion are [PENDING] until Phase A/B/C
MD runs finish. Nothing in brackets is a real number. Formatted for Reinvention:
an International Journal of Undergraduate Research (see docs/01 for the format
rules and their provenance; [VERIFY] against the live style guide before
submission). Word count excludes title page, abstract, acknowledgements,
references, appendices per Reinvention's rule.
---

# An ensemble MM/GBSA screen of single-residue fibrin "knob A" peptide variants at fibrinogen hole a: a validated computational pipeline and a candidate-generating saturation scan

**Author:** [NAME] -- [Department], [University] -- [email]

## Abstract (100-200 words)

[PENDING -- written last, after Results/Discussion are final.]

**Keywords:** fibrin, fibrinogen, molecular dynamics, MM/GBSA, knob-hole interaction, peptide design, hemostasis

---

## Introduction

Fibrin clot formation begins when thrombin cleaves fibrinopeptides A and B from
circulating fibrinogen, exposing new N-termini -- "knobs A and B" -- on the
central (E) region of the resulting fibrin monomer (Weisel and Litvinov, 2017).
These knobs dock into complementary "holes a and b" on the outer (D) regions of
neighbouring fibrin(ogen) molecules, nucleating the half-staggered protofibrils
that assemble into the fibrin mesh of a clot (Weisel and Litvinov, 2017). The
short synthetic peptide Gly-Pro-Arg-Pro (GPRP), modelled on the new N-terminus
exposed at knob A, reproduces this interaction well enough to competitively
block fibrin polymerisation in solution and was the ligand used to solve the
first crystal structures of the hole-a pocket (Laudano and Doolittle, 1978;
Pratt et al., 1997; Spraggon et al., 1997).

Because knob-hole docking is the first committed step of clot assembly, it is
an attractive target for peptide-based topical hemostats -- materials designed
to accelerate mesh formation at a bleeding site, such as a nosebleed, without
relying on the patient's own coagulation cascade (see e.g. patent
US8513380B2, cited for context, not as a source of quantitative data). A
"tighter-binding knob A" is one plausible lever for such a material: a peptide
that engages hole a more avidly than the native sequence could, in principle,
compete more effectively for the pocket or template fibrin-like assembly more
efficiently. Whether that intuition holds for any particular substitution is
an empirical question, and it is one that structure-based computation can
address cheaply, before any peptide is synthesised.

Two bodies of prior work bound this project. First, the hole-a pocket and its
determinants are unusually well characterised structurally: multiple crystal
structures of fibrinogen fragment D in complex with GPRP and related peptides
are available (Pratt et al., 1997; Spraggon et al., 1997; Everse et al., 1998),
and site-directed mutagenesis has identified γAsp364 as electrostatically
important for knob-A binding, with alanine, histidine and valine substitutions
at that position progressively weakening peptide protection from plasmin
digestion in solution even though the alanine variant's crystal structure
retains near-native knob-A contacts (Bowley et al., 2008). Second, molecular
simulation of this exact interaction has already been done at least twice --
coarse structural/kinetic modelling of designed knob mimics by surface
plasmon resonance plus simulation (Stabenfeldt et al., 2010), and all-atom
steered-MD unbinding simulations of the GPRP:a and GHRP:b complexes under
varying pH and temperature (Kononova et al., 2013) -- but neither study
reports a systematic, replicate-resolved free-energy ranking of single-residue
substitutions of the knob-A peptide itself. That is the gap this project
targets.

This project asks two questions. First (H1), can an inexpensive, CPU-only,
truncated-receptor MM/GBSA protocol reproduce three independent, already-
published determinants of hole-a binding -- that the GPR tripeptide binds more
weakly than the GPRP tetrapeptide (Laudano and Doolittle, 1978), that
substituting the pocket's γAsp364 weakens binding (Bowley et al., 2008), and
that substituting *distant*, non-contacting pocket residues does not -- well
enough that its outputs on untested variants should be believed at all?
Second (H2), if so, does a saturation scan of single substitutions at each of
GPRP's four positions identify candidate variants predicted to bind more
tightly than the native sequence, and do those candidates survive independent
checks (more replicates, a larger receptor shell, alternative solvation
models)? Full, numeric, pre-registered decision rules for both questions --
fixed in writing before any production simulation was examined -- are given in
Methods and in the repository's `docs/02_study_design.md`.

This is explicitly a computational screening study. No peptide reported here
has been synthesised or tested experimentally; every number is a model output
with its own, only partly quantified, systematic error (Methods,
Limitations). The paper's contribution is not a claim that any specific
variant *is* a better hemostatic knob, but (a) a reusable, freely available,
single-CPU pipeline for this class of problem, (b) an explicit, falsifiable
check of whether that pipeline's outputs track known biology before it is
trusted on anything new, and (c) a ranked, statistically qualified list of
candidates for anyone who wants to test them.

## Methods

**Overview.** All computation used a 2024-vintage snapshot of the RCSB Protein
Data Bank (PDB), AMBER ff14SB in implicit solvent, and OpenMM 8.6 running on
4 CPU cores with no GPU. Every setting below was fixed in writing (frozen
protocol, `docs/02_study_design.md`) before any comparison of interest
(validation panel, saturation scan) was examined; the one deviation made
*after* seeing pilot data, and why, is logged there. All code, raw per-replicate
data, and figures are at [repository URL -- GitHub, see Data Availability].

**Structural template.** The receptor and reference ligand came from PDB entry
2FIB, a 2.01 A crystal structure of the human fibrinogen γ-module (residues
143-411) bound to GPRP (Pratt et al., 1997). A pocket model was built by
keeping every receptor residue with a heavy atom within 6 A of the template
GPRP, plus the structure's bound Ca²⁺ ion and its direct ligands whenever any
of them fell within that shell; the resulting fragments were capped with
acetyl (ACE) and N-methylamide (NME) groups built from the coordinates of the
real, now-omitted neighbouring residues, so that no cut carried a spurious
free charge. The same procedure, with the shell radius as the only free
parameter, was re-used at 7 A (for the validation panel, Methods below) and
8 A (for a sensitivity check on the scan's candidates).

**Peptide variants.** Each 3-4-residue peptide was built by substituting
residues onto the template GPRP backbone position-by-position (so that
backbone geometry at every position is inherited from the solved structure,
not re-built from an idealised geometry), free at the N-terminus (as the new
N-terminus exposed by thrombin cleavage) and capped with N-methylamide at the
C-terminus so the peptide mimics a continuing chain rather than a free acid.
Missing side-chain atoms and all hydrogens were added with PDBFixer at pH 7.4.

**Force field and solvent.** Protein and peptide used AMBER ff14SB (Maier et
al., 2015); the Ca²⁺ ion used the AMBER/OpenMM Li-Merz parameter set
distributed with OpenMM's `amber14` force-field files. Solvent was treated
implicitly with the OBC-II generalized-Born model (Onufriev et al., 2004) as
implemented natively in OpenMM (`GBSAOBCForce`; interior/exterior dielectric
1/78.5; no explicit salt term). This choice trades the greater physical realism
of explicit solvent for roughly two orders of magnitude more sampling per CPU-
hour, which is what makes a several-hundred-variant scan tractable on 4 cores
(Results, Table S1); the trade-off, and what it costs in accuracy, is
discussed explicitly under Limitations.

**Restraints and dynamics.** To keep the truncated receptor stable without a
full-protein simulation, every receptor Cα atom (and the Ca²⁺ ion) was held by
a harmonic restraint (1 kcal mol⁻¹ A⁻²); side chains and the entire peptide were
otherwise free. A flat-bottom restraint on the peptide backbone centroid (6 A
radius, inert inside that radius) prevented the unbound peptide from drifting
away during the longer anneal step without biasing its pose within the pocket.
Each system was energy-minimised, annealed through 150-450-300 K (7 x 2 ps, to
independently randomise side-chain and peptide starting conformations across
replicates), then simulated with 4-fs Langevin dynamics (300 K, hydrogen mass
repartitioned to 3 amu, all bonds to hydrogen constrained) for 100 ps of
discarded equilibration followed by 100 ps of production, saved every 5 ps (20
frames). The equilibration length was set to 100 ps based on a pilot
convergence check (four independent 300-ps runs of the reference peptide
showed the first 100 ps still relaxing by ~2.4 kcal/mol relative to the 100-300
ps window; `docs/02_study_design.md`, deviation D1) -- i.e. it was lengthened
*before* any production comparison was examined, in response to a diagnostic
check, not a result.

**Scoring.** Binding was scored by single-trajectory MM/GBSA (Genheden and
Ryde, 2015; Hou et al., 2011): for every production frame, dG = G(complex) -
G(receptor) - G(peptide), each term evaluated with the same OBC-II model on
the same frame (no separate receptor/peptide trajectories, no entropy term).
The unit of statistical analysis is one independent MD replicate's mean dG
over its 20 production frames -- not an individual frame, because frames
within one replicate are autocorrelated -- so every reported n is a number of
independent restarts (distinct random seeds, independent velocity
reassignment and annealing), not a frame count.

**Hypotheses and pre-registered decision rules (abbreviated; full text and
exact thresholds in `docs/02_study_design.md`).**
*H1 (does the pipeline reproduce known biology?)*: (a) the tripeptide GPR
scores worse than GPRP (n=3 vs 6, one-sided Welch's t, alpha=0.05); (b) in
silico Ala/His/Val substitution of γAsp364 (7 A shell) each weakens binding by
>=3.0 kcal/mol relative to the same-shell wild type (n=3 vs 3, one-sided
Welch's t, alpha=0.05) -- a magnitude chosen as roughly 4-8x the reference
peptide's between-replicate noise observed in pilot runs (0.4-0.8 kcal/mol),
not a calibrated error bar; (c) Ala substitution of four neutral, non-
contacting second-shell residues (Phe293, Ser300, Met336, Asn337; selected by
the rule "side chain >=6 A from the peptide," before any result was seen) each
changes dG by <3.0 kcal/mol (specificity control). All three components must
hold for the pipeline to be called consistent with known biology; any
saturation-scan finding is reported as *unvalidated* otherwise.
*H2 (does a saturation scan find real candidates?)*: single substitutions at
each of GPRP's four positions (76 variants, n=2 each) are screened; a
*stage-1 candidate* needs mean ddG <= -3.0 kcal/mol with both replicates
<= -2.0 and backbone RMSD < 3.0 A from the crystal pose; a *confirmed
candidate* (Phase C) needs that result to hold at n=6 with a
Benjamini-Hochberg-adjusted p < 0.05 across all tested variants, a
qualitatively consistent (<= -2.0 kcal/mol) result in the independent 8 A
pocket model, and an unchanged sign on rescoring the same trajectories with
the solute dielectric raised to 2 and 4 and with the independent GBn2
generalized-Born model (Methods above).

**Statistics.** Between-group comparisons used Welch's unequal-variance t-test
on replicate means; multiple comparisons across the 76-variant scan were
corrected with the Benjamini-Hochberg procedure. Two diagnostics were
pre-specified to check for a trivial confound before interpreting any scan
result: (i) the Spearman correlation of a substitution's predicted ddG with
its side-chain heavy-atom count (an MM/GBSA size/polarity bias is a known
failure mode of this method, Hou et al., 2011; Genheden and Ryde, 2015), at
each peptide position and pooled; (ii) the gas-phase and solvation components
of ddG, reported separately, since the two can trade off in ways that a total-
energy number alone can hide.

**Software.** OpenMM 8.6 (Eastman et al., 2024) with PDBFixer; structures
fetched from the wwPDB AWS Open Data snapshot (2024-01-01) rather than
rcsb.org, which was not reachable from the analysis environment; AutoDock
Vina and Open Babel were installed for an earlier-considered docking-based
arm of this project but are not used in the final pipeline (see Discussion of
alternative approaches considered, `docs/01`). All analysis code, the frozen
study design, and every raw per-replicate result are openly available (Data
Availability).

---

## Results

[PENDING -- awaiting Phase A (validation panel, H1) and Phase B/C (saturation
scan and confirmation, H2) completion. Planned structure:
1. Pipeline throughput and pilot convergence (one short paragraph + Table S1).
2. H1 panel: table of GPR vs GPRP, the three Asp364 mutants vs WT, the four
   neutral controls vs WT, with means, n, ddG, p-values, and the explicit
   pass/fail call for each of the three H1 components -- Figure 2.
3. If and only if H1 is consistent: H2 saturation-scan heatmap (4 positions x
   19 substitutions, Figure 3), stage-1 candidates, size-bias and gas/solvation
   diagnostics, and Phase C confirmation table (Figure 4 / Table 2).
4. If H1 is *not* fully consistent: H2 results are still reported, explicitly
   labelled unvalidated, with the failing H1 component stated plainly.]

## Discussion

[PENDING -- to cover: (i) interpretation of the H1 outcome and what it does
and does not license claiming about H2; (ii) for any confirmed H2 candidate,
what a reasonable next experimental step would look like (e.g. SPR or ITC
against recombinant fragment D, as in Stabenfeldt et al., 2010); (iii)
Limitations: implicit solvent and no entropy term (Genheden and Ryde, 2015),
a single static crystal template and restrained backbone, one force field and
one GB flavour as primary (partially addressed by the GBn2/dielectric
sensitivity check on confirmed candidates only), no claim whatsoever about
hemostatic efficacy, wound biology, or in vivo behaviour of any peptide; (iv)
why the two other candidate project ideas considered (coarse-grained RADA16
self-assembly under pH/ionic variation; docking of plant polyphenols against
plasmin) were not pursued -- summarised from `docs/01_journal_venue_requirements.md`.]

## Acknowledgements

[PENDING]

## Data Availability

All code, the pre-registered study design, raw per-replicate results (JSON),
trajectories, and figures are openly available at
[GitHub repository URL -- ethanyyy123/nosebleed, branch/release to be pinned
at submission].

## References

*(Reinvention's own Harvard variant; [VERIFY] exact style against the live
guide before submission. Full verification status and exact wording checked
for each entry is in `docs/03_reference_verification.md`.)*

Bowley, S.R., Merenbloom, B.K., Okumura, N., Betts, L., Heroux, A., Gorkun,
O.V. and Lord, S.T. (2008) 'Polymerization-defective fibrinogen variant
γD364A binds knob "A" peptide mimic', *Biochemistry*, 47(33), pp. 8607-8613.

Eastman, P. et al. (2024) 'OpenMM 8: molecular dynamics simulation with
machine learning potentials', *Journal of Physical Chemistry B*, 128(1), pp.
109-116. doi: 10.1021/acs.jpcb.3c06662.

Everse, S.J., Spraggon, G., Veerapandian, L., Riley, M. and Doolittle, R.F.
(1998) 'Crystal structure of fragment double-D from human fibrin with two
different bound ligands', *Biochemistry*, 37(24), pp. 8637-8642. doi:
10.1021/bi9804129.

Genheden, S. and Ryde, U. (2015) 'The MM/PBSA and MM/GBSA methods to estimate
ligand-binding affinities', *Expert Opinion on Drug Discovery*, 10(5), pp.
449-461. doi: 10.1517/17460441.2015.1032936.

Hou, T., Wang, J., Li, Y. and Wang, W. (2011) 'Assessing the performance of
the MM/PBSA and MM/GBSA methods. 1. The accuracy of binding free energy
calculations based on molecular dynamics simulations', *Journal of Chemical
Information and Modeling*, 51(1), pp. 69-82. doi: 10.1021/ci100275a.

Kononova, O., Litvinov, R.I., Zhmurov, A., Alekseenko, A., Cheng, C.H.,
Agarwal, S., Marx, K.A., Weisel, J.W. and Barsegov, V. (2013) 'Molecular
mechanisms, thermodynamics, and dissociation kinetics of knob-hole
interactions in fibrin', *Journal of Biological Chemistry*, 288, pp.
22681-22692. doi: 10.1074/jbc.M113.472365.

Laudano, A.P. and Doolittle, R.F. (1978) 'Synthetic peptide derivatives that
bind to fibrinogen and prevent the polymerization of fibrin monomers',
*Proceedings of the National Academy of Sciences USA*, 75(7), pp. 3085-3089.
doi: 10.1073/pnas.75.7.3085.

Maier, J.A., Martinez, C., Kasavajhala, K., Wickstrom, L., Hauser, K.E. and
Simmerling, C. (2015) 'ff14SB: improving the accuracy of protein side chain
and backbone parameters from ff99SB', *Journal of Chemical Theory and
Computation*, 11(8), pp. 3696-3713. doi: 10.1021/acs.jctc.5b00255.

Onufriev, A., Bashford, D. and Case, D.A. (2004) 'Exploring protein native
states and large-scale conformational changes with a modified generalized
Born model', *Proteins*, 55, pp. 383-394.

Pratt, K.P., Cote, H.C.F., Chung, D.W., Stenkamp, R.E. and Davie, E.W. (1997)
'The primary fibrin polymerization pocket: three-dimensional structure of a
30-kDa C-terminal gamma chain fragment complexed with the peptide
Gly-Pro-Arg-Pro', *Proceedings of the National Academy of Sciences USA*,
94(14), pp. 7176-7181. doi: 10.1073/pnas.94.14.7176.

Spraggon, G., Everse, S.J. and Doolittle, R.F. (1997) 'Crystal structures of
fragment D from human fibrinogen and its crosslinked counterpart from
fibrin', *Nature*, 389, pp. 455-462. doi: 10.1038/38947.

Stabenfeldt, S.E., Gossett, J.J. and Barker, T.H. (2010) 'Building better
fibrin knob mimics: an investigation of synthetic fibrin knob peptide
structures in solution and their dynamic binding with fibrinogen/fibrin
holes', *Blood*, 116(8), pp. 1352-1359. doi: 10.1182/blood-2009-11-251801.

Weisel, J.W. and Litvinov, R.I. (2017) 'Fibrin formation, structure and
properties', in Parry, D.A.D. and Squire, J.M. (eds.) *Fibrous Proteins:
Structures and Mechanisms*, Subcellular Biochemistry, 82. Cham: Springer, pp.
405-456. doi: 10.1007/978-3-319-49674-0_13.
