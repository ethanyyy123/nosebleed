# Reference verification log

The analysis sandbox could not open any journal, PubMed or PMC page, so
bibliographic details were checked against **web-search result summaries**
only (no paper full text was read). "Confirmed" below means the field appeared
in a search result for that paper in this session. Anything marked
**check** came from the author's recall or appeared in only one source and
**must be checked against the paper itself before submission**. Numbers taken
from literature (affinities, rankings) are listed separately in
`data/literature_values.csv` with the same flags.

| # | Reference (as cited in the manuscript) | Confirmed | Check |
|---|---|---|---|
| Tunkel | Tunkel DE, et al. Clinical practice guideline: nosebleed (epistaxis). *Otolaryngol Head Neck Surg* 2020;162(1 Suppl):S1-S38. | journal, year, volume, page range, first author | the "at least 60% of people in the US" figure appeared in AAO-HNSF press material, not the guideline text -- confirm in the guideline |
| Weisel | Weisel JW, Litvinov RI. Fibrin formation, structure and properties. *Subcell Biochem* 2017;82:405-456. doi:10.1007/978-3-319-49674-0_13 | all fields | -- |
| Laudano | Laudano AP, Doolittle RF. Synthetic peptide derivatives that bind to fibrinogen and prevent the polymerization of fibrin monomers. *Proc Natl Acad Sci USA* 1978;75(7):3085-3089. doi:10.1073/pnas.75.7.3085 | title, journal, year, volume/issue, first page, DOI, authors' surnames | last page (3089); the "association constant ~5x10^4 M^-1" and "GPR binds less tightly" statements came from a search summary |
| Spraggon | Spraggon G, Everse SJ, Doolittle RF. Crystal structures of fragment D from human fibrinogen and its crosslinked counterpart from fibrin. *Nature* 1997;389:455-462. doi:10.1038/38947 | all fields | -- |
| Pratt | Pratt KP, Cote HCF, Chung DW, Stenkamp RE, Davie EW. The primary fibrin polymerization pocket: three-dimensional structure of a 30-kDa C-terminal gamma chain fragment complexed with the peptide Gly-Pro-Arg-Pro. *Proc Natl Acad Sci USA* 1997;94(14):7176-7181. doi:10.1073/pnas.94.14.7176 | title, journal, year, volume/issue, pages, DOI, 2FIB at 2.01 A | author list (from recall) |
| Betts | Betts L, Merenbloom BK, Lord ST. The structure of fibrinogen fragment D with the 'A' knob peptide GPRVVE. *J Thromb Haemost* 2006;4(5):1139-1141. doi:10.1111/j.1538-7836.2006.01902.x | journal, volume, pages, DOI, PDB 2FFD (2.89 A), authors' surnames | author order/initials; a first draft wrongly said *Acta Cryst F* -- corrected |
| Stabenfeldt | Stabenfeldt SE, Gossett JJ, Barker TH. Building better fibrin knob mimics: an investigation of synthetic fibrin knob peptide structures in solution and their dynamic binding with fibrinogen/fibrin holes. *Blood* 2010;116(8):1352-1359. doi:10.1182/blood-2009-11-251801 | all fields (PMC2938242) | SPR numbers (GPRPFPAC ~1.1 uM, GPRPAAC ~8.5 uM, GPRP ~4x tighter than GPRVVE) are from a search summary; the summary says "second binding site" -- read the paper's tables before quoting any KD |
| Kononova | Kononova O, Litvinov RI, Zhmurov A, Alekseenko A, Cheng CH, Agarwal S, Marx KA, Weisel JW, Barsegov V. Molecular mechanisms, thermodynamics, and dissociation kinetics of knob-hole interactions in fibrin. *J Biol Chem* 2013;288:22681-22692. doi:10.1074/jbc.M113.472365 | all fields | the arXiv version swaps the hole labels; cite the JBC version |
| Asquith | Asquith NL, Duval C, Zhmurov A, Baker SR, McPherson HR, Domingues MM, Connell SDA, Barsegov V, Ariens RAS. Fibrin protofibril packing and clot stability are enhanced by extended knob-hole interactions and catch-slip bonds. *Blood Adv* 2022;6(13):4015-4027. doi:10.1182/bloodadvances.2022006977 | all fields | "GPRP binds hole a with KD ~25 uM (double-D) vs ~5.8 uM for larger fragments" came from a search summary -- confirm |
| Bowley | Bowley SR, Merenbloom BK, Okumura N, Betts L, Heroux A, Gorkun OV, Lord ST. Polymerization-defective fibrinogen variant gammaD364A binds knob 'A' peptide mimic. *Biochemistry* 2008;47(33):8607-8613. | all fields; PDB 3BVH 2.6 A | findings (D364A weakens, D364H/V abolish GPRP protection from plasmin) are from abstract-level summaries |
| Berman | Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. *Nucleic Acids Res* 2000;28(1):235-242. doi:10.1093/nar/28.1.235 | all fields | -- |
| Maier | Maier JA, Martinez C, Kasavajhala K, Wickstrom L, Hauser KE, Simmerling C. ff14SB: improving the accuracy of protein side chain and backbone parameters from ff99SB. *J Chem Theory Comput* 2015;11(8):3696-3713. doi:10.1021/acs.jctc.5b00255 | all fields | -- |
| Eastman | Eastman P, et al. OpenMM 8: molecular dynamics simulation with machine learning potentials. *J Phys Chem B* 2024;128(1):109-116. doi:10.1021/acs.jpcb.3c06662 | all fields | -- |
| Onufriev | Onufriev A, Bashford D, Case DA. Exploring protein native states and large-scale conformational changes with a modified generalized Born model. *Proteins* 2004;55:383-394. | title, journal, year, volume, pages | DOI |
| Genheden | Genheden S, Ryde U. The MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities. *Expert Opin Drug Discov* 2015;10(5):449-461. doi:10.1517/17460441.2015.1032936 | all fields | -- |
| Hou | Hou T, Wang J, Li Y, Wang W. Assessing the performance of the MM/PBSA and MM/GBSA methods. 1. The accuracy of binding free energy calculations based on molecular dynamics simulations. *J Chem Inf Model* 2011;51(1):69-82. doi:10.1021/ci100275a | all fields | -- |
| Trott | Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J Comput Chem* 2010;31(2):455-461. doi:10.1002/jcc.21334 | all fields | -- |

## Claims deliberately **not** made
* No claim that a designed peptide *is* a better hemostat, binder, or inhibitor.
* No quantitative comparison of computed numbers to experimental KD values
  (the primary tables were not readable here).
* Nasal-mucosal pH/ion values are not used, because no verified source for them
  was obtained.
