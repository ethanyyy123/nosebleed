"""make_jobs.py -- write the phased job lists for the study (frozen protocol).

Protocol (all jobs): 4-fs Langevin MD, 300 K, anneal 7 x 2 ps (up to 450 K),
100 ps equilibration (discarded), 100 ps production sampled every 5 ps
(20 frames), C-alpha restraint 1 kcal/mol/A^2, peptide leash 6 A.
Phase A: reference + validation panel.   Phase B: saturation scan (2 seeds).
Phase C (written later by 06_analyze.py): confirmation of candidates.
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROTO = dict(equil_ps=100, prod_ps=100, anneal_ps=2)
AA = "ACDEFGHIKLMNPQRSTVWY"
WT = "GPRP"


def job(label, seq, seed, shell=6.0, rec_mut=None, **kw):
    j = dict(label=label, sequence=seq, shell=shell, seed=seed, **PROTO, **kw)
    if rec_mut:
        j["rec_mut"] = {str(k): v for k, v in rec_mut.items()}
    return j


def phase_a():
    jobs = []
    for s in range(1, 7):
        jobs.append(job("GPRP", WT, s))                       # reference, n = 6
    for s in range(1, 4):
        jobs.append(job("GPR", "GPR", s))                     # H1a
    for s in range(1, 4):                                     # H1b/c at 7 A shell
        jobs.append(job("GPRP", WT, s, shell=7.0))
        for name, mut in (("D364A", {364: "A"}), ("D364H", {364: "H"}), ("D364V", {364: "V"})):
            jobs.append(job(f"GPRP_{name}", WT, s, shell=7.0, rec_mut=mut))
    for s in (1, 2):
        for name, mut in (("F293A", {293: "A"}), ("S300A", {300: "A"}),
                          ("M336A", {336: "A"}), ("N337A", {337: "A"})):
            jobs.append(job(f"GPRP_{name}", WT, s, shell=7.0, rec_mut=mut))
    return jobs


def phase_b(positions=(4, 2, 3, 1), seeds=(1, 2)):
    jobs = []
    for pos in positions:
        for aa in AA:
            if aa == WT[pos - 1]:
                continue
            seq = WT[:pos - 1] + aa + WT[pos:]
            for s in seeds:
                jobs.append(job(seq, seq, s))
    return jobs


if __name__ == "__main__":
    os.makedirs(os.path.join(ROOT, "data", "jobs"), exist_ok=True)
    for name, jobs in (("phaseA", phase_a()), ("phaseB", phase_b())):
        path = os.path.join(ROOT, "data", "jobs", f"{name}.json")
        json.dump(jobs, open(path, "w"), indent=0)
        print(name, len(jobs), "jobs ->", path)
