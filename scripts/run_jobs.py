"""run_jobs.py -- run MD + single-trajectory MM/GBSA jobs in parallel (resumable).

A job is a dict:
  {"label": "GPRP", "sequence": "GPRP", "rec_mut": {"364": "A"}, "shell": 6.0,
   "seed": 1, "equil_ps": 30, "prod_ps": 100, "anneal_ps": 2}
Each finished job writes data/results/raw/<label>__s<shell>__r<seed>.json
(+ a compressed trajectory in data/trajectories/, which is git-ignored).

usage: python scripts/run_jobs.py jobs.json [--workers 4]
"""
import json
import multiprocessing as mp
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "data", "results", "raw")
TRAJ = os.path.join(ROOT, "data", "trajectories")


def job_id(j):
    return f"{j['label']}__s{j['shell']:g}__r{j['seed']}"


def run_one(job):
    import numpy as np
    sys.path.insert(0, HERE)
    import knobhole_lib as kl
    from openmm import unit

    jid = job_id(job)
    out = os.path.join(RAW, jid + ".json")
    if os.path.exists(out):
        return jid, "skipped"
    t0 = time.time()
    try:
        rec_mut = {int(k): v for k, v in job.get("rec_mut", {}).items()}
        m = kl.build_model(job["sequence"], shell=job["shell"], rec_mutations=rec_mut)
        frames, meta = kl.run_md(m, seed=job["seed"], equil_ps=job.get("equil_ps", 30),
                                 prod_ps=job.get("prod_ps", 100), save_ps=5, threads=1,
                                 k_rest=job.get("k_rest", 1.0), anneal=job.get("anneal", True),
                                 anneal_ps=job.get("anneal_ps", 2))
        sc = kl.GBSAScorer(m).score(frames)

        # ---- structural descriptors relative to the crystal pose
        ref = np.array(m.positions.value_in_unit(unit.nanometer))
        atoms = list(m.topology.atoms())
        bb = [a.index for a in atoms[m.n_rec:] if a.name in ("N", "CA", "C") and a.residue.name != "NME"]
        rmsd = [10 * float(np.sqrt(((f[bb] - ref[bb]) ** 2).sum(1).mean())) for f in frames]
        d364 = [a.index for a in atoms[:m.n_rec] if a.residue.id == "364" and a.name in ("OD1", "OD2")]
        n1 = [a.index for a in atoms[m.n_rec:] if a.name == "N" and a.residue.index == atoms[m.n_rec].residue.index][0]
        dist_n_d364 = ([10 * float(np.linalg.norm(f[n1] - f[d364].mean(0))) for f in frames]
                       if d364 else None)
        pep = np.arange(m.n_rec, len(atoms))
        heavy_rec = np.array([i for i in range(m.n_rec) if atoms[i].element.symbol != "H"])
        heavy_pep = np.array([i for i in pep if atoms[i].element.symbol != "H"])
        contacts = []
        for f in frames:
            d = np.linalg.norm(f[heavy_pep][:, None] - f[heavy_rec][None], axis=2)
            contacts.append(int((d < 0.4).sum()))
        res = dict(job=job, jid=jid, n_atoms=len(atoms), meta=meta,
                   dg=sc["dg"].tolist(), dg_gas=sc["dg_gas"].tolist(), dg_gb=sc["dg_gb"].tolist(),
                   mean_dg=float(sc["dg"].mean()), sd_dg=float(sc["dg"].std(ddof=1)),
                   mean_gas=float(sc["dg_gas"].mean()), mean_gb=float(sc["dg_gb"].mean()),
                   half1=float(sc["dg"][:len(frames) // 2].mean()),
                   half2=float(sc["dg"][len(frames) // 2:].mean()),
                   rmsd_bb_A=rmsd, mean_rmsd_bb_A=float(np.mean(rmsd)),
                   dist_Nterm_D364_A=dist_n_d364, contacts_lt4A=contacts,
                   wall_total_s=round(time.time() - t0, 1))
        os.makedirs(TRAJ, exist_ok=True)
        np.savez_compressed(os.path.join(TRAJ, jid + ".npz"), frames=frames,
                            n_rec=m.n_rec)
        if job.get("save_pdb"):
            os.makedirs(os.path.join(ROOT, "data", "results", "pdb"), exist_ok=True)
            kl.write_pdb(m, frames[-1], os.path.join(ROOT, "data", "results", "pdb", jid + ".pdb"))
        with open(out, "w") as fh:
            json.dump(res, fh)
        return jid, f"ok  dG={res['mean_dg']:.1f}  rmsd={res['mean_rmsd_bb_A']:.2f}  {res['wall_total_s']}s"
    except Exception:
        return jid, "ERROR " + traceback.format_exc().replace("\n", " | ")[-400:]


if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True)
    jobs = json.load(open(sys.argv[1]))
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 4
    todo = [j for j in jobs if not os.path.exists(os.path.join(RAW, job_id(j) + ".json"))]
    print(f"{len(jobs)} jobs, {len(todo)} to run on {workers} workers", flush=True)
    ctx = mp.get_context("spawn")
    t0 = time.time()
    with ctx.Pool(workers) as pool:
        for k, (jid, status) in enumerate(pool.imap_unordered(run_one, todo), 1):
            print(f"[{k}/{len(todo)} {time.time() - t0:.0f}s] {jid}: {status}", flush=True)
