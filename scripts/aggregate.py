"""aggregate.py -- collect raw per-job JSON into replicate-level tables.

Unit of analysis = one independent MD replicate (mean MM/GBSA over its
production frames).  Per-variant statistics are computed across replicates, never
across frames, because frames within a replicate are autocorrelated.
"""
from __future__ import annotations

import glob
import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "results", "raw")


def load_replicates(pattern="*", raw_dir=RAW):
    rows = []
    for f in sorted(glob.glob(os.path.join(raw_dir, pattern + ".json"))):
        r = json.load(open(f))
        j = r["job"]
        rows.append(dict(
            jid=r["jid"], label=j["label"], sequence=j["sequence"],
            rec_mut=",".join(f"{k}{v}" for k, v in sorted(j.get("rec_mut", {}).items())),
            shell=j["shell"], seed=j["seed"], prod_ps=r["meta"]["prod_ps"],
            dg=r["mean_dg"], dg_gas=r["mean_gas"], dg_gb=r["mean_gb"],
            sd_frames=r["sd_dg"], half1=r["half1"], half2=r["half2"],
            rmsd=r["mean_rmsd_bb_A"], max_rmsd=max(r["rmsd_bb_A"]),
            d_n_asp364=(np.mean(r["dist_Nterm_D364_A"]) if r["dist_Nterm_D364_A"] else np.nan),
            contacts=float(np.mean(r["contacts_lt4A"])), wall_s=r["wall_total_s"],
        ))
    return pd.DataFrame(rows)


def variant_table(df, by=("label", "sequence", "rec_mut", "shell")):
    g = df.groupby(list(by))
    t = g.agg(n=("dg", "size"), dg_mean=("dg", "mean"), dg_sd=("dg", lambda x: x.std(ddof=1)),
              gas=("dg_gas", "mean"), gb=("dg_gb", "mean"), rmsd=("rmsd", "mean"),
              max_rmsd=("max_rmsd", "max"), contacts=("contacts", "mean")).reset_index()
    t["dg_sem"] = t["dg_sd"] / np.sqrt(t["n"])
    return t


def welch(a, b):
    """Welch t statistic, df and two-sided p for two small samples."""
    from scipy import stats
    a, b = np.asarray(a, float), np.asarray(b, float)
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return float(t), float(p)


def bh(pvals):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n)
    prev = 1.0
    for rank, idx in zip(range(n, 0, -1), order[::-1]):
        prev = min(prev, p[idx] * n / rank)
        adj[idx] = prev
    return adj
