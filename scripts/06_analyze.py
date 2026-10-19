"""06_analyze.py -- pre-specified analyses (see docs/02_study_design.md).

usage:
  python scripts/06_analyze.py h1        # consistency panel (Phase A)
  python scripts/06_analyze.py scan      # saturation scan statistics + stage-1 candidates
  python scripts/06_analyze.py confirm   # Phase C confirmation of candidates
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import aggregate as ag  # noqa: E402

ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "data", "results")
THR = 3.0            # practical-significance threshold (kcal/mol), fixed in advance
WT = "GPRP"
SC_HEAVY = dict(G=0, A=1, S=2, C=2, V=3, T=3, P=3, D=4, N=4, I=4, L=4, M=4, E=5, Q=5,
                K=5, H=6, F=7, R=7, Y=8, W=10)


def one_sided(a, b, alt):
    t, p = stats.ttest_ind(a, b, equal_var=False, alternative=alt)
    return float(p)


def h1_panel(df):
    rows = []
    ref6 = df[(df.label == "GPRP") & (df.shell == 6)].dg.values
    ref7 = df[(df.label == "GPRP") & (df.shell == 7)].dg.values
    # H1a
    g = df[(df.label == "GPR") & (df.shell == 6)].dg.values
    if len(g) >= 2 and len(ref6) >= 2:
        d = g.mean() - ref6.mean()
        p = one_sided(g, ref6, "greater")
        rows.append(dict(test="H1a", comparison="GPR vs GPRP (6 A)", n_test=len(g), n_ref=len(ref6),
                         mean_test=g.mean(), mean_ref=ref6.mean(), ddG=d, p_one_sided=p,
                         criterion="ddG>0 and p<0.05", passed=bool(d > 0 and p < 0.05)))
    # H1b
    for m in ("D364A", "D364H", "D364V"):
        x = df[(df.label == f"GPRP_{m}") & (df.shell == 7)].dg.values
        if len(x) >= 2 and len(ref7) >= 2:
            d = x.mean() - ref7.mean()
            p = one_sided(x, ref7, "greater")
            rows.append(dict(test="H1b", comparison=f"{m} vs WT (7 A)", n_test=len(x), n_ref=len(ref7),
                             mean_test=x.mean(), mean_ref=ref7.mean(), ddG=d, p_one_sided=p,
                             criterion=f"ddG>={THR} and p<0.05", passed=bool(d >= THR and p < 0.05)))
    # H1c
    for m in ("F293A", "S300A", "M336A", "N337A"):
        x = df[(df.label == f"GPRP_{m}") & (df.shell == 7)].dg.values
        if len(x) >= 1 and len(ref7) >= 2:
            d = x.mean() - ref7.mean()
            rows.append(dict(test="H1c", comparison=f"{m} vs WT (7 A)", n_test=len(x), n_ref=len(ref7),
                             mean_test=x.mean(), mean_ref=ref7.mean(), ddG=d, p_one_sided=np.nan,
                             criterion=f"|ddG|<{THR}", passed=bool(abs(d) < THR)))
    t = pd.DataFrame(rows)
    # secondary / descriptive
    sec = []
    wt7 = ref7.mean() if len(ref7) else np.nan
    a = {m: df[(df.label == f"GPRP_{m}") & (df.shell == 7)].dg.mean() - wt7 for m in ("D364A", "D364H", "D364V")}
    sec.append(dict(item="ddG order D364A vs D364H/V (assay: partial vs none)", values=a))
    for lab in ("GPRV", "GHRP"):
        x = df[(df.label == lab) & (df.shell == 6)].dg.values
        if len(x):
            sec.append(dict(item=f"ddG {lab} vs GPRP", values={"ddG": float(x.mean() - ref6.mean()), "n": len(x)}))
    return t, sec


def scan_table(df, ref_label="GPRP"):
    d6 = df[df.shell == 6]
    ref = d6[d6.label == ref_label].dg.values
    base = ref.mean()
    rows = []
    for seq, g in d6[(d6.label == d6.sequence) & (d6.sequence != ref_label)].groupby("sequence"):
        if len(seq) != 4:
            continue
        v = g.dg.values
        pos = [i for i in range(4) if seq[i] != WT[i]]
        if len(pos) != 1:
            continue
        i = pos[0]
        _, p = stats.ttest_ind(v, ref, equal_var=False) if len(v) >= 2 else (np.nan, np.nan)
        rows.append(dict(sequence=seq, position=i + 1, aa=seq[i], n=len(v), dg_mean=v.mean(),
                         dg_sd=v.std(ddof=1) if len(v) > 1 else np.nan, ddG=v.mean() - base,
                         ddG_min=v.min() - base, ddG_max=v.max() - base,
                         gas=g.dg_gas.mean() - d6[d6.label == ref_label].dg_gas.mean(),
                         solv=g.dg_gb.mean() - d6[d6.label == ref_label].dg_gb.mean(),
                         rmsd=g.rmsd.mean(), max_rmsd=g.max_rmsd.max(), p_welch=p,
                         sc_heavy=SC_HEAVY[seq[i]], sc_heavy_wt=SC_HEAVY[WT[i]]))
    t = pd.DataFrame(rows)
    if len(t):
        t["p_bh"] = ag.bh(t.p_welch.fillna(1.0).values)
        t["stage1_candidate"] = (t.ddG <= -THR) & (t.ddG_max <= -2.0) & (t.rmsd < 3.0)
    return t, base, len(ref)


def diagnostics(t):
    out = {}
    for pos, g in t.groupby("position"):
        if len(g) >= 5:
            rho, p = stats.spearmanr(g.sc_heavy, g.ddG)
            out[f"P{pos}"] = dict(n=len(g), rho_size_vs_ddG=float(rho), p=float(p))
    if len(t) >= 8:
        rho, p = stats.spearmanr(t.sc_heavy - t.sc_heavy_wt, t.ddG)
        out["pooled"] = dict(n=len(t), rho_dsize_vs_ddG=float(rho), p=float(p))
        rho, p = stats.spearmanr(t.gas, t.ddG)
        out["gas_vs_total"] = dict(rho=float(rho), p=float(p))
        rho, p = stats.spearmanr(t.solv, t.ddG)
        out["solv_vs_total"] = dict(rho=float(rho), p=float(p))
    return out


def main(mode):
    df = ag.load_replicates("*")
    df = df[~df.label.str.startswith("PILOT")]
    df.round(4).to_csv(os.path.join(RES, "replicates_all.csv"), index=False)
    if mode == "h1":
        t, sec = h1_panel(df)
        t.round(4).to_csv(os.path.join(RES, "h1_panel.csv"), index=False)
        print(t.round(3).to_string(index=False))
        print(json.dumps(sec, indent=1, default=float))
        if len(t):
            ok = {k: bool(t[t.test == k].passed.all()) for k in ("H1a", "H1b", "H1c") if (t.test == k).any()}
            print("component outcomes:", ok, "| all consistent:", all(ok.values()) and len(ok) == 3)
    elif mode == "scan":
        t, base, nref = scan_table(df)
        t.round(4).sort_values("ddG").to_csv(os.path.join(RES, "scan_ddG.csv"), index=False)
        print(f"reference GPRP mean dG = {base:.2f} (n={nref}); variants scored: {len(t)}")
        print(t.sort_values("ddG").head(15).round(2).to_string(index=False))
        print(t.sort_values("ddG").tail(5).round(2).to_string(index=False))
        cand = t[t.stage1_candidate].sequence.tolist()
        print("stage-1 candidates:", cand)
        json.dump(cand, open(os.path.join(RES, "stage1_candidates.json"), "w"))
        diag = diagnostics(t)
        json.dump(diag, open(os.path.join(RES, "scan_diagnostics.json"), "w"), indent=1)
        print(json.dumps(diag, indent=1))
    elif mode == "confirm":
        cand = json.load(open(os.path.join(RES, "stage1_candidates.json")))
        # BH is applied across ALL scanned variants (candidates now carry n = 6)
        t_all, base, nref = scan_table(df)
        ref8 = df[(df.label == "GPRP") & (df.shell == 8)].dg.values
        rows = []
        for seq in cand:
            r = t_all[t_all.sequence == seq].iloc[0]
            x8 = df[(df.label == seq) & (df.shell == 8)].dg.values
            ddg8 = (x8.mean() - ref8.mean()) if len(x8) and len(ref8) else np.nan
            rows.append(dict(sequence=seq, position=int(r.position), aa=r.aa, n6=int(r.n),
                             ddG6=r.ddG, sd6=r.dg_sd, p_welch=r.p_welch, p_bh_all=r.p_bh,
                             n8=len(x8), ddG8=ddg8, rmsd=r.rmsd,
                             confirmed_6A=bool(r.ddG <= -THR and r.p_bh < 0.05 and r.n >= 6),
                             confirmed_8A=bool(ddg8 <= -2.0) if not np.isnan(ddg8) else False))
        t = pd.DataFrame(rows)
        t.round(4).to_csv(os.path.join(RES, "confirmation.csv"), index=False)
        print(t.round(3).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1])
