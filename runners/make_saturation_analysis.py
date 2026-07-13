"""Analysis + figures for the saturation matrices (alpha/beta/tails/corr/foldsdeep/
budgetfull/bigd). Robust to partial data.

    /tmp/gsv_venv/bin/python runners/make_saturation_analysis.py

Writes results/analysis/sat_*.csv, SATURATION.md, and results/figures/sat_*.pdf.
"""
import os, glob, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXP = os.path.join(ROOT, "results", "experiments")
OUTA = os.path.join(ROOT, "results", "analysis")
OUTF = os.path.join(ROOT, "results", "figures")
os.makedirs(OUTA, exist_ok=True); os.makedirs(OUTF, exist_ok=True)
FORM = {"paper_so": "SO", "paper_ro_ellipsoid": "RO", "paper_saa": "SAA",
        "paper_dro_wasserstein": "Wasserstein", "paper_fast": "FAST"}
PROP = ["UG", "NGS", "UNGS"]


def load(matrix):
    rows = []
    for f in sorted(glob.glob(os.path.join(EXP, matrix, "*_summary.json"))):
        j = json.load(open(f))
        dp = j.get("dgp_params") or {}
        for m, s in j["summaries"].items():
            if isinstance(s, dict) and "coverage" in s:
                rows.append({"config": j["config"], "n": j["n"], "split": j["split"], "d": j["d"],
                             "alpha": j.get("alpha"), "beta": j.get("beta"),
                             "df": dp.get("df"), "corr": dp.get("corr"),
                             "method": m, "coverage": s["coverage"], "mean_obj": s["mean_obj"],
                             "target": s.get("target"), "meets": s.get("meets_target")})
    return pd.DataFrame(rows)


def _line(df, xcol, title, fname, methods=PROP, target_from="col", xlog=False):
    cfgs = [c for c in FORM if c in df.config.unique()]
    if not cfgs:
        return
    fig, axes = plt.subplots(1, len(cfgs), figsize=(4.3 * len(cfgs), 4), squeeze=False)
    for ax, cfg in zip(axes[0], cfgs):
        g = df[df.config == cfg]
        for m in methods:
            gm = g[g.method == m].sort_values(xcol)
            if gm.empty:
                continue
            ax.plot(gm[xcol], gm.coverage, marker="o", ms=4, label=m)
        if target_from == "col" and "target" in g and g.target.notna().any():
            gt = g[g.method == methods[0]].sort_values(xcol)
            ax.plot(gt[xcol], gt.target, "k--", lw=1, label="target 1-beta")
        else:
            ax.axhline(0.95, ls="--", color="r", lw=1, label="target 0.95")
        if xlog:
            ax.set_xscale("log")
        ax.set_xlabel(xcol); ax.set_ylabel("coverage"); ax.set_title(FORM.get(cfg, cfg)); ax.legend(fontsize=7)
    fig.suptitle(title); fig.tight_layout(); fig.savefig(os.path.join(OUTF, fname), bbox_inches="tight"); plt.close(fig)
    print(f"wrote {fname}")


def sweep_section(md, matrix, xcol, title):
    df = load(matrix)
    if df.empty:
        return None
    df.to_csv(os.path.join(OUTA, f"sat_{matrix}.csv"), index=False)
    md.append(f"\n## {title} ({matrix})\n")
    for cfg in [c for c in FORM if c in df.config.unique()]:
        piv = df[(df.config == cfg) & (df.method.isin(PROP))].pivot_table(index=xcol, columns="method", values="coverage")
        md.append(f"{FORM.get(cfg,cfg)} coverage vs {xcol}:\n```\n" + piv.round(3).to_string() + "\n```\n")
    return df


def main():
    md = ["# Saturation campaign — analysis\n"]

    # alpha
    df = sweep_section(md, "alpha", "alpha", "Tolerance sweep (1-alpha)")
    if df is not None:
        _line(df, "alpha", "Coverage vs tolerance parameter alpha (target 0.95)", "sat_alpha.pdf", target_from="fixed")

    # beta calibration
    df = sweep_section(md, "beta", "beta", "Target-confidence sweep / calibration (1-beta)")
    if df is not None:
        # calibration: coverage vs nominal target = 1-beta
        df2 = df.copy(); df2["nominal"] = 1 - df2["beta"]
        cfgs = [c for c in FORM if c in df2.config.unique()]
        fig, axes = plt.subplots(1, len(cfgs), figsize=(4.3 * len(cfgs), 4), squeeze=False)
        for ax, cfg in zip(axes[0], cfgs):
            g = df2[df2.config == cfg]
            for m in PROP:
                gm = g[g.method == m].sort_values("nominal")
                if not gm.empty:
                    ax.plot(gm.nominal, gm.coverage, marker="o", ms=4, label=m)
            lim = [g.nominal.min(), 1.0]
            ax.plot(lim, lim, "k--", lw=1, label="nominal = empirical")
            ax.set_xlabel("nominal target 1-beta"); ax.set_ylabel("empirical coverage")
            ax.set_title(FORM.get(cfg, cfg)); ax.legend(fontsize=7)
        fig.suptitle("Calibration: empirical coverage vs nominal 1-beta (points on/above dashed = valid/conservative)")
        fig.tight_layout(); fig.savefig(os.path.join(OUTF, "sat_beta_calibration.pdf"), bbox_inches="tight"); plt.close(fig)
        print("wrote sat_beta_calibration.pdf")

    # tails
    df = sweep_section(md, "tails", "df", "Heavy-tail sweep (multivariate-t degrees of freedom)")
    if df is not None:
        _line(df, "df", "Coverage vs t degrees of freedom (smaller df = heavier tails; target 0.95)", "sat_tails.pdf", target_from="fixed")

    # corr
    df = sweep_section(md, "corr", "corr", "Correlation sweep (Gaussian coordinate correlation)")
    if df is not None:
        _line(df[df.d == 10], "corr", "Coverage vs correlation, d=10 (target 0.95)", "sat_corr.pdf", target_from="fixed")

    # foldsdeep: coverage vs K per scheme
    fd = load("foldsdeep")
    if not fd.empty:
        fd.to_csv(os.path.join(OUTA, "sat_foldsdeep.csv"), index=False)
        md.append("\n## Deep fold sweep: coverage vs K (foldsdeep)\n")
        # parse K from method names CV{K}/BS{K}/Sec{K}
        def parse(m):
            for s in ("CV", "BS", "Sec"):
                if m.startswith(s) and m[len(s):].isdigit():
                    return s, int(m[len(s):])
            return None, None
        fd[["scheme", "K"]] = fd["method"].apply(lambda m: pd.Series(parse(m)))
        so = fd[(fd.config == "paper_so") & fd.K.notna()]
        for cfg in [c for c in FORM if c in fd.config.unique()]:
            g = fd[(fd.config == cfg) & fd.K.notna()]
            if g.empty:
                continue
            piv = g.pivot_table(index=["n"], columns=["scheme", "K"], values="coverage")
            md.append(f"{FORM.get(cfg,cfg)} coverage by (scheme, K) and N:\n```\n" + piv.round(3).to_string() + "\n```\n")
        # figure: coverage vs K for SO at largest N, per scheme, with proposed reference lines
        if not so.empty:
            nmax = so.n.max()
            g = so[so.n == nmax]
            fig, ax = plt.subplots(figsize=(6, 4.2))
            for scheme in ("CV", "BS", "Sec"):
                gs = g[g.scheme == scheme].sort_values("K")
                if not gs.empty:
                    ax.plot(gs.K, gs.coverage, marker="o", label=scheme)
            prop = fd[(fd.config == "paper_so") & (fd.n == nmax) & (fd.method.isin(PROP))]
            for m in PROP:
                v = prop[prop.method == m].coverage
                if len(v):
                    ax.axhline(float(v.iloc[0]), ls=":", lw=1, label=f"{m} (proposed)")
            ax.axhline(0.95, color="r", ls="--", lw=1, label="target")
            ax.set_xlabel("fold count K"); ax.set_ylabel("coverage")
            ax.set_title(f"SO N={int(nmax)}: coverage vs fold count K"); ax.legend(fontsize=7)
            fig.savefig(os.path.join(OUTF, "sat_foldsdeep_so.pdf"), bbox_inches="tight"); plt.close(fig)
            print("wrote sat_foldsdeep_so.pdf")

    # bigd: coverage vs N at extreme d
    bd = load("bigd")
    if not bd.empty:
        bd.to_csv(os.path.join(OUTA, "sat_bigd.csv"), index=False)
        md.append("\n## Extreme dimension (bigd): coverage vs N at d in {100,200}\n")
        for cfg in [c for c in FORM if c in bd.config.unique()]:
            piv = bd[(bd.config == cfg) & (bd.method == "NGS")].pivot_table(index="d", columns="n", values="coverage")
            md.append(f"{FORM.get(cfg,cfg)} NGS coverage:\n```\n" + piv.round(3).to_string() + "\n```\n")

    # budgetfull merged into best-n1 (Wasserstein d=20 + high-D SO/RO)
    bf = load("budgetfull")
    if not bf.empty:
        bf.to_csv(os.path.join(OUTA, "sat_budgetfull.csv"), index=False)
        md.append("\n## Extended best-n1 (budgetfull): feasible Phase-1 band\n")
        for (cfg, n, d), g in bf[bf.method == "NGS"].groupby(["config", "n", "d"]):
            band = sorted(g[g.coverage >= 0.95].split.tolist())
            md.append(f"- {FORM.get(cfg,cfg)} n={n} d={d}: NGS feasible band = "
                      f"{','.join(f'{int(s*100)}%' for s in band) or 'none'}\n")

    with open(os.path.join(OUTA, "SATURATION.md"), "w") as f:
        f.write("\n".join(md))
    print("wrote results/analysis/SATURATION.md")


if __name__ == "__main__":
    main()
