"""Analysis + figures for the comprehensive matrices (budget / ndgrid / folds).

    /tmp/gsv_venv/bin/python runners/make_comprehensive_analysis.py

Writes results/analysis/comp_*.csv, results/analysis/COMPREHENSIVE.md, and
results/figures/comp_*.pdf. Robust to partial data (skips absent matrices).
Target coverage = 1 - beta = 0.95.
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
TARGET = 0.95
FORM = {"paper_so": "SO", "paper_ro_ellipsoid": "RO", "paper_saa": "SAA",
        "paper_dro_wasserstein": "Wasserstein", "paper_fast": "FAST"}


def load(matrix):
    rows = []
    for f in sorted(glob.glob(os.path.join(EXP, matrix, "*_summary.json"))):
        j = json.load(open(f))
        for m, s in j["summaries"].items():
            if isinstance(s, dict) and "coverage" in s:
                rows.append({"config": j["config"], "n": j["n"], "split": j["split"], "d": j["d"],
                             "reps": j["reps"], "method": m, "coverage": s["coverage"],
                             "coverage_lo": s["coverage_lo"], "coverage_hi": s["coverage_hi"],
                             "mean_obj": s["mean_obj"], "mean_oracle_gap": s.get("mean_oracle_gap"),
                             "meets": s["coverage"] >= TARGET})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- #
def analyze_budget(md):
    df = load("budget")
    if df.empty:
        return
    df.to_csv(os.path.join(OUTA, "comp_budget.csv"), index=False)
    md.append("## Best Phase-1 fraction n1/n vs N and D (budget matrix)\n")
    md.append("CI-AWARE criterion: a fraction is FEASIBLE only if its 95% Wilson lower bound "
              "clears 0.95 (coverage_lo >= 0.95), not merely the point estimate — this avoids the "
              "selection optimism the audit flagged. The RECOMMENDED n1 is the CI-feasible fraction "
              "with the best objective (ties -> smallest n1). `borderline` counts fractions whose "
              "point coverage >= 0.95 but whose CI still crosses 0.95.\n")
    recs = []; n_border = 0
    for cfg in sorted(df.config.unique()):
        for method in ["UG", "NGS", "UNGS"]:
            sub = df[(df.config == cfg) & (df.method == method)]
            for (n, d), g in sub.groupby(["n", "d"]):
                g = g.sort_values("split")
                feas = g[g.coverage_lo >= TARGET]                  # CI-strong feasibility
                point_only = g[(g.coverage >= TARGET) & (g.coverage_lo < TARGET)]
                n_border += len(point_only)
                band = sorted(feas.split.tolist())
                if not feas.empty:
                    best = feas.sort_values(["mean_obj", "split"]).iloc[0]
                    rec_s = float(best.split); rec_cov = float(best.coverage); rec_lo = float(best.coverage_lo)
                else:
                    rec_s = np.nan; rec_cov = float(g.coverage.max()); rec_lo = float(g.coverage_lo.max())
                recs.append({"formulation": FORM.get(cfg, cfg), "method": method, "n": n, "d": d,
                             "feasible_band_ci": ",".join(f"{int(s*100)}%" for s in band) or "none",
                             "n_feasible_ci": len(band), "recommended_n1": rec_s,
                             "cov_at_rec": rec_cov, "cov_lo_at_rec": rec_lo,
                             "n_point_only_borderline": len(point_only)})
    md.append(f"Borderline (point>=0.95 but CI crosses) across budget cells: {n_border}.\n")
    rec_df = pd.DataFrame(recs)
    rec_df.to_csv(os.path.join(OUTA, "comp_budget_best_n1.csv"), index=False)
    # trend summary for UG and NGS
    for method in ["UG", "NGS"]:
        md.append(f"\n### Recommended n1 (best feasible), {method}\n")
        piv = rec_df[rec_df.method == method].pivot_table(
            index=["formulation", "d"], columns="n", values="recommended_n1")
        md.append("```\n" + piv.round(2).to_string() + "\n```\n")
    md.append("Reading: entries are the recommended Phase-1 fraction; `NaN` = no tested n1 met "
              "target (data too small). Trends to note: how the recommended fraction moves as N "
              "grows (more data -> smaller n1 can suffice) and as D grows (higher dim -> need more "
              "Phase-1 data to estimate the path).\n")

    # figures: coverage heatmap over (n1 x N) at d=10, per formulation & method=NGS
    _budget_heatmaps(df, "NGS", "d", 10, "n", "comp_budget_heat_n1_vs_N_d10.pdf",
                     "Coverage over (Phase-1 fraction x N) at d=10 (NGS)")
    _budget_heatmaps(df, "NGS", "n", 500, "d", "comp_budget_heat_n1_vs_D_n500.pdf",
                     "Coverage over (Phase-1 fraction x d) at N=500 (NGS)")


def _budget_heatmaps(df, method, fix_col, fix_val, var_col, fname, title):
    sub = df[(df.method == method) & (df[fix_col] == fix_val)]
    cfgs = [c for c in ["paper_so", "paper_ro_ellipsoid", "paper_saa", "paper_dro_wasserstein"]
            if c in sub.config.unique()]
    if not cfgs:
        return
    fig, axes = plt.subplots(1, len(cfgs), figsize=(4.2 * len(cfgs), 4), squeeze=False)
    for ax, cfg in zip(axes[0], cfgs):
        g = sub[sub.config == cfg]
        piv = g.pivot_table(index="split", columns=var_col, values="coverage")
        im = ax.imshow(piv.values, aspect="auto", origin="lower", vmin=0.5, vmax=1.0, cmap="RdYlGn")
        ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns)
        ax.set_yticks(range(len(piv.index))); ax.set_yticklabels([f"{int(s*100)}%" for s in piv.index])
        ax.set_xlabel(var_col.upper()); ax.set_ylabel("Phase-1 fraction"); ax.set_title(FORM.get(cfg, cfg))
        for (i, si) in enumerate(piv.index):
            for (jj, cj) in enumerate(piv.columns):
                v = piv.values[i, jj]
                if np.isfinite(v):
                    ax.text(jj, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                            color="black" if v >= TARGET else "white")
    fig.colorbar(im, ax=axes[0].tolist(), shrink=0.8, label="coverage")
    fig.suptitle(title + "  (green >= 0.95 target)")
    fig.savefig(os.path.join(OUTF, fname), bbox_inches="tight"); plt.close(fig)
    print(f"wrote {fname}")


# ----------------------------------------------------------------------------- #
def analyze_ndgrid(md):
    df = load("ndgrid")
    if df.empty:
        return
    df.to_csv(os.path.join(OUTA, "comp_ndgrid.csv"), index=False)
    md.append("\n## Feasibility & optimality vs N and D (ndgrid matrix, split=0.5)\n")
    for method in ["UG", "NGS"]:
        md.append(f"\n### {method} coverage over (N x D)\n")
        for cfg in sorted(df.config.unique()):
            g = df[(df.config == cfg) & (df.method == method)]
            if g.empty:
                continue
            piv = g.pivot_table(index="d", columns="n", values="coverage")
            md.append(f"{FORM.get(cfg,cfg)}:\n```\n" + piv.round(3).to_string() + "\n```\n")
    # heatmaps: coverage over (N x D) for NGS, per formulation
    cfgs = [c for c in ["paper_so", "paper_ro_ellipsoid", "paper_saa", "paper_dro_wasserstein", "paper_fast"]
            if c in df.config.unique()]
    fig, axes = plt.subplots(1, len(cfgs), figsize=(4.0 * len(cfgs), 4), squeeze=False)
    for ax, cfg in zip(axes[0], cfgs):
        g = df[(df.config == cfg) & (df.method == "NGS")]
        piv = g.pivot_table(index="d", columns="n", values="coverage")
        im = ax.imshow(piv.values, aspect="auto", origin="lower", vmin=0.5, vmax=1.0, cmap="RdYlGn")
        ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns, rotation=45)
        ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
        ax.set_xlabel("N"); ax.set_ylabel("d"); ax.set_title(FORM.get(cfg, cfg))
    fig.colorbar(im, ax=axes[0].tolist(), shrink=0.8, label="coverage (NGS)")
    fig.suptitle("Coverage over (N x d), NGS (green >= 0.95)")
    fig.savefig(os.path.join(OUTF, "comp_ndgrid_coverage_heat.pdf"), bbox_inches="tight"); plt.close(fig)
    print("wrote comp_ndgrid_coverage_heat.pdf")


# ----------------------------------------------------------------------------- #
def analyze_folds(md):
    df = load("folds")
    if df.empty:
        return
    df.to_csv(os.path.join(OUTA, "comp_folds.csv"), index=False)
    md.append("\n## CV / bootstrap at K in {3,5,10} vs proposed (folds matrix, split=0.7)\n")
    for cfg in sorted(df.config.unique()):
        sub = df[df.config == cfg]
        methods = [m for m in ["UG", "NGS", "UNGS", "CV3", "CV5", "CV10", "BS3", "BS5", "BS10",
                               "Sec3", "Sec5", "Sec10", "benchmark"] if m in sub.method.unique()]
        piv = sub.pivot_table(index="n", columns="method", values="coverage")[[m for m in methods if m in sub.method.unique()]]
        md.append(f"\n### {FORM.get(cfg,cfg)} coverage vs N (target 0.95)\n```\n" + piv.round(3).to_string() + "\n```\n")
    # figure: coverage vs N, CV/BS by fold-count vs proposed, per formulation
    cfgs = sorted(df.config.unique())
    fig, axes = plt.subplots(1, len(cfgs), figsize=(5.0 * len(cfgs), 4.2), squeeze=False)
    styles = {"UG": ("C0", "-", "o"), "NGS": ("C1", "-", "s"), "UNGS": ("C2", "-", "^"),
              "CV3": ("C5", ":", "v"), "CV5": ("C5", "--", "v"), "CV10": ("C5", "-", "v"),
              "BS3": ("C6", ":", "P"), "BS5": ("C6", "--", "P"), "BS10": ("C6", "-", "P")}
    for ax, cfg in zip(axes[0], cfgs):
        sub = df[df.config == cfg]
        for m, (col, ls, mk) in styles.items():
            g = sub[sub.method == m].sort_values("n")
            if g.empty:
                continue
            ax.plot(g.n, g.coverage, color=col, ls=ls, marker=mk, ms=4, label=m)
        ax.axhline(TARGET, ls="--", color="r", lw=1)
        ax.set_xlabel("N"); ax.set_ylabel("coverage"); ax.set_title(FORM.get(cfg, cfg)); ax.legend(fontsize=6, ncol=2)
    fig.suptitle("CV/BS at K=3(dotted)/5(dashed)/10(solid) vs proposed (target 0.95 dashed red)")
    fig.savefig(os.path.join(OUTF, "comp_folds_coverage.pdf"), bbox_inches="tight"); plt.close(fig)
    print("wrote comp_folds_coverage.pdf")


def main():
    md = ["# Comprehensive experiments — analysis\n",
          f"Config: alpha=0.10 (tolerance 90%), beta=0.05 (target coverage {TARGET}); Gaussian DGP.\n"]
    analyze_budget(md)
    analyze_ndgrid(md)
    analyze_folds(md)
    with open(os.path.join(OUTA, "COMPREHENSIVE.md"), "w") as f:
        f.write("\n".join(md))
    print("wrote results/analysis/COMPREHENSIVE.md")


if __name__ == "__main__":
    main()
