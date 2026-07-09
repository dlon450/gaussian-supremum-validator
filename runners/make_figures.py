"""Aggregate experiment summaries -> one machine-readable table + paper figures.

    /tmp/gsv_venv/bin/python runners/make_figures.py

Reads results/experiments/<matrix>/<cell>_summary.json (whatever exists) and writes
results/figures/all_results.csv plus, where the data is present:
  fig_dimension.pdf     UG coverage & objective advantage vs d (dimension-free headline)
  fig_coverage_vs_n.pdf coverage vs n per formulation, with the 1-beta target line
  fig_gap_vs_n.pdf      objective gap to the path oracle vs n
  fig_robust.pdf        coverage under heavy-tailed t (where GS can beat UG)
Every figure is backed by rows in all_results.csv.
"""
import os, sys, glob, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXP = os.path.join(ROOT, "results", "experiments")
OUT = os.path.join(ROOT, "results", "figures")
os.makedirs(OUT, exist_ok=True)

STYLE = {"UG": ("C0", "o", "UG (focal)"), "NGS": ("C1", "s", "NGS"), "UNGS": ("C2", "^", "UNGS"),
         "NV": ("C3", "x", "NV (naive)"), "benchmark": ("C4", "d", "benchmark"),
         "CV": ("C5", "v", "CV"), "BS": ("C6", "P", "bootstrap"), "Sectioning": ("C7", "*", "sectioning")}


def load():
    rows = []
    for f in glob.glob(os.path.join(EXP, "*", "*_summary.json")):
        j = json.load(open(f))
        matrix = os.path.basename(os.path.dirname(f))
        for m, s in j["summaries"].items():
            if not isinstance(s, dict) or "coverage" not in s:
                continue
            rows.append({"matrix": matrix, "config": j["config"], "formulation": None,
                         "n": j["n"], "split": j["split"], "d": j["d"], "reps": j["reps"],
                         "method": m, **{k: s.get(k) for k in
                         ["coverage", "coverage_lo", "coverage_hi", "meets_target", "mean_obj",
                          "obj_se", "mean_oracle_gap", "mean_excess_s", "failure_rate", "n_reps"]}})
    return pd.DataFrame(rows)


def fig_dimension(df):
    d = df[df.matrix == "dim"].sort_values("d")
    if d.empty:
        return
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for m in ["NV", "UG", "NGS", "UNGS"]:
        sub = d[d.method == m]
        if sub.empty:
            continue
        col, mk, lab = STYLE[m]
        ax[0].plot(sub.d, sub.coverage, marker=mk, color=col, label=lab)
    target = 0.95
    ax[0].axhline(target, ls="--", color="k", lw=1, label=f"target {target}")
    ax[0].set_xlabel("dimension d"); ax[0].set_ylabel("empirical coverage"); ax[0].set_title("Coverage stays at target across d")
    ax[0].set_xscale("log"); ax[0].legend(fontsize=8)
    # objective advantage benchmark - UG
    ug = d[d.method == "UG"].set_index("d")["mean_obj"]
    bench = d[d.method == "benchmark"].set_index("d")["mean_obj"]
    adv = (bench - ug).dropna()
    if not adv.empty:
        ax[1].plot(adv.index, adv.values, marker="o", color="C4")
        ax[1].set_xlabel("dimension d"); ax[1].set_ylabel("benchmark objective $-$ UG objective")
        ax[1].set_title("UG's advantage over the benchmark grows with d"); ax[1].set_xscale("log")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_dimension.pdf"), bbox_inches="tight"); plt.close(fig)
    print("wrote fig_dimension.pdf")


def _cov_vs_n(df, matrix, fname, title):
    d = df[df.matrix == matrix]
    if d.empty:
        return
    forms = sorted(d.config.unique())
    fig, axes = plt.subplots(1, len(forms), figsize=(5 * len(forms), 4), squeeze=False)
    for ax, cfg in zip(axes[0], forms):
        sub = d[d.config == cfg]
        for m in ["NV", "UG", "NGS", "UNGS", "benchmark", "CV", "BS", "Sectioning"]:
            ms = sub[sub.method == m].sort_values("n")
            if ms.empty:
                continue
            col, mk, lab = STYLE[m]
            ax.plot(ms.n, ms.coverage, marker=mk, color=col, label=lab)
            ax.fill_between(ms.n, ms.coverage_lo, ms.coverage_hi, color=col, alpha=0.12)
        ax.axhline(0.95, ls="--", color="k", lw=1)
        ax.set_xlabel("data size n"); ax.set_ylabel("coverage"); ax.set_title(cfg); ax.legend(fontsize=7)
    fig.suptitle(title); fig.tight_layout(); fig.savefig(os.path.join(OUT, fname), bbox_inches="tight"); plt.close(fig)
    print(f"wrote {fname}")


def fig_gap_vs_n(df):
    d = df[df.matrix == "nsweep"]
    if d.empty:
        return
    forms = sorted(d.config.unique())
    fig, axes = plt.subplots(1, len(forms), figsize=(5 * len(forms), 4), squeeze=False)
    for ax, cfg in zip(axes[0], forms):
        sub = d[d.config == cfg]
        for m in ["UG", "NGS", "UNGS"]:
            ms = sub[sub.method == m].sort_values("n")
            if ms.empty:
                continue
            col, mk, lab = STYLE[m]
            ax.plot(ms.n, ms.mean_oracle_gap, marker=mk, color=col, label=lab)
        ax.set_xlabel("data size n"); ax.set_ylabel("objective gap to path oracle"); ax.set_title(cfg); ax.legend(fontsize=7)
    fig.suptitle("Objective gap to oracle shrinks with n"); fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_gap_vs_n.pdf"), bbox_inches="tight"); plt.close(fig)
    print("wrote fig_gap_vs_n.pdf")


def main():
    df = load()
    if df.empty:
        print("no results yet"); return
    df.to_csv(os.path.join(OUT, "all_results.csv"), index=False)
    print(f"wrote all_results.csv ({len(df)} rows, matrices={sorted(df.matrix.unique())})")
    fig_dimension(df)
    _cov_vs_n(df, "nsweep", "fig_coverage_vs_n.pdf", "Coverage vs n (target 0.95 dashed)")
    _cov_vs_n(df, "robust", "fig_robust.pdf", "Coverage under heavy-tailed t (target 0.95 dashed)")
    fig_gap_vs_n(df)


if __name__ == "__main__":
    main()
