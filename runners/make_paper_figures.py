"""Generate the figures the manuscript \\includegraphics but that are missing
from the repo, directly from results/experiments/<matrix>/*_summary.json.

    /tmp/gsv_venv/bin/python runners/make_paper_figures.py

Writes to results/figures/ (paper-named files):
  ro_ellipsoid_comparison.pdf  coverage & objective vs dimension d and vs n (RO), incl. SCA benchmark
  so_all.pdf                   SO: UG/NGS/UNGS/NV coverage & objective vs Phase-1 fraction n1/n
  saa_all.pdf                  SAA: same, vs n1/n
  so_delta_objective.pdf       SO: objective increase of CV/BS/Sectioning over UG (left) & NGS (right) vs n
  saa_delta_objective.pdf      SAA: same
  existing_coverage.pdf        SO/SAA/Wasserstein: coverage of ALL methods vs n at split=0.7 (the §6.5 comparison)

Every figure is backed by rows in results/figures/all_results.csv (written by make_figures.py).
The target line is 1 - beta = 0.95 (the author-confirmed paper config; see gsv/config.py).
"""
import os, glob, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXP = os.path.join(ROOT, "results", "experiments")
OUT = os.path.join(ROOT, "results", "figures")
os.makedirs(OUT, exist_ok=True)

TARGET = 0.95  # 1 - beta (paper config); the chance-constraint tolerance is 1 - alpha = 0.90


def _legend(ax, **kw):
    """Add a legend only if there are labeled artists (avoids empty-panel warnings)."""
    if ax.get_legend_handles_labels()[0]:
        ax.legend(**kw)

STYLE = {
    "UG":   ("C0", "o", "uni. Gaussian"),
    "NGS":  ("C1", "s", "norm. GS"),
    "UNGS": ("C2", "^", "unnorm. GS"),
    "NV":   ("C3", "x", "plain average"),
    "benchmark": ("C4", "d", "benchmark"),
    "CV":   ("C5", "v", "cross validation"),
    "BS":   ("C6", "P", "bootstrapping"),
    "Sectioning": ("C7", "*", "sectioning"),
}


def load_matrix(matrix):
    """Return list of (meta, method, summary_dict) for one matrix."""
    rows = []
    for f in sorted(glob.glob(os.path.join(EXP, matrix, "*_summary.json"))):
        j = json.load(open(f))
        meta = {"config": j["config"], "n": j["n"], "split": j["split"], "d": j["d"],
                "reps": j["reps"], "elapsed_s": j.get("elapsed_s")}
        for m, s in j["summaries"].items():
            if isinstance(s, dict) and "coverage" in s:
                rows.append((meta, m, s))
    return rows


def _by(rows, config=None, method=None, x="n"):
    """Filter rows and return sorted (xvals, dict-of-arrays) for a config/method."""
    sel = [(mt, s) for (mt, m, s) in rows
           if (config is None or mt["config"] == config) and (method is None or m == method)]
    sel.sort(key=lambda t: t[0][x])
    xs = np.array([mt[x] for mt, _ in sel], float)
    return sel, xs


# --------------------------------------------------------------------------- #
def fig_ro_comparison():
    dim = load_matrix("dim")
    nsw = [r for r in load_matrix("nsweep") if r[0]["config"] == "paper_ro_ellipsoid"]
    if not dim and not nsw:
        return
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    # Row 0: vs dimension d (dim matrix, RO n=500)
    for m in ["NV", "UG", "NGS", "UNGS", "benchmark"]:
        sel, xs = _by(dim, method=m, x="d")
        if len(xs):
            col, mk, lab = STYLE[m]
            ax[0, 0].plot(xs, [s["coverage"] for _, s in sel], marker=mk, color=col, label=lab)
            ax[0, 1].plot(xs, [s["mean_obj"] for _, s in sel], marker=mk, color=col, label=lab)
    ax[0, 0].axhline(TARGET, ls="--", color="k", lw=1, label=f"target {TARGET}")
    ax[0, 0].set(xscale="log", xlabel="dimension $d$", ylabel="empirical coverage",
                 title="RO coverage stays $\\geq$ target across $d$")
    ax[0, 1].set(xscale="log", xlabel="dimension $d$", ylabel="mean objective",
                 title="RO objective vs $d$ (lower is better)")
    _legend(ax[0, 0], fontsize=7); _legend(ax[0, 1], fontsize=7)
    # Row 1: vs data size n (nsweep matrix, RO d=10)
    for m in ["NV", "UG", "NGS", "UNGS", "benchmark"]:
        sel, xs = _by(nsw, method=m, x="n")
        if len(xs):
            col, mk, lab = STYLE[m]
            ax[1, 0].plot(xs, [s["coverage"] for _, s in sel], marker=mk, color=col, label=lab)
            lo = [s["coverage_lo"] for _, s in sel]; hi = [s["coverage_hi"] for _, s in sel]
            ax[1, 0].fill_between(xs, lo, hi, color=col, alpha=0.12)
            ax[1, 1].plot(xs, [s["mean_obj"] for _, s in sel], marker=mk, color=col, label=lab)
    ax[1, 0].axhline(TARGET, ls="--", color="k", lw=1)
    ax[1, 0].set(xlabel="data size $n$", ylabel="empirical coverage", title="RO coverage vs $n$")
    ax[1, 1].set(xlabel="data size $n$", ylabel="mean objective", title="RO objective vs $n$")
    _legend(ax[1, 0], fontsize=7); _legend(ax[1, 1], fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "ro_ellipsoid_comparison.pdf"), bbox_inches="tight")
    plt.close(fig); print("wrote ro_ellipsoid_comparison.pdf")


def fig_split(config, fname, title):
    """n1/n2 budgeting: coverage & objective of the proposed validators vs Phase-1 fraction."""
    rows = [r for r in load_matrix("split") if r[0]["config"] == config]
    if not rows:
        return
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for m in ["NV", "UG", "NGS", "UNGS"]:
        sel, xs = _by(rows, method=m, x="split")
        if not len(xs):
            continue
        col, mk, lab = STYLE[m]
        ax[0].plot(xs * 100, [s["coverage"] for _, s in sel], marker=mk, color=col, label=lab)
        lo = [s["coverage_lo"] for _, s in sel]; hi = [s["coverage_hi"] for _, s in sel]
        ax[0].fill_between(xs * 100, lo, hi, color=col, alpha=0.12)
        ax[1].plot(xs * 100, [s["mean_obj"] for _, s in sel], marker=mk, color=col, label=lab)
    ax[0].axhline(TARGET, ls="--", color="r", lw=1, label=f"target {TARGET}")
    ax[0].axvline(70, ls=":", color="gray", lw=1)
    ax[0].set(xlabel="% of data in Phase 1 ($n_1/n$)", ylabel="empirical coverage",
              title=f"{title}: coverage vs Phase-1 fraction")
    ax[1].set(xlabel="% of data in Phase 1 ($n_1/n$)", ylabel="mean objective",
              title=f"{title}: objective vs Phase-1 fraction")
    _legend(ax[0], fontsize=8); _legend(ax[1], fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, fname), bbox_inches="tight")
    plt.close(fig); print(f"wrote {fname}")


def fig_delta_objective(config, fname, title):
    """Objective increase of existing schemes over UG (left) and NGS (right) vs n."""
    rows = [r for r in load_matrix("existing") if r[0]["config"] == config]
    if not rows:
        return
    def obj_map(method):
        sel, xs = _by(rows, method=method, x="n")
        return {mt["n"]: s["mean_obj"] for mt, s in sel}
    base = {"UG": obj_map("UG"), "NGS": obj_map("NGS")}
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for j, ref in enumerate(["UG", "NGS"]):
        for m in ["CV", "BS", "Sectioning", "UNGS"]:
            om = obj_map(m)
            ns = sorted(set(om) & set(base[ref]))
            if not ns:
                continue
            col, mk, lab = STYLE[m]
            ax[j].plot(ns, [om[n] - base[ref][n] for n in ns], marker=mk, color=col, label=lab)
        ax[j].axhline(0, ls="--", color="k", lw=1)
        ax[j].set(xlabel="data size $n$", ylabel=f"mean objective $-$ {STYLE[ref][2]} objective",
                  title=f"{title}: extra objective over {STYLE[ref][2]}")
        _legend(ax[j], fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, fname), bbox_inches="tight")
    plt.close(fig); print(f"wrote {fname}")


def fig_existing_coverage():
    rows = load_matrix("existing")
    if not rows:
        return
    configs = sorted({mt["config"] for mt, _, _ in rows})
    fig, axes = plt.subplots(1, len(configs), figsize=(5.2 * len(configs), 4.2), squeeze=False)
    for ax, cfg in zip(axes[0], configs):
        sub = [r for r in rows if r[0]["config"] == cfg]
        for m in ["NV", "UG", "NGS", "UNGS", "CV", "BS", "Sectioning", "benchmark"]:
            sel, xs = _by(sub, method=m, x="n")
            if not len(xs):
                continue
            col, mk, lab = STYLE[m]
            ax.plot(xs, [s["coverage"] for _, s in sel], marker=mk, color=col, label=lab)
        ax.axhline(TARGET, ls="--", color="r", lw=1)
        ax.set(xlabel="data size $n$", ylabel="empirical coverage",
               title=f"{cfg} (split=0.7)"); _legend(ax, fontsize=6)
    fig.suptitle("Comparison with existing validation schemes: coverage vs $n$ (target 0.95 dashed)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "existing_coverage.pdf"), bbox_inches="tight")
    plt.close(fig); print("wrote existing_coverage.pdf")


def main():
    fig_ro_comparison()
    fig_split("paper_so", "so_all.pdf", "SO")
    fig_split("paper_saa", "saa_all.pdf", "SAA")
    fig_delta_objective("paper_so", "so_delta_objective.pdf", "SO")
    fig_delta_objective("paper_saa", "saa_delta_objective.pdf", "SAA")
    fig_existing_coverage()


if __name__ == "__main__":
    main()
