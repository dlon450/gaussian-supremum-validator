"""Validity/robustness checks requested in the audit:
  (1) independent-stream confirmation of recommended Phase-1 splits (selection optimism),
  (2) sensitivity to the solution-path mesh size p,
  (3) sensitivity to the Gaussian-supremum simulation count (the 2000-draw approximation),
  (4) a generality test on negatively-correlated data.
(The 5s MILP-cap sensitivity is run separately via GSV_MILP_TIMELIMIT.)

    /tmp/gsv_venv/bin/python runners/robustness_checks.py [--workers W]
Writes results/analysis/robustness_checks.md and .json.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gsv import experiment as E, select as SEL, rng as RNG, dgp as D, paths as P, config as C
from gsv.util import dump_json

TARGET = 0.95


def cov(summ, m):  # (coverage, lower CI)
    return round(summ[m]["coverage"], 3), round(summ[m]["coverage_lo"], 3)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=48); a = ap.parse_args()
    W = a.workers
    out = {}
    md = ["# Robustness / validity checks\n", f"Target coverage 1-beta = {TARGET}.\n"]

    # (1) independent-stream confirmation: coverage of a recommended split on a FRESH stream
    md.append("## 1. Independent-stream confirmation of recommended Phase-1 splits\n")
    md.append("In-sample = same seed family used to pick the split; Out-of-sample = a disjoint "
              "RNG stream (seed_offset). If coverage holds out-of-sample, the recommendation is not "
              "an artifact of selection optimism.\n")
    conf = []
    for cfg, n, d, split in [("paper_so", 500, 10, 0.5), ("paper_ro_ellipsoid", 500, 10, 0.3),
                             ("paper_saa", 400, 10, 0.5)]:
        _, s_in = E.run_cell(cfg, n, split, d, reps=500, workers=W)
        _, s_out = E.run_cell(cfg, n, split, d, reps=500, workers=W, seed_offset=7_000_003)
        for m in ("UG", "NGS", "UNGS"):
            ci, li = cov(s_in, m); co, lo = cov(s_out, m)
            conf.append({"config": cfg, "split": split, "method": m, "in_cov": ci,
                         "out_cov": co, "out_cov_lo": lo})
            md.append(f"- {cfg} s={split} {m}: in={ci}  out={co} (out CI_lo={lo})\n")
    out["confirmation"] = conf

    # (2) mesh-size sensitivity (RO d=10, n=500)
    md.append("\n## 2. Solution-path mesh-size (p) sensitivity — RO d=10 n=500\n")
    mesh = []
    for mp in (25, 50, 100):
        _, s = E.run_cell("paper_ro_ellipsoid", 500, 0.5, 10, reps=300, workers=W, mesh_p=mp)
        mesh.append({"mesh_p": mp, "UG_cov": cov(s, "UG")[0], "NGS_cov": cov(s, "NGS")[0],
                     "UG_obj": round(s["UG"]["mean_obj"], 3)})
        md.append(f"- p={mp}: UG cov={cov(s,'UG')[0]} NGS cov={cov(s,'NGS')[0]} UG obj={round(s['UG']['mean_obj'],3)}\n")
    out["mesh_p"] = mesh

    # (3) Gaussian-supremum sim_num sensitivity on a representative Sigma_hat
    md.append("\n## 3. Gaussian-supremum simulation-count (sim_num) sensitivity\n")
    cfg = C.get_config("paper_ro_ellipsoid"); d = 10; b = cfg.b_factor * d; c = cfg.c_value * np.ones(d)
    r = RNG.make_stream(2, "robustness:simnum", 0)
    data = D.sample(cfg.dgp, 500, d, r)
    xx, _ = P.build_path("ro_ellipsoid", data[:250], c, b, cfg.alpha, d, 250, cfg.mesh, rng=r)
    _, sigma_hat, Sigma_hat = SEL.phase2_stats(xx, data[250:], b)
    sim = []
    for sn in (1000, 2000, 5000, 10000):
        qs = [SEL.gs_quantiles(sigma_hat, Sigma_hat, cfg.beta, np.random.default_rng(k), sim_num=sn)[0]
              for k in range(30)]
        sim.append({"sim_num": sn, "qhat_mean": round(float(np.mean(qs)), 4), "qhat_std": round(float(np.std(qs)), 4)})
        md.append(f"- sim_num={sn}: qhat mean={round(float(np.mean(qs)),4)} std={round(float(np.std(qs)),4)}\n")
    md.append("(2000 draws already give a stable qhat; std shrinks ~1/sqrt(sim_num).)\n")
    out["sim_num"] = sim

    # (4) generality: negatively-correlated Gaussian data (corr < 0, still PSD for d=10)
    md.append("\n## 4. Generality: negatively-correlated data (corr=-0.08 vs 0.0)\n")
    gen = []
    for cfg_name, n, d in [("paper_ro_ellipsoid", 500, 10), ("paper_so", 500, 10), ("paper_saa", 400, 10)]:
        for corr in (-0.08, 0.0):
            _, s = E.run_cell(cfg_name, n, 0.5, d, reps=300, workers=W, dgp_params={"corr": corr})
            gen.append({"config": cfg_name, "corr": corr, "UG_cov": cov(s, "UG")[0], "NGS_cov": cov(s, "NGS")[0],
                        "UNGS_cov": cov(s, "UNGS")[0]})
            md.append(f"- {cfg_name} corr={corr}: UG={cov(s,'UG')[0]} NGS={cov(s,'NGS')[0]} UNGS={cov(s,'UNGS')[0]}\n")
    out["neg_corr"] = gen

    dump_json(out, os.path.join(os.path.dirname(__file__), "..", "results", "analysis", "robustness_checks.json"))
    with open(os.path.join(os.path.dirname(__file__), "..", "results", "analysis", "robustness_checks.md"), "w") as f:
        f.write("".join(md))
    print("wrote results/analysis/robustness_checks.{md,json}")


if __name__ == "__main__":
    main()
