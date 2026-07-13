"""Experiment orchestrator: one replication -> metrics, parallel-ready.

`run_replication` is a top-level, picklable function that executes ONE Monte-Carlo
replication of one (config, n, split, d) cell and returns one raw record per
method. It is deterministic in its ``(config, n, split, d, rep)`` key via
``gsv.rng`` streams, so results are independent of worker count / completion
order (reproducible parallelism). `run_cell` runs many replications (serially or
via a process pool) and aggregates with `gsv.metrics`.

Two-phase validators (UG focal, NGS, UNGS, NV) + the formulation benchmark are
evaluated for every formulation. The existing schemes (CV, bootstrap, sectioning)
are evaluated for the simple-mesh formulations (SO, SAA, Wasserstein), matching
the paper's Section 6.4 comparison, at the same total data budget.
"""

from __future__ import annotations

import time
from dataclasses import replace
import numpy as np

from . import config as C, rng as RNG, dgp as D, oracle as O, paths as P, select as SEL, metrics as M
from . import validators as V

TWO_PHASE = ["NV", "UG", "NGS", "UNGS"]
EXISTING = ["CV", "BS", "Sectioning"]
EXISTING_FORMULATIONS = {"so", "saa", "dro_wasserstein"}


def _existing_grid(formulation, alpha):
    """(delta placeholder, cv_grid, bs_grid, solver) for the existing-scheme validators."""
    if formulation == "so":
        grid = lambda ne: np.arange(1, ne + 1)
        return grid, grid, V_solver_so
    if formulation == "saa":
        grid = lambda ne: np.arange(np.ceil(ne * (1 - alpha)).astype(int), ne + 1)
        return grid, grid, V_solver_saa
    if formulation == "dro_wasserstein":
        return None, None, V_solver_wass
    raise ValueError(formulation)


# module-level solver adapters (picklable) with the CCP signature solve(para,c,b,data,alpha,d,n)
def V_solver_so(para, c, b, data, alpha, d, n):
    from .formulations import CCP_SO
    return CCP_SO(para, c, b, data, alpha, d, n)

def V_solver_saa(para, c, b, data, alpha, d, n):
    from .formulations import CCP_SAA
    return CCP_SAA(para, c, b, data, alpha, d, n)

def V_solver_wass(para, c, b, data, alpha, d, n):
    from .formulations import CCP_DRO_wasserstein
    return CCP_DRO_wasserstein(para, c, b, data, alpha, d, n)


def run_replication(cfg_name, n, split, d, rep, include_existing=False, folds=None, mesh_p=None,
                    alpha=None, beta=None, dgp_params=None):
    """Execute one replication; return list of raw record dicts (one per method).

    ``folds`` (tuple of ints, e.g. (3, 5, 10)) sweeps the CV/bootstrap/sectioning
    fold-count K, recording methods ``CV{K}``/``BS{K}``/``Sec{K}`` with per-scheme
    deterministic seeding (each independent of the others). ``folds=None`` keeps the
    legacy single-K=10 behavior (methods CV/BS/Sectioning) used by the ``existing``
    matrix. ``mesh_p`` overrides the solution-path grid size (caps path solves for
    wide sweeps of expensive formulations). ``alpha``/``beta`` override the
    chance-constraint tolerance / target confidence; ``dgp_params`` (dict, may include
    ``kind``) overrides the data-generating process (e.g. ``{'df': 3}`` for heavier
    tails, ``{'corr': 0.5}`` for correlated coordinates). Overrides are keyed into the
    RNG stream so each swept setting is an independent, reproducible experiment."""
    from .config import DGP
    _overridden = (alpha is not None) or (beta is not None) or bool(dgp_params)
    cfg = C.get_config(cfg_name)
    if mesh_p is not None:
        cfg = replace(cfg, mesh=replace(cfg.mesh, p=mesh_p))
    if alpha is not None or beta is not None:
        cfg = replace(cfg, alpha=alpha if alpha is not None else cfg.alpha,
                      beta=beta if beta is not None else cfg.beta)
    if dgp_params:
        kind = dgp_params.get("kind", cfg.dgp.kind)
        params = {**cfg.dgp.params, **{k: v for k, v in dgp_params.items() if k != "kind"}}
        cfg = replace(cfg, dgp=DGP(kind, params))
    alpha, beta = cfg.alpha, cfg.beta
    b = cfg.b_factor * d
    c = cfg.c_value * np.ones(d)
    target = 1.0 - beta
    mu_true, sigma_true = D.moments(cfg.dgp, d)
    n1 = max(2, round(split * n)); n2 = n - n1
    # override tag keeps swept settings on distinct RNG streams (reproducible, independent)
    otag = f":a{alpha}:be{beta}:{cfg.dgp.kind}:{sorted(cfg.dgp.params.items())}" if _overridden else ""
    key = f"{cfg_name}:n{n}:s{split}:d{d}{otag}"
    r = RNG.make_stream(cfg.base_seed, key, rep)

    def true_feas(x):
        if cfg.dgp.kind == "gaussian":
            return O.gaussian_feasibility(x, mu_true, sigma_true, b)
        return O.large_sample_feasibility(x, eval_sample, b)

    records = []
    t0 = time.time()
    try:
      data = D.sample(cfg.dgp, n, d, r)
      # closed-form oracle for Gaussian; large independent sample otherwise (200k -> MC error ~1e-3)
      eval_sample = None if cfg.dgp.kind == "gaussian" else D.sample(cfg.dgp, min(cfg.eval_samples, 200_000), d, r)
      perm = r.permutation(n)
      phase1, phase2 = data[perm[:n1]], data[perm[n1:]]

      # ---- two-phase validators + benchmark ----
      xx, s_values = P.build_path(cfg.formulation, phase1, c, b, alpha, d, n1, cfg.mesh, rng=r)
      objectives = c @ xx
      P_hat, sigma_hat, Sigma_hat = SEL.phase2_stats(xx, phase2, b)
      picks = SEL.select_all(P_hat, sigma_hat, Sigma_hat, objectives, alpha=alpha, beta=beta, n2=n2, rng=r)
      feas_path = (O.gaussian_feasibility(xx, mu_true, sigma_true, b) if cfg.dgp.kind == "gaussian"
                   else O.large_sample_feasibility(xx, eval_sample, b))
      po = O.path_oracle(s_values, objectives, feas_path, target=1 - alpha)

      def record(method, x, s_sel):
          feas = float(true_feas(x)); obj = float(c @ x)
          rec = {"method": method, "config": cfg_name, "formulation": cfg.formulation,
                 "n": n, "n1": n1, "n2": n2, "d": d, "split": split, "dgp": cfg.dgp.kind, "rep": rep,
                 "feasible": 1.0 if feas >= 1 - alpha else 0.0, "true_feas": feas, "objective": obj,
                 "selected_s": s_sel}
          if po["any_feasible"]:
              rec["oracle_gap"] = obj - po["best_obj"]
              rec["rel_oracle_gap"] = (obj - po["best_obj"]) / (abs(po["best_obj"]) or 1.0)
              rec["excess_s"] = (s_sel - po["min_feasible_s"]) if np.isfinite(s_sel) else np.nan
          records.append(rec)

      for m in TWO_PHASE:
          idx = picks[m][0]
          record(m, xx[:, idx], float(s_values[idx]))

      if cfg.benchmark:
          bench = P.benchmark_solution(cfg.benchmark, c=c, b=b, alpha=alpha, d=d,
                                       mu_true=mu_true, sigma_true=sigma_true, full_data=data, n1=n1)
          if bench is not None and np.all(np.isfinite(bench)):
              record("benchmark", np.asarray(bench, float), np.nan)

      # ---- existing schemes (CV / bootstrap / sectioning) on simple-mesh formulations ----
      if include_existing and cfg.formulation in EXISTING_FORMULATIONS:
          cv_grid, bs_grid, solve = _existing_grid(cfg.formulation, alpha)
          delta = P.WASSERSTEIN_GRID if cfg.formulation == "dro_wasserstein" else np.arange(1, n + 1)
          if folds is None:
              # legacy path (single K=10, names CV/BS/Sectioning): one seed per rep, sequential.
              K = 10
              np.random.seed(RNG.stable_hash(key + f":{rep}", bits=31))
              schemes = [
                  ("CV", lambda K=K: V.cross_validation(solve, delta, c, b, alpha, n, d, data, beta, K, param_grid=cv_grid)),
                  ("BS", lambda K=K: V.bootstrapping(solve, delta, c, b, alpha, n, d, data, beta, K, param_grid=bs_grid)),
                  ("Sectioning", lambda K=K: V.sectioning(solve, (cv_grid(n1) if cv_grid else delta), c, b, alpha, n, d, data, beta, n1, n2, K)),
              ]
          else:
              # fold-sweep path: methods CV{K}/BS{K}/Sec{K}, each independently seeded
              # (result of one scheme/K is independent of which others are run).
              schemes = []
              for K in folds:
                  schemes += [
                      (f"CV{K}", lambda K=K: V.cross_validation(solve, delta, c, b, alpha, n, d, data, beta, K, param_grid=cv_grid)),
                      (f"BS{K}", lambda K=K: V.bootstrapping(solve, delta, c, b, alpha, n, d, data, beta, K, param_grid=bs_grid)),
                      (f"Sec{K}", lambda K=K: V.sectioning(solve, (cv_grid(n1) if cv_grid else delta), c, b, alpha, n, d, data, beta, n1, n2, K)),
                  ]
          for name, fn in schemes:
              if folds is not None:
                  np.random.seed(RNG.stable_hash(key + f":{rep}:{name}", bits=31))  # per-scheme reproducibility
              try:
                  record(name, fn(), np.nan)   # per-scheme failure is isolated, not silently dropped for others
              except Exception as e:
                  records.append({"method": name, "config": cfg_name, "formulation": cfg.formulation,
                                  "n": n, "d": d, "split": split, "rep": rep, "failed": 1.0, "error": str(e)[:150]})
    except Exception as e:
        # whole-replication failure (e.g., solver license/size limit): record, do not crash the run
        return [{"method": "__error__", "config": cfg_name, "formulation": cfg.formulation,
                 "n": n, "split": split, "d": d, "rep": rep, "failed": 1.0, "error": str(e)[:200]}]

    for rec in records:
        rec["runtime_s"] = time.time() - t0
    return records


def run_cell(cfg_name, n, split, d, reps, include_existing=False, workers=1, folds=None, mesh_p=None,
             alpha=None, beta=None, dgp_params=None):
    """Run ``reps`` replications of one cell; return (raw_records, summaries_dict).

    ``folds``/``mesh_p``/``alpha``/``beta``/``dgp_params`` are forwarded to
    :func:`run_replication` (fold-count sweep; path grid-size; tolerance/confidence/
    DGP overrides)."""
    all_rows = []
    args = (include_existing, folds, mesh_p, alpha, beta, dgp_params)
    if workers and workers > 1:
        import concurrent.futures as cf
        with cf.ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(run_replication, cfg_name, n, split, d, rep, *args): rep
                    for rep in range(reps)}
            for f in cf.as_completed(futs):
                try:
                    all_rows.extend(f.result())
                except Exception as e:  # worker died (should be rare; run_replication self-guards)
                    all_rows.append({"method": "__error__", "rep": futs[f], "failed": 1.0, "error": str(e)[:200]})
    else:
        for rep in range(reps):
            all_rows.extend(run_replication(cfg_name, n, split, d, rep, *args))

    cfg = C.get_config(cfg_name)
    target = 1.0 - (beta if beta is not None else cfg.beta)
    n_errors = sum(1 for r in all_rows if r["method"] == "__error__")
    method_rows = [r for r in all_rows if r["method"] != "__error__"]
    methods = sorted({r["method"] for r in method_rows},
                     key=lambda m: (TWO_PHASE + ["benchmark"] + EXISTING).index(m) if m in (TWO_PHASE + ["benchmark"] + EXISTING) else 99)
    summaries = {}
    for m in methods:
        rows = [r for r in method_rows if r["method"] == m and not r.get("failed")]
        if not rows:
            summaries[m] = {"method": m, "n_reps": 0, "n_failed": sum(1 for r in method_rows if r["method"] == m)}
            continue
        s = M.summarize_method(
            m, feasible_flags=[r["feasible"] for r in rows], objectives=[r["objective"] for r in rows],
            target=target, oracle_gaps=[r.get("oracle_gap", np.nan) for r in rows],
            rel_oracle_gaps=[r.get("rel_oracle_gap", np.nan) for r in rows],
            excess_s=[r.get("excess_s", np.nan) for r in rows],
            selected_s=[r.get("selected_s", np.nan) for r in rows]).to_dict()
        s["n_failed"] = sum(1 for r in method_rows if r["method"] == m and r.get("failed"))
        summaries[m] = s
    summaries["_n_replication_errors"] = n_errors
    return all_rows, summaries
