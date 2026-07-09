"""Statistical reporting for the evaluation criterion.

Everything here operates on per-replication arrays (one entry per Monte-Carlo
replication) and returns point estimates **with uncertainty**, as the brief
requires: every coverage estimate carries a Monte-Carlo confidence interval, and
objective/gap estimates carry standard errors. Common-random-number (paired)
differences are supported for like-for-like method comparison.

Primary criterion: empirical feasibility coverage vs. the target ``1 - beta``
(attain/near, not exceed). Secondary: lowest objective among near-target methods.
Plus objective gap and excess conservativeness vs. the path oracle.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
from scipy.stats import norm, t as student_t

__all__ = [
    "coverage_ci", "mean_se", "paired_diff", "MethodSummary", "summarize_method",
]


def coverage_ci(feasible_flags, level: float = 0.95, method: str = "wilson") -> tuple[float, float, float]:
    """Coverage point estimate and CI from 0/1 feasibility flags.

    Returns ``(p_hat, lo, hi)``. Wilson (default) is well-behaved near 0/1 and for
    small samples; ``method='normal'`` gives the Wald interval.
    """
    flags = np.asarray(feasible_flags, dtype=float)
    n = flags.size
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = float(flags.mean())
    z = float(norm.ppf(0.5 * (1.0 + level)))
    if method == "normal":
        half = z * np.sqrt(p * (1.0 - p) / n)
        return (p, max(0.0, p - half), min(1.0, p + half))
    # Wilson score interval
    denom = 1.0 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z / denom) * np.sqrt(p * (1.0 - p) / n + z * z / (4 * n * n))
    return (p, max(0.0, center - half), min(1.0, center + half))


def mean_se(values) -> tuple[float, float]:
    """Sample mean and standard error of the mean (ddof=1)."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    n = v.size
    if n == 0:
        return (np.nan, np.nan)
    if n == 1:
        return (float(v[0]), np.nan)
    return (float(v.mean()), float(v.std(ddof=1) / np.sqrt(n)))


def paired_diff(a, b, level: float = 0.95) -> dict:
    """Paired (common-random-number) mean difference ``a - b`` with a CI.

    Uses a Student-t interval (exact under normality, sensible for small N;
    ->normal as N grows). NaN pairs are dropped jointly.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    diff = a - b
    diff = diff[np.isfinite(diff)]
    n = diff.size
    if n < 2:
        return {"mean_diff": float(diff.mean()) if n else np.nan, "se": np.nan,
                "lo": np.nan, "hi": np.nan, "n": n}
    m = float(diff.mean())
    se = float(diff.std(ddof=1) / np.sqrt(n))
    tcrit = float(student_t.ppf(0.5 * (1.0 + level), df=n - 1))
    return {"mean_diff": m, "se": se, "lo": m - tcrit * se, "hi": m + tcrit * se, "n": n}


@dataclass
class MethodSummary:
    method: str
    n_reps: int
    target: float                 # 1 - beta
    coverage: float
    coverage_lo: float
    coverage_hi: float
    coverage_gap: float           # coverage - target (negative => under target)
    meets_target: bool            # point coverage >= target (paper's criterion; see coverage_lo/hi for uncertainty)
    significantly_below: bool     # coverage CI upper bound < target (coverage is statistically short of target)
    mean_obj: float
    obj_se: float
    mean_oracle_gap: float        # mean(selected_obj - path_oracle_obj)
    oracle_gap_se: float
    mean_rel_oracle_gap: float
    mean_excess_s: float          # mean(selected_s - min_feasible_s)
    selected_s_std: float
    failure_rate: float

    def to_dict(self) -> dict:
        return asdict(self)


def summarize_method(method: str, *, feasible_flags, objectives, target: float,
                     oracle_gaps=None, rel_oracle_gaps=None, excess_s=None,
                     selected_s=None, failures=None, level: float = 0.95) -> MethodSummary:
    """Aggregate one method's per-replication arrays into a reporting record."""
    flags = np.asarray(feasible_flags, dtype=float)
    n = flags.size
    p, lo, hi = coverage_ci(flags, level=level)
    mobj, seobj = mean_se(objectives)
    mgap, segap = mean_se(oracle_gaps) if oracle_gaps is not None else (np.nan, np.nan)
    mrel, _ = mean_se(rel_oracle_gaps) if rel_oracle_gaps is not None else (np.nan, np.nan)
    mexc, _ = mean_se(excess_s) if excess_s is not None else (np.nan, np.nan)
    s_std = float(np.nanstd(np.asarray(selected_s, dtype=float), ddof=1)) if selected_s is not None and n > 1 else np.nan
    fail = float(np.mean(np.asarray(failures, dtype=float))) if failures is not None else 0.0
    return MethodSummary(
        method=method, n_reps=n, target=target,
        coverage=p, coverage_lo=lo, coverage_hi=hi, coverage_gap=p - target,
        meets_target=bool(p >= target), significantly_below=bool(hi < target),
        mean_obj=mobj, obj_se=seobj,
        mean_oracle_gap=mgap, oracle_gap_se=segap, mean_rel_oracle_gap=mrel,
        mean_excess_s=mexc, selected_s_std=s_std, failure_rate=fail,
    )
