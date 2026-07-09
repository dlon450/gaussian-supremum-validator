"""Structured experiment configuration and a named registry.

Replaces the parameters that were hard-coded at module scope in each driver.
An :class:`ExperimentConfig` fully specifies one experiment; the registry holds
both **legacy** configs (the exact values the committed drivers used, for the
numerical-equivalence gate) and **paper** configs (the author-confirmed values
for the revised experiments).

Author-confirmed decisions baked into the ``paper_*`` configs (2026-07-07):
  * (alpha, beta) = (0.10, 0.05)  -> tolerance 1-alpha = 90%, target 1-beta = 95%
  * N_reps = 1000 (full); pilots override to a small value
  * baseline DGP = true multivariate Gaussian N(mu, Sigma)
  * robustness DGP = heavy-tailed multivariate t
  * featured DRO = Wasserstein (moment-DRO kept only as a degeneracy remark)
  * oracle = large-independent-sample (and exact closed-form under Gaussian)

This module is pure Python (no numpy import) so it can be parsed/validated
anywhere; grid materialization lives in ``gsv.paths`` and sampling in the DGP
layer.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any

__all__ = [
    "DGP", "MeshSpec", "SolverOpts", "ExperimentConfig",
    "REGISTRY", "get_config", "list_configs", "pilot",
]

FORMULATIONS = {"ro_ellipsoid", "dro_moment", "so", "saa", "dro_wasserstein", "fast"}
VALIDATORS = {"UG", "NGS", "UNGS", "NV", "CV", "BS", "Sectioning"}
RNG_MODES = {"stream", "legacy"}
DGP_KINDS = {"gaussian", "half_normal", "multivariate_t"}


@dataclass(frozen=True)
class DGP:
    """Data-generating process. ``params`` carries kind-specific settings.

    gaussian / half_normal: {mu_scale, var_diag, corr}; multivariate_t adds {df}.
    """
    kind: str
    params: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        assert self.kind in DGP_KINDS, f"unknown DGP kind {self.kind!r}"
        if self.kind == "multivariate_t":
            assert float(self.params.get("df", 0)) > 2, "multivariate_t needs df > 2 for finite variance"


@dataclass(frozen=True)
class MeshSpec:
    """Declarative spec for the conservativeness-parameter grid {s_1,...,s_p}.

    ``kind`` is interpreted by ``gsv.paths`` (some kinds are data-dependent, e.g.
    the RO order-statistic mesh). ``p`` is the target number of grid points.
    """
    kind: str
    p: int | None = None
    params: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        if self.p is not None:
            assert self.p >= 2, "mesh needs at least 2 points"


@dataclass(frozen=True)
class SolverOpts:
    output_flag: int = 0
    threads: int = 1          # one solver thread per worker (see plan §7)
    time_limit: float | None = None
    mip_gap: float | None = None


@dataclass(frozen=True)
class ExperimentConfig:
    name: str
    formulation: str
    validators: tuple[str, ...]
    benchmark: str | None          # formulation-specific calibration, e.g. "SCA", "chi2", "SO_all", "FAST", or None
    alpha: float                   # chance-constraint tolerance is 1 - alpha
    beta: float                    # target feasibility coverage is 1 - beta
    d: int
    b_factor: float                # b = b_factor * d
    c_value: float                 # c = c_value * ones(d)
    n_grid: tuple[int, ...]
    split_grid: tuple[float, ...]  # Phase-1 fractions n1/n to sweep
    mesh: MeshSpec
    n_reps: int
    base_seed: int
    rng_mode: str
    dgp: DGP
    solver: SolverOpts = field(default_factory=SolverOpts)
    eval_samples: int = 1_000_000  # large-sample oracle / feasibility evaluation
    ub: float = 1.0                # box upper bound on x (0 <= x <= ub)

    def validate(self) -> "ExperimentConfig":
        assert self.formulation in FORMULATIONS, f"unknown formulation {self.formulation!r}"
        assert set(self.validators) <= VALIDATORS, f"unknown validators {set(self.validators) - VALIDATORS}"
        assert 0.0 < self.alpha < 1.0, "alpha must be in (0,1)"
        assert 0.0 < self.beta < 1.0, "beta must be in (0,1)"
        assert self.d >= 1 and self.n_reps >= 1
        assert all(n >= 2 for n in self.n_grid), "each n must be >= 2"
        assert all(0.0 < f < 1.0 for f in self.split_grid), "split fractions must be in (0,1)"
        assert self.rng_mode in RNG_MODES, f"unknown rng_mode {self.rng_mode!r}"
        self.dgp.validate()
        self.mesh.validate()
        return self


# --------------------------------------------------------------------------- #
# Shared DGP presets
# --------------------------------------------------------------------------- #
# Half-normal moments (what the committed code actually samples): |N(0,1)| has
# mean sqrt(2/pi) and variance 1 - 2/pi. The Gaussian baseline reuses these first
# two moments so the threshold b = 0.833 d stays comparably calibrated.
_HALF_NORMAL = DGP("half_normal", {"mu_scale": 1.0, "var_diag": 1.0, "corr": 0.0})
_GAUSSIAN = DGP("gaussian", {"mu_scale": 0.7978845608, "var_diag": 0.3633802276, "corr": 0.0})
_MVT = DGP("multivariate_t", {"mu_scale": 0.7978845608, "var_diag": 0.3633802276, "corr": 0.0, "df": 4})


def _legacy(name, formulation, mesh, n_grid, benchmark=None, dim_validators=("Sectioning", "NGS", "UNGS", "NV", "UG")):
    """Legacy config: exact committed values (half-normal DGP, alpha=0.05, beta=0.10,
    N=100, global-RNG seed=2, split sweep 10%..90%)."""
    return ExperimentConfig(
        name=name, formulation=formulation, validators=tuple(dim_validators), benchmark=benchmark,
        alpha=0.05, beta=0.10, d=10, b_factor=0.833, c_value=-1.0,
        n_grid=tuple(n_grid), split_grid=tuple(round(i / 10, 1) for i in range(1, 10)),
        mesh=mesh, n_reps=100, base_seed=2, rng_mode="legacy", dgp=_HALF_NORMAL,
    ).validate()


def _paper(name, formulation, mesh, n_grid, benchmark=None, split_grid=None,
           validators=("UG", "NGS", "UNGS", "NV", "CV", "BS", "Sectioning"), dgp=_GAUSSIAN):
    """Paper config: author-confirmed values (Gaussian DGP, alpha=0.10, beta=0.05, N=1000)."""
    return ExperimentConfig(
        name=name, formulation=formulation, validators=tuple(validators), benchmark=benchmark,
        alpha=0.10, beta=0.05, d=10, b_factor=0.833, c_value=-1.0,
        n_grid=tuple(n_grid),
        split_grid=tuple(split_grid) if split_grid else tuple(round(i / 10, 1) for i in range(1, 10)),
        mesh=mesh, n_reps=1000, base_seed=2, rng_mode="stream", dgp=dgp,
    ).validate()


# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #
REGISTRY: dict[str, ExperimentConfig] = {}


def _register(cfg: ExperimentConfig) -> None:
    REGISTRY[cfg.name] = cfg


# Legacy configs (reproduce committed results; meshes exactly as committed).
_register(_legacy("legacy_ro_ellipsoid", "ro_ellipsoid",
                  MeshSpec("linear_scaled", p=25, params={"coeff": 45}), [100, 200, 300, 400, 500]))
_register(_legacy("legacy_dro_moment", "dro_moment",
                  MeshSpec("linear_scaled", p=25, params={"coeff": 20}), [100, 200]))
_register(_legacy("legacy_so", "so",
                  MeshSpec("so_counts"), [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
                  dim_validators=("CV", "BS")))
_register(_legacy("legacy_saa", "saa",
                  MeshSpec("saa_counts"), [500], dim_validators=("CV", "BS")))
_register(_legacy("legacy_dro_wasserstein", "dro_wasserstein",
                  MeshSpec("wasserstein_geom"), [20, 40, 60, 80, 100, 120, 140, 160, 180, 200],
                  dim_validators=("NGS", "UNGS", "NV", "UG")))

# Paper configs (author-confirmed; corrected meshes; Gaussian baseline).
_register(_paper("paper_ro_ellipsoid", "ro_ellipsoid",
                MeshSpec("ro_orderstat", p=25, params={"offset": 20}),
                [100, 200, 300, 400, 500, 600, 800, 1000, 1500, 2000], benchmark="SCA"))
_register(_paper("paper_fast", "fast",
                MeshSpec("fast_segment", p=11), [200, 500], benchmark="FAST",
                validators=("UG", "NGS", "UNGS", "NV")))
_register(_paper("paper_dro_wasserstein", "dro_wasserstein",
                MeshSpec("wasserstein_geom"), [100, 200, 300, 400, 500], benchmark=None))
_register(_paper("paper_so", "so",
                MeshSpec("so_counts"), [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
                benchmark="SO_all",
                split_grid=(0.75,) + tuple(round(i / 10, 1) for i in range(1, 10))))
_register(_paper("paper_saa", "saa",
                MeshSpec("saa_offset", p=25),   # step defaults to alpha/p so s spans up to alpha (all-constraints)
                [100, 150, 200, 250, 300, 350, 400, 450, 500], benchmark=None))
_register(_paper("paper_dro_moment", "dro_moment",
                MeshSpec("dro_chi2", p=25, params={"scale": 1.5}), [200, 500], benchmark="chi2",
                validators=("UG", "NGS", "UNGS", "NV")))

# Robustness variant (heavy-tailed): featured formulations under multivariate t.
_register(_paper("robust_ro_ellipsoid_t", "ro_ellipsoid",
                MeshSpec("ro_orderstat", p=25, params={"offset": 20}),
                [200, 500, 1000], benchmark="SCA", dgp=_MVT))
_register(_paper("robust_dro_wasserstein_t", "dro_wasserstein",
                MeshSpec("wasserstein_geom"), [100, 200, 300, 400, 500], dgp=_MVT))


def get_config(name: str) -> ExperimentConfig:
    if name not in REGISTRY:
        raise KeyError(f"unknown config {name!r}; available: {sorted(REGISTRY)}")
    return REGISTRY[name]


def list_configs() -> list[str]:
    return sorted(REGISTRY)


def pilot(name: str, n_reps: int = 20, n_grid: tuple[int, ...] | None = None) -> ExperimentConfig:
    """Return a small-scale copy of a config for fast pilot/dev runs."""
    cfg = get_config(name)
    return replace(cfg, name=f"pilot_{name}", n_reps=n_reps,
                   n_grid=tuple(n_grid) if n_grid else cfg.n_grid,
                   eval_samples=min(cfg.eval_samples, 100_000)).validate()
