"""Data-generating processes: true moments and sampling.

Supports the baseline Gaussian, the legacy half-normal (what the committed code
sampled, ``|N(0,1)|`` per coordinate), and a heavy-tailed multivariate-t for the
robustness experiment. All three are parameterized to share the same first two
moments (``mu``, ``Sigma``) so the threshold ``b`` stays comparably calibrated.

``sample`` accepts any RNG exposing ``normal`` / ``multivariate_normal`` /
``chisquare`` — i.e. both a ``numpy.random.Generator`` (stream mode) and the
global ``numpy.random`` module (legacy mode). For legacy bit-for-bit
reproduction, half-normal sampling uses ``rng.normal(size=(n, d))`` exactly as
the committed drivers did.
"""

from __future__ import annotations

import numpy as np

from .config import DGP

__all__ = ["moments", "sample"]

_HALF_NORMAL_MEAN = np.sqrt(2.0 / np.pi)
_HALF_NORMAL_VAR = 1.0 - 2.0 / np.pi


def moments(dgp: DGP, d: int) -> tuple[np.ndarray, np.ndarray]:
    """Return the true mean vector (d,) and covariance matrix (d, d)."""
    if dgp.kind == "half_normal":
        mu = np.full(d, _HALF_NORMAL_MEAN)
        Sigma = _HALF_NORMAL_VAR * np.eye(d)
        return mu, Sigma
    mu_scale = float(dgp.params.get("mu_scale", 0.0))
    var_diag = float(dgp.params.get("var_diag", 1.0))
    corr = float(dgp.params.get("corr", 0.0))
    mu = np.full(d, mu_scale)
    Sigma = var_diag * ((1.0 - corr) * np.eye(d) + corr * np.ones((d, d)))
    return mu, Sigma


def sample(dgp: DGP, n: int, d: int, rng) -> np.ndarray:
    """Draw ``n`` i.i.d. observations of shape (n, d) from ``dgp``."""
    if dgp.kind == "half_normal":
        # matches committed drivers: abs(np.random.normal(size=(n, d)))
        return np.abs(rng.normal(size=(n, d)))

    mu, Sigma = moments(dgp, d)
    if dgp.kind == "gaussian":
        return rng.multivariate_normal(mu, Sigma, size=n)

    if dgp.kind == "multivariate_t":
        df = float(dgp.params["df"])
        # multivariate t with mean mu and covariance Sigma:
        # X = mu + Z / sqrt(W/df), Z ~ N(0, Sigma * (df-2)/df), W ~ chi2(df)
        scale = Sigma * (df - 2.0) / df
        z = rng.multivariate_normal(np.zeros(d), scale, size=n)
        w = rng.chisquare(df, size=n)
        return mu + z / np.sqrt(w / df)[:, None]

    raise ValueError(f"unknown DGP kind {dgp.kind!r}")
