"""Centralized random-number generation for reproducible (and parallel) experiments.

Two backends are provided, deliberately kept separate:

* **Stream RNG** (:func:`make_stream`) — the default for *new* experiments. Each
  Monte-Carlo replication gets its own ``numpy.random.Generator`` seeded from a
  ``SeedSequence`` derived deterministically from ``(base_seed, config_id,
  replication_id)``. Because each stream is derived from its *own* id (not from a
  shared spawn counter), the streams are independent **and** the result of any
  replication does not depend on how many workers run or in what order they
  finish — the property we need for reproducible parallelism.

* **Legacy RNG** (:func:`legacy_seed`) — reproduces the *committed* results, which
  were generated with NumPy's global legacy ``RandomState`` (``np.random.seed`` +
  ``np.random.normal/permutation/randint/multivariate_normal``). ``Generator``
  and ``RandomState`` produce different numbers from the same seed, so exact
  bit-for-bit reproduction of the legacy CSVs must use this backend. The existing
  drivers/validators already call the global functions, so legacy mode simply
  seeds the global state and runs them unchanged.

Do **not** mix the two within one experiment: a config declares its RNG mode
(see :mod:`gsv.config`).
"""

from __future__ import annotations

import hashlib

import numpy as np

__all__ = ["stable_hash", "make_seed_sequence", "make_stream", "legacy_seed"]


def stable_hash(text: str, bits: int = 63) -> int:
    """Deterministic non-negative integer hash of ``text``.

    Python's built-in ``hash`` is salted per process, so it cannot key a
    reproducible seed. We use a truncated SHA-256 instead, which is stable across
    processes and runs.
    """
    digest = hashlib.sha256(text.encode("utf-8")).digest()
    return int.from_bytes(digest, "big") % (1 << bits)


def make_seed_sequence(base_seed: int, config_id: str, replication_id: int) -> np.random.SeedSequence:
    """Build the ``SeedSequence`` for one (config, replication) cell.

    The entropy tuple ``[base_seed, stable_hash(config_id), replication_id]`` fully
    determines the stream, independently of any other cell — so streams are both
    reproducible and order-independent across parallel workers.
    """
    return np.random.SeedSequence([int(base_seed), stable_hash(config_id), int(replication_id)])


def make_stream(base_seed: int, config_id: str, replication_id: int) -> np.random.Generator:
    """Return the independent ``Generator`` for one (config, replication) cell."""
    return np.random.default_rng(make_seed_sequence(base_seed, config_id, replication_id))


def legacy_seed(seed: int) -> None:
    """Seed NumPy's global legacy RNG, matching how the committed drivers ran.

    Use only for legacy-equivalence runs that reuse the existing driver/validator
    code paths (which call ``np.random.*`` directly).
    """
    np.random.seed(int(seed))
