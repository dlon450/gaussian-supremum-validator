"""Serialization + provenance helpers for reproducible result artifacts.

``sanitize`` converts NaN/Inf (which Python's ``json`` emits as the non-standard
tokens ``NaN``/``Infinity``, rejected by strict JSON parsers) into ``null`` so the
result files are valid JSON. ``dump_json``/``dump_jsonl`` write sanitized, and
stamp a ``_provenance`` block (git commit, package versions, timestamp) so every
result file records exactly how it was produced.
"""
from __future__ import annotations

import json
import math
import os
import subprocess


def sanitize(obj):
    """Recursively replace NaN/Inf with None so the result is strict-JSON-valid."""
    if isinstance(obj, float):
        return None if (math.isnan(obj) or math.isinf(obj)) else obj
    if isinstance(obj, dict):
        return {k: sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [sanitize(v) for v in obj]
    return obj


def provenance(timestamp: str | None = None) -> dict:
    """Git commit + key package versions + timestamp for a result file."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    def _git(*a):
        try:
            return subprocess.check_output(["git", "-C", root, *a], stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return None
    vers = {}
    for pkg in ("numpy", "scipy", "cvxpy"):
        try:
            vers[pkg] = __import__(pkg).__version__
        except Exception:
            vers[pkg] = None
    try:
        import gurobipy
        vers["gurobipy"] = ".".join(map(str, gurobipy.gurobi.version()))
    except Exception:
        vers["gurobipy"] = None
    return {"git_commit": _git("rev-parse", "HEAD"), "git_dirty": bool(_git("status", "--porcelain")),
            "versions": vers, "timestamp": timestamp}


def dump_json(obj, path, timestamp: str | None = None, add_provenance: bool = True):
    """Write ``obj`` as strict JSON (NaN->null), with a ``_provenance`` block."""
    if add_provenance and isinstance(obj, dict) and "_provenance" not in obj:
        obj = {**obj, "_provenance": provenance(timestamp)}
    with open(path, "w") as f:
        json.dump(sanitize(obj), f, indent=2, allow_nan=False)


def dump_jsonl(rows, path):
    """Write an iterable of dict rows as strict JSON lines (NaN->null)."""
    with open(path, "w") as f:
        for r in rows:
            f.write(json.dumps(sanitize(r), allow_nan=False) + "\n")
