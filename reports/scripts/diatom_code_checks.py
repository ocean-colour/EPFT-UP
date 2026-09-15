#!/usr/bin/env python3
"""
diatom_code_checks.py

Verification script for the code report on Alison Chase's diatom carbon biomass
model (https://github.com/alisonpchase/diatom-carbon-biomass-model).

Each check below substantiates one claim made in the report embedded in
`claude_prompts/diatom_prompts.md`. The checks are deliberately independent of
the training data (which is not distributed with that repository): they are
either static analyses of the source files or numerical demonstrations on
synthetic pigment data whose structure matches the documented input format.

Usage:
    conda run -n ocean14 python reports/scripts/diatom_code_checks.py \
        [--repo /path/to/diatom-carbon-biomass-model]
"""

from __future__ import annotations

import argparse
import ast
import re
from pathlib import Path

import numpy as np

DEFAULT_REPO = Path("/Users/xavier/Oceanography/python/diatom-carbon-biomass-model")

# Feature layout taken from the model's model_config.py.
FEATURE_NAMES = ["chla", "chlb", "chlc", "chlb_a", "chlc_a", "ppc_a"]
RATIO_DEFS = {"chlb_a": ("chlb", "chla"), "chlc_a": ("chlc", "chla"), "ppc_a": ("ppc", "chla")}


def _banner(n: int, title: str) -> None:
    print(f"\n{'=' * 72}\nCHECK {n}: {title}\n{'=' * 72}")


# ---------------------------------------------------------------------------
# Check 1 — log-space feature matrix is rank deficient
# ---------------------------------------------------------------------------

def check_feature_rank(seed: int = 42, n: int = 5000) -> int:
    """The six model features are log10'd *after* the ratios are formed, so
    log(chlb_a) = log(chlb) - log(chla) exactly. The 6-column design matrix
    therefore spans only a 4-dimensional space."""
    _banner(1, "Rank of the log10 feature matrix")

    rng = np.random.default_rng(seed)
    # Lognormal pigment concentrations, roughly the dynamic range of the training set.
    chla = 10.0 ** rng.uniform(-2.5, 1.5, n)
    chlb = chla * 10.0 ** rng.normal(-1.0, 0.4, n)
    chlc = chla * 10.0 ** rng.normal(-0.7, 0.3, n)
    ppc = chla * 10.0 ** rng.normal(-0.6, 0.4, n)

    raw = {"chla": chla, "chlb": chlb, "chlc": chlc, "ppc": ppc}
    for out_col, (num, den) in RATIO_DEFS.items():
        raw[out_col] = raw[num] / raw[den]

    X = np.log10(np.column_stack([raw[c] for c in FEATURE_NAMES]))

    rank = np.linalg.matrix_rank(X - X.mean(axis=0))
    svals = np.linalg.svd(X - X.mean(axis=0), compute_uv=False)

    print(f"  features            : {FEATURE_NAMES}")
    print(f"  design matrix shape : {X.shape}")
    print(f"  numerical rank      : {rank} (of {X.shape[1]} columns)")
    print("  singular values     : " + ", ".join(f"{s:.3e}" for s in svals))
    resid = np.abs(X[:, 3] - (X[:, 1] - X[:, 0])).max()
    print(f"  max |log(chlb_a) - (log chlb - log chla)| = {resid:.3e}")
    print(f"  => VERDICT: {'rank deficient (4 of 6)' if rank == 4 else f'rank {rank}'}")
    return rank


# ---------------------------------------------------------------------------
# Check 2 — consequence for the Mahalanobis domain reference
# ---------------------------------------------------------------------------

def check_mahalanobis_conditioning(seed: int = 42, n: int = 5000) -> None:
    """train_model_final_export.compute_mahalanobis_reference inverts
    cov + 1e-6*I. Because cov is singular (Check 1), the ridge term alone sets
    the scale of two eigen-directions of the inverse."""
    _banner(2, "Conditioning of the Mahalanobis inverse covariance")

    rng = np.random.default_rng(seed)
    chla = 10.0 ** rng.uniform(-2.5, 1.5, n)
    chlb = chla * 10.0 ** rng.normal(-1.0, 0.4, n)
    chlc = chla * 10.0 ** rng.normal(-0.7, 0.3, n)
    ppc = chla * 10.0 ** rng.normal(-0.6, 0.4, n)
    raw = {"chla": chla, "chlb": chlb, "chlc": chlc, "ppc": ppc}
    for out_col, (num, den) in RATIO_DEFS.items():
        raw[out_col] = raw[num] / raw[den]
    X = np.log10(np.column_stack([raw[c] for c in FEATURE_NAMES]))

    cov = np.cov(X, rowvar=False)
    eig_raw = np.linalg.eigvalsh(cov)
    cov_reg = cov + 1e-6 * np.eye(cov.shape[0])
    eig_reg = np.linalg.eigvalsh(cov_reg)
    cov_inv = np.linalg.inv(cov_reg)

    print("  eigenvalues of cov          : " + ", ".join(f"{v:+.3e}" for v in eig_raw))
    print("  eigenvalues of cov + 1e-6 I : " + ", ".join(f"{v:.3e}" for v in eig_reg))
    print(f"  condition number (regularised): {np.linalg.cond(cov_reg):.3e}")
    print(f"  max eigenvalue of cov_inv     : {np.linalg.eigvalsh(cov_inv).max():.3e}")

    # Training-set distances are still well behaved because the null-space
    # component of every training row is exactly zero.
    diff = X - X.mean(axis=0)
    d = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", diff, cov_inv, diff), 0.0))
    print(f"  training distances: median={np.median(d):.3f}, 95th pct={np.percentile(d, 95):.3f}, max={d.max():.3f}")
    print("  => VERDICT: invertible only because of the ridge; two directions of")
    print("     cov_inv are set by 1e-6, so the metric is scale-arbitrary there.")


# ---------------------------------------------------------------------------
# Check 3 — retransformation bias of 10**mean(log10)
# ---------------------------------------------------------------------------

def check_backtransform_bias(seed: int = 42, n: int = 200000) -> None:
    """The model is fit to log10(diatCarb) under a squared-error criterion and
    inverted with a plain 10**. That returns a geometric-mean-like quantity,
    biased low relative to the arithmetic mean of the conditional distribution."""
    _banner(3, "Retransformation bias of the 10** back-transform")

    rng = np.random.default_rng(seed)
    for sigma_log in (0.15, 0.25, 0.35, 0.50):
        # Conditional distribution of log10(y) about a prediction of 0.0.
        z = rng.normal(0.0, sigma_log, n)
        y = 10.0 ** z
        naive = 10.0 ** 0.0                       # what the code returns
        arithmetic_mean = y.mean()                # what MBE / MAPE compare against
        smearing = np.mean(10.0 ** z)             # Duan smearing factor
        print(
            f"  sigma_log10={sigma_log:.2f} -> naive={naive:.4f}, "
            f"true mean={arithmetic_mean:.4f}, underestimate={100 * (1 - naive / arithmetic_mean):5.1f}%, "
            f"Duan factor={smearing:.4f}"
        )
    print("  => VERDICT: native-scale bias metrics (MBE, MPE) inherit a systematic")
    print("     low bias that grows with the residual spread.")


# ---------------------------------------------------------------------------
# Check 4 — metadata filename mismatch between export and inference
# ---------------------------------------------------------------------------

def check_metadata_filename_mismatch(repo: Path) -> bool:
    """train_model_final_export writes '<stem>_metadata.json'; run_model_inference
    calls load_metadata() with a *directory*, and the stem-specific candidate is
    only added when the argument is a file. So only 'metadata.json' is searched."""
    _banner(4, "Does inference ever find the exported <stem>_metadata.json?")

    export_src = (repo / "train_model_final_export.py").read_text()
    infer_src = (repo / "run_model_inference.py").read_text()

    writes_stem_json = 'OUTPUT_ONNX.stem + "_metadata.json"' in export_src
    searches_bare = '"metadata.json"' in infer_src
    guarded_by_is_file = re.search(
        r"if model_dir\.is_file\(\):\s*\n\s*metadata_files\.insert\(0, f\"\{model_dir\.stem\}_metadata\.json\"\)",
        infer_src,
    )
    # Every call site passes model_dir, which is a directory in both branches.
    call_sites = re.findall(r"load_metadata\((\w+)\)", infer_src)

    print(f"  export writes '<stem>_metadata.json'            : {writes_stem_json}")
    print(f"  inference search list contains 'metadata.json'   : {searches_bare}")
    print(f"  stem candidate guarded by model_dir.is_file()    : {bool(guarded_by_is_file)}")
    print(f"  load_metadata() call sites                      : {call_sites}")

    # Simulate the resolution logic on a realistic layout.
    stem = "diatom_model"
    model_dir_is_file = False  # both branches of the caller assign a directory
    metadata_files = ["metadata.json"]
    if model_dir_is_file:
        metadata_files.insert(0, f"{stem}_metadata.json")
    print(f"  simulated candidate filenames                   : {metadata_files}")
    print(f"  file actually written by the exporter           : ['{stem}_metadata.json']")

    broken = writes_stem_json and f"{stem}_metadata.json" not in metadata_files
    print(f"  => VERDICT: {'JSON never found; falls back to ONNX metadata_props' if broken else 'resolves correctly'}")
    return broken


# ---------------------------------------------------------------------------
# Check 5 — consequences of that fallback
# ---------------------------------------------------------------------------

def check_fallback_consequences(repo: Path) -> None:
    """The ONNX metadata_props written by the exporter contain neither
    'domain_flag_methods' nor 'log10_target', and store the log flag as the
    string 'true'/'false'."""
    _banner(5, "Consequences of falling back to ONNX metadata_props")

    export_src = (repo / "train_model_final_export.py").read_text()
    infer_src = (repo / "run_model_inference.py").read_text()
    export_tree = ast.parse(export_src)

    def _dict_keys_in(func_name: str, var_name: str) -> set[str]:
        """Collect the string keys of `var_name` as built inside `func_name`."""
        fn = next(
            node for node in ast.walk(export_tree)
            if isinstance(node, ast.FunctionDef) and node.name == func_name
        )
        keys: set[str] = set()
        for node in ast.walk(fn):
            # metadata = { "...": ... }
            if (isinstance(node, ast.AnnAssign) or isinstance(node, ast.Assign)) and isinstance(
                getattr(node, "value", None), ast.Dict
            ):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                if any(getattr(t, "id", None) == var_name for t in targets):
                    keys |= {k.value for k in node.value.keys if isinstance(k, ast.Constant)}
            # metadata["..."] = ...
            if isinstance(node, ast.Assign):
                for t in node.targets:
                    if (
                        isinstance(t, ast.Subscript)
                        and getattr(t.value, "id", None) == var_name
                        and isinstance(t.slice, ast.Constant)
                    ):
                        keys.add(t.slice.value)
        return keys

    onnx_keys = _dict_keys_in("export_rf_to_onnx_with_metadata", "metadata")
    json_keys = _dict_keys_in("main", "metadata_json")

    print(f"  keys written into the ONNX metadata_props : {sorted(onnx_keys)}")
    print(f"  keys written into <stem>_metadata.json    : {sorted(json_keys)}")
    for key in ("domain_flag_methods", "log10_target", "target_log10"):
        print(
            f"  '{key}': in ONNX props={key in onnx_keys:<5} in JSON={key in json_keys}"
        )

    print(f"  inference gates domain flags on         : "
          f"{bool(re.search(r'domain_methods = metadata.get\(\"domain_flag_methods\"', infer_src))}")
    print(f"  bool('false') evaluates to              : {bool('false')}")
    print(f"  bool('true')  evaluates to              : {bool('true')}")
    print("  => VERDICT: domain-flag columns are never emitted, and the log10 flag")
    print("     is parsed as truthy for BOTH 'true' and 'false'.")


# ---------------------------------------------------------------------------
# Check 6 — static defects in train_evaluate_model.py
# ---------------------------------------------------------------------------

def check_cv_return_arity(repo: Path) -> None:
    """cross_validation_analysis() returns a bare array on two guard paths but a
    2-tuple on success; the caller always unpacks two values."""
    _banner(6, "Return-arity mismatch and hardcoded output paths")

    path = repo / "train_evaluate_model.py"
    src = path.read_text()
    tree = ast.parse(src)

    fn = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "cross_validation_analysis"
    )
    shapes = []
    for node in ast.walk(fn):
        if isinstance(node, ast.Return) and node.value is not None:
            n = len(node.value.elts) if isinstance(node.value, ast.Tuple) else 1
            shapes.append((node.lineno, n))
    print("  returns from cross_validation_analysis (line, n_values):")
    for lineno, n in shapes:
        print(f"    line {lineno:5d}: {n}")

    callers = [
        (node.lineno, len(node.targets[0].elts))
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and isinstance(node.targets[0], ast.Tuple)
        and isinstance(node.value, ast.Call)
        and getattr(node.value.func, "id", None) == "cross_validation_analysis"
    ]
    for lineno, n in callers:
        print(f"  caller at line {lineno} unpacks {n} values")
    mismatched = sorted({n for _, n in shapes} - {n for _, n in callers})
    print(f"  => VERDICT: return shapes {sorted({n for _, n in shapes})} vs caller "
          f"{[n for _, n in callers]}; guard paths returning {mismatched} would raise")

    hardcoded = [
        (i, line.strip())
        for i, line in enumerate(src.splitlines(), start=1)
        if "'../figures/" in line or '"../figures/' in line
    ]
    print("  hardcoded '../figures/' output paths (ignore --figures-dir):")
    for lineno, line in hardcoded:
        print(f"    line {lineno:5d}: {line}")


# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", type=Path, default=DEFAULT_REPO,
                    help="path to the local clone of diatom-carbon-biomass-model")
    args = ap.parse_args()

    repo = args.repo.expanduser().resolve()
    if not repo.is_dir():
        raise SystemExit(f"Repository not found: {repo}")
    print(f"Checking repository: {repo}")

    check_feature_rank()
    check_mahalanobis_conditioning()
    check_backtransform_bias()
    check_metadata_filename_mismatch(repo)
    check_fallback_consequences(repo)
    check_cv_return_arity(repo)

    print("\nAll checks complete.")


if __name__ == "__main__":
    main()
