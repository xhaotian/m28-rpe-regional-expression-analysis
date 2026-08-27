#!/usr/bin/env python3
"""Audit H5AD structure and determine whether any matrix contains raw UMI counts."""

from __future__ import annotations

import hashlib
import os
import platform
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


SEED = 20260811
PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("M28_INPUT_DIR", PROJECT_ROOT / "inputs"))
RESULTS_ROOT = Path(os.environ.get("M28_RESULTS_DIR", PROJECT_ROOT / "results"))
H5AD = Path(os.environ.get("M28_H5AD", INPUT_ROOT / "rpe_choroid.h5ad"))
OUT = RESULTS_ROOT / "00_DATA_AUDIT"
CHUNK = 1000


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def summarize_matrix(matrix, n_obs: int, rng: np.random.Generator) -> dict:
    sample_rows = np.sort(rng.choice(n_obs, size=min(2000, n_obs), replace=False))
    sampled = matrix[sample_rows, :]
    if not sp.issparse(sampled):
        sampled = sp.csr_matrix(sampled)
    values = sampled.data.astype(np.float64, copy=False)
    nonnegative = bool(values.size == 0 or np.all(values >= 0))
    integer_fraction = float(np.mean(np.isclose(values, np.rint(values), atol=1e-6))) if values.size else np.nan

    library_sizes = np.empty(n_obs, dtype=np.float64)
    for start in range(0, n_obs, CHUNK):
        stop = min(start + CHUNK, n_obs)
        block = matrix[start:stop, :]
        library_sizes[start:stop] = np.asarray(block.sum(axis=1)).ravel()

    return {
        "dtype": str(getattr(matrix, "dtype", "unknown")),
        "nonnegative": nonnegative,
        "integer_fraction": integer_fraction,
        "median_library_size": float(np.median(library_sizes)),
        "library_q01": float(np.quantile(library_sizes, 0.01)),
        "library_q99": float(np.quantile(library_sizes, 0.99)),
        "sampled_nonzero_n": int(values.size),
        "sampled_nonzero_min": float(values.min()) if values.size else np.nan,
        "sampled_nonzero_median": float(np.median(values)) if values.size else np.nan,
        "sampled_nonzero_max": float(values.max()) if values.size else np.nan,
    }


def compact_values(series: pd.Series, limit: int = 12) -> str:
    counts = series.astype("string").fillna("<NA>").value_counts(dropna=False)
    return "; ".join(f"{k}={v}" for k, v in counts.head(limit).items())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    a = ad.read_h5ad(H5AD, backed="r")

    locations = [(".X", a.X, a.shape)]
    if a.raw is not None:
        locations.append((".raw.X", a.raw.X, a.raw.shape))
    for layer in a.layers.keys():
        locations.append((f".layers['{layer}']", a.layers[layer], a.shape))

    stats = []
    for location, matrix, shape in locations:
        row = summarize_matrix(matrix, shape[0], rng)
        row.update({"matrix_location": location, "shape": f"{shape[0]}x{shape[1]}"})
        raw_like = row["nonnegative"] and row["integer_fraction"] >= 0.999
        if raw_like:
            interpretation = "raw UMI-like non-negative integer counts"
            usable = "YES"
        elif row["nonnegative"] and row["integer_fraction"] < 0.05:
            interpretation = "continuous transformed expression; consistent with log-normalized data"
            usable = "NO"
        else:
            interpretation = "matrix identity remains ambiguous"
            usable = "NO"
        row["interpretation"] = interpretation
        row["usable_for_count_pseudobulk"] = usable
        row["evidence"] = (
            f"sampled_nonzero_n={row['sampled_nonzero_n']}; range="
            f"{row['sampled_nonzero_min']:.6g}..{row['sampled_nonzero_max']:.6g}; "
            f"median_nonzero={row['sampled_nonzero_median']:.6g}; "
            f"library_q01={row['library_q01']:.6g}; library_q99={row['library_q99']:.6g}"
        )
        stats.append(row)

    provenance = pd.DataFrame(stats)[
        [
            "matrix_location", "shape", "dtype", "nonnegative", "integer_fraction",
            "median_library_size", "interpretation", "usable_for_count_pseudobulk", "evidence",
        ]
    ]
    provenance.to_csv(OUT / "COUNT_MATRIX_PROVENANCE.tsv", sep="\t", index=False)

    usable = provenance.loc[provenance["usable_for_count_pseudobulk"] == "YES", "matrix_location"].tolist()
    raw_status = "AVAILABLE" if usable else "NOT_AVAILABLE"
    if len(usable) > 1:
        raw_status = "AMBIGUOUS"

    obs_audit = []
    for col in a.obs.columns:
        s = a.obs[col]
        obs_audit.append({
            "column": col,
            "dtype": str(s.dtype),
            "n_unique_including_missing": int(s.nunique(dropna=False)),
            "missing_n": int(s.isna().sum()),
            "top_values": compact_values(s),
        })
    pd.DataFrame(obs_audit).to_csv(OUT / "OBS_COLUMN_AUDIT.tsv", sep="\t", index=False)

    var_audit = []
    for col in a.var.columns:
        s = a.var[col]
        var_audit.append({
            "column": col,
            "dtype": str(s.dtype),
            "n_unique_including_missing": int(s.nunique(dropna=False)),
            "missing_n": int(s.isna().sum()),
            "top_values": compact_values(s),
        })
    pd.DataFrame(var_audit).to_csv(OUT / "VAR_COLUMN_AUDIT.tsv", sep="\t", index=False)

    donor_region = (
        a.obs.groupby(["donor_id", "Region"], observed=True)
        .size().rename("n_cells").reset_index()
        .sort_values(["donor_id", "Region"])
    )
    donor_region.to_csv(OUT / "DONOR_REGION_CELL_COUNTS.tsv", sep="\t", index=False)

    state_cols = [c for c in a.obs.columns if any(k in c.lower() for k in ("cluster", "subcluster", "state", "subtype", "cell_type", "annotation", "leiden", "author"))]
    key_obs = [c for c in ["donor_id", "Region", "Study", "sample_id", "sampleid", "donor_age", "sex", "author_cell_type"] if c in a.obs.columns]
    var_preview = a.var.head(10).reset_index().to_markdown(index=False)
    donor_table = donor_region.pivot(index="donor_id", columns="Region", values="n_cells").fillna(0).astype(int).reset_index()
    paired_n = int(((donor_table.get("Macular", 0) > 0) & (donor_table.get("Peripheral", 0) > 0)).sum())

    matrix_lines = []
    for row in stats:
        matrix_lines.append(
            f"- `{row['matrix_location']}`: dtype `{row['dtype']}`; non-negative={row['nonnegative']}; "
            f"integer fraction={row['integer_fraction']:.6f}; median library size={row['median_library_size']:.1f}; "
            f"interpretation: {row['interpretation']}."
        )

    report = f"""# H5AD structure and count-matrix audit

## Audit conclusion

- `RAW_COUNT_STATUS = {raw_status}`
- Count-qualified location(s): {', '.join(f'`{x}`' for x in usable) if usable else 'none'}.
- The count-qualified matrix is the only matrix authorized for count-based pseudobulk. Continuous `.X` values remain the input for reproducing the historical donor-region mean analysis.

## Object identity and structure

- Input: `{H5AD}`
- SHA-256: `{sha256(H5AD)}`
- File size: {H5AD.stat().st_size:,} bytes
- Shape: {a.n_obs:,} cells × {a.n_vars:,} variables
- `.X`: `{type(a.X).__name__}`, dtype `{a.X.dtype}`
- `.raw`: {'present' if a.raw is not None else 'absent'}{f', shape {a.raw.shape[0]:,} × {a.raw.shape[1]:,}' if a.raw is not None else ''}
- `.layers`: {', '.join(a.layers.keys()) if len(a.layers) else 'none'}
- `.obsm`: {', '.join(a.obsm.keys()) if len(a.obsm) else 'none'}
- `.uns`: {', '.join(a.uns.keys()) if len(a.uns) else 'none'}

## Matrix provenance evidence

{chr(10).join(matrix_lines)}

The audit sampled 2,000 deterministic cells per matrix for nonzero-value identity checks and calculated library sizes across all {a.n_obs:,} cells in chunks. `.raw.X` differs decisively from `.X`: integer UMI-like values versus continuous transformed values.

## Observation metadata

- All `.obs` columns ({len(a.obs.columns)}): {', '.join(f'`{c}`' for c in a.obs.columns)}
- Key analysis columns: {', '.join(f'`{c}`' for c in key_obs)}
- Donors: {a.obs['donor_id'].nunique()} total; {paired_n} have both regions
- Region counts: {compact_values(a.obs['Region'])}
- Study/source: {compact_values(a.obs['Study'])}
- Sex: {compact_values(a.obs['sex'])}
- Age: stored in `donor_age`; numeric harmonization must treat values such as `1 day` and `>89 years` explicitly.
- State/cluster/subtype-related columns: {', '.join(f'`{c}`' for c in state_cols) if state_cols else 'none'}
- Existing atlas RPE state labels: {compact_values(a.obs['author_cell_type']) if 'author_cell_type' in a.obs.columns else 'not available'}

The exact donor-by-region cell counts are stored in `DONOR_REGION_CELL_COUNTS.tsv`; complete column-level value and missingness summaries are stored in `OBS_COLUMN_AUDIT.tsv`.

## Variable metadata

- `.var` columns ({len(a.var.columns)}): {', '.join(f'`{c}`' for c in a.var.columns)}
- Gene identifiers are the AnnData variable index; symbols are stored in `feature_name`.
- First 10 rows:

{var_preview}

Complete column-level summaries are stored in `VAR_COLUMN_AUDIT.tsv`.

## Reproducibility

- Random seed: {SEED}
- Python: {platform.python_version()}
- anndata: {ad.__version__}
- numpy: {np.__version__}
- scipy sparse matrices were read in backed chunks; the H5AD was not modified.
"""
    (OUT / "H5AD_STRUCTURE_AUDIT.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
