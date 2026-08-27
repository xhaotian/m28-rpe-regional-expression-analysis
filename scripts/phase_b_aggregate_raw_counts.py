#!/usr/bin/env python3
"""Aggregate audited raw UMI counts to donor-by-region pseudobulk samples."""

import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("M28_INPUT_DIR", PROJECT_ROOT / "inputs"))
RESULTS_ROOT = Path(os.environ.get("M28_RESULTS_DIR", PROJECT_ROOT / "results"))
H5AD = Path(os.environ.get("M28_H5AD", INPUT_ROOT / "rpe_choroid.h5ad"))
OUT = RESULTS_ROOT / "01_DONOR_ROBUSTNESS"
CHUNK = 500


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(H5AD, backed="r")
    if a.raw is None:
        raise RuntimeError("Audited .raw matrix is absent")

    counts = a.obs.groupby(["donor_id", "Region"], observed=True).size().unstack(fill_value=0)
    paired = sorted(counts.index[(counts.get("Macular", 0) > 0) & (counts.get("Peripheral", 0) > 0)])
    if len(paired) != 15:
        raise RuntimeError(f"Expected 15 paired donors, found {len(paired)}")

    samples = [(donor, region) for donor in paired for region in ("Peripheral", "Macular")]
    sample_names = [f"{donor}__{region}" for donor, region in samples]
    sample_index = {key: i for i, key in enumerate(samples)}
    obs_keys = list(zip(a.obs["donor_id"].astype(str), a.obs["Region"].astype(str)))
    group_codes = np.array([sample_index.get(key, -1) for key in obs_keys], dtype=np.int16)

    summed = np.zeros((len(samples), a.raw.n_vars), dtype=np.int64)
    for start in range(0, a.n_obs, CHUNK):
        stop = min(start + CHUNK, a.n_obs)
        block = a.raw.X[start:stop, :]
        if not sp.issparse(block):
            block = sp.csr_matrix(block)
        codes = group_codes[start:stop]
        for code in np.unique(codes[codes >= 0]):
            summed[code, :] += np.asarray(block[codes == code, :].sum(axis=0)).ravel().astype(np.int64)

    var = a.raw.var.copy()
    gene_ids = var.index.astype(str)
    gene_symbols = var["feature_name"].astype(str).to_numpy()
    out = pd.DataFrame(summed.T, columns=sample_names)
    out.insert(0, "gene", gene_symbols)
    out.insert(0, "gene_id", gene_ids)
    out.to_csv(OUT / "RAWCOUNT_PSEUDOBULK_COUNTS.tsv.gz", sep="\t", index=False, compression="gzip")

    rows = []
    for i, (donor, region) in enumerate(samples):
        mask = group_codes == i
        unique = a.obs.loc[mask]
        rows.append({
            "sample": sample_names[i],
            "donor_id": donor,
            "Region": region,
            "n_cells": int(mask.sum()),
            "library_size": int(summed[i, :].sum()),
            "Study": "|".join(sorted(unique["Study"].astype(str).unique())),
            "donor_age": "|".join(sorted(unique["donor_age"].astype(str).unique())),
            "sex": "|".join(sorted(unique["sex"].astype(str).unique())),
        })
    meta = pd.DataFrame(rows)
    meta.to_csv(OUT / "RAWCOUNT_PSEUDOBULK_SAMPLE_METADATA.tsv", sep="\t", index=False)

    if int(summed.sum()) != int(meta["library_size"].sum()):
        raise RuntimeError("Count conservation check failed")
    if (summed < 0).any():
        raise RuntimeError("Negative aggregated counts detected")

    qc = pd.DataFrame({
        "metric": ["paired_donors", "pseudobulk_samples", "genes", "total_aggregated_umis", "minimum_sample_library", "maximum_sample_library"],
        "value": [len(paired), len(samples), a.raw.n_vars, int(summed.sum()), int(meta.library_size.min()), int(meta.library_size.max())],
    })
    qc.to_csv(OUT / "RAWCOUNT_PSEUDOBULK_AGGREGATION_QC.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
