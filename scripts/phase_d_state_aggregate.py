#!/usr/bin/env python3
"""Audit atlas RPE states and aggregate raw counts plus normalized state means."""

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
DISCOVERY_RESULTS = Path(os.environ.get("M28_DISCOVERY_RESULTS", INPUT_ROOT / "discovery_paired_results.csv"))
OUT = RESULTS_ROOT / "02_RPE_STATE_DECOMPOSITION"
CHUNK = 500


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(H5AD, backed="r")
    states = sorted(a.obs["author_cell_type"].astype(str).unique())
    region_counts = a.obs.groupby(["donor_id", "Region"], observed=True).size().unstack(fill_value=0)
    donors = sorted(region_counts.index[(region_counts.Macular > 0) & (region_counts.Peripheral > 0)])
    groups = [(d, r, s) for d in donors for r in ("Peripheral", "Macular") for s in states]
    group_names = [f"{d}__{r}__{s.replace('+', 'pos')}" for d, r, s in groups]
    group_index = {g: i for i, g in enumerate(groups)}
    keys = list(zip(a.obs.donor_id.astype(str), a.obs.Region.astype(str), a.obs.author_cell_type.astype(str)))
    codes = np.array([group_index.get(k, -1) for k in keys], dtype=np.int16)

    raw_sum = np.zeros((len(groups), a.n_vars), dtype=np.int64)
    x_sum = np.zeros((len(groups), a.n_vars), dtype=np.float64)
    n_cells = np.bincount(codes[codes >= 0], minlength=len(groups)).astype(int)
    for start in range(0, a.n_obs, CHUNK):
        stop = min(start + CHUNK, a.n_obs)
        raw = a.raw.X[start:stop, :]
        x = a.X[start:stop, :]
        if not sp.issparse(raw): raw = sp.csr_matrix(raw)
        if not sp.issparse(x): x = sp.csr_matrix(x)
        cc = codes[start:stop]
        for code in np.unique(cc[cc >= 0]):
            take = cc == code
            raw_sum[code] += np.asarray(raw[take].sum(axis=0)).ravel().astype(np.int64)
            x_sum[code] += np.asarray(x[take].sum(axis=0)).ravel()

    symbols = a.raw.var["feature_name"].astype(str).to_numpy()
    gene_ids = a.raw.var.index.astype(str).to_numpy()
    raw_out = pd.DataFrame(raw_sum.T, columns=group_names)
    raw_out.insert(0, "gene", symbols)
    raw_out.insert(0, "gene_id", gene_ids)
    raw_out.to_csv(OUT / "STATE_RAWCOUNT_PSEUDOBULK_COUNTS.tsv.gz", sep="\t", index=False, compression="gzip")

    meta_rows = []
    for i, (d, r, s) in enumerate(groups):
        subset = a.obs.loc[codes == i]
        meta_rows.append({
            "sample": group_names[i], "donor_id": d, "Region": r, "state": s,
            "n_cells": int(n_cells[i]), "library_size": int(raw_sum[i].sum()),
            "Study": "|".join(sorted(subset.Study.astype(str).unique())) if len(subset) else "",
        })
    meta = pd.DataFrame(meta_rows)
    meta.to_csv(OUT / "STATE_PSEUDOBULK_SAMPLE_METADATA.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "STATE_CELL_COUNT_BY_DONOR_REGION.tsv", sep="\t", index=False)

    old = pd.read_csv(DISCOVERY_RESULTS)
    cand = old.loc[(old.paired_t_fdr < 0.05) & (old.mean_mac_minus_per > 0), "gene"].tolist()
    selected = list(dict.fromkeys(cand))
    symbol_to_idx = {g: i for i, g in enumerate(symbols)}
    missing = [g for g in selected if g not in symbol_to_idx]
    if missing:
        raise RuntimeError(f"Candidate symbols absent from H5AD: {missing}")
    idx = np.array([symbol_to_idx[g] for g in selected])
    means = np.full((len(groups), len(selected)), np.nan)
    present = n_cells > 0
    means[present] = x_sum[present][:, idx] / n_cells[present, None]
    long = pd.DataFrame(means, columns=selected)
    long.insert(0, "n_cells", n_cells)
    long.insert(0, "state", [x[2] for x in groups])
    long.insert(0, "Region", [x[1] for x in groups])
    long.insert(0, "donor_id", [x[0] for x in groups])
    long.to_csv(OUT / "STATE_CANDIDATE82_DONOR_REGION_MEANS.tsv.gz", sep="\t", index=False, compression="gzip")

    feasibility = []
    for state in states:
        sm = meta[meta.state == state].pivot(index="donor_id", columns="Region", values="n_cells").fillna(0)
        state_total = int((a.obs.author_cell_type.astype(str) == state).sum())
        state_obs = a.obs.loc[a.obs.author_cell_type.astype(str) == state]
        donor_fraction = state_obs.donor_id.value_counts(normalize=True)
        study_fraction = state_obs.Study.value_counts(normalize=True)
        for threshold in (10, 20):
            eligible = (sm.Macular >= threshold) & (sm.Peripheral >= threshold)
            feasibility.append({
                "state": state, "threshold_cells_per_donor_region": threshold,
                "paired_donor_n": int(eligible.sum()), "formal_testing_eligible": bool(eligible.sum() >= 8),
                "total_cells": state_total, "largest_donor_fraction": float(donor_fraction.max()),
                "largest_donor": donor_fraction.idxmax(), "Chen_lab_fraction": float(study_fraction.get("Chen_lab", 0)),
                "Sanes_GSE236566_fraction": float(study_fraction.get("Sanes_GSE236566", 0)),
            })
    feas = pd.DataFrame(feasibility)
    feas.to_csv(OUT / "STATE_FEASIBILITY.tsv", sep="\t", index=False)

    comp = meta.pivot_table(index=["donor_id", "Region"], columns="state", values="n_cells", fill_value=0)
    comp = comp.div(comp.sum(axis=1), axis=0).reset_index()
    comp.to_csv(OUT / "STATE_DONOR_COMPOSITION.tsv", sep="\t", index=False)

    primary = feas[(feas.threshold_cells_per_donor_region == 20) & feas.formal_testing_eligible].state.tolist()
    sensitivity = feas[(feas.threshold_cells_per_donor_region == 10) & feas.formal_testing_eligible].state.tolist()
    audit = f"""# RPE state-label audit

## Decision

The analysis uses the original atlas annotation column `author_cell_type`. It contains five biologically named RPE states that match the expected atlas vocabulary: {', '.join(f'`{s}`' for s in states)}. Region-blind reclustering was therefore not performed; doing so would discard a traceable author annotation without a demonstrated deficiency.

## Annotation inventory

- Cells: {a.n_obs:,} RPE cells.
- State counts: {a.obs.author_cell_type.astype(str).value_counts().to_dict()}.
- Primary within-state rule: at least 8 paired donors with ≥20 cells per donor×region×state.
- Sensitivity rule: at least 8 paired donors with ≥10 cells per donor×region×state.
- Primary-eligible states: {', '.join(primary) if primary else 'none'}.
- ≥10-cell sensitivity-eligible states: {', '.join(sensitivity) if sensitivity else 'none'}.

## Donor and source balance

`STATE_FEASIBILITY.tsv` records, for every state and threshold, paired-donor availability, the largest single-donor contribution, and study-source fractions. `STATE_DONOR_COMPOSITION.tsv` and `STATE_CELL_COUNT_BY_DONOR_REGION.tsv` preserve the complete donor composition structure.

The state labels are highly region-imbalanced for several states. In particular, sparse macular recovery limits formal paired within-state testing outside the dominant LGI1+ state. This is a property of the observed state-by-region sampling and is not repaired by imputing absent state cells.

## Consequence for interpretation

Formal within-state inference is restricted to states meeting the prespecified donor/cell gates. States failing the gate remain sample-structure observations. Full five-state donor-level decomposition is reported only where every state is observed in both regions; a prespecified transparent coarse sensitivity (`LGI1+` versus all other RPE states) is additionally used to avoid treating unobserved conditional state means as zero.
"""
    (OUT / "RPE_STATE_LABEL_AUDIT.md").write_text(audit, encoding="utf-8")


if __name__ == "__main__":
    main()
