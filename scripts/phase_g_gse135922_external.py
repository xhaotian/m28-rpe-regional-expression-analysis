#!/usr/bin/env python3
"""Audit and test predefined M28 genes in GSE135922 RPE cells.

The GEO matrices contain author-provided normalized expression and the
author's final cluster labels.  No attempt is made to reconstruct raw counts.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("M28_INPUT_DIR", PROJECT_ROOT / "inputs"))
RESULTS_ROOT = Path(os.environ.get("M28_RESULTS_DIR", PROJECT_ROOT / "results"))
OUT = RESULTS_ROOT / "04_GSE135922_EXTERNAL"
DATA = Path(os.environ.get("M28_GSE135922_DIR", INPUT_ROOT / "GSE135922" / "processed_tsv"))
CANDIDATES = Path(os.environ.get("M28_DISCOVERY_RESULTS", INPUT_ROOT / "discovery_paired_results.csv"))
CORE5 = ["WFDC1", "PDE3B", "SULF1", "COL8A1", "ABCA1"]


def metadata(path: Path) -> dict:
    name = path.name
    donor = int(re.search(r"donor_(\d+)", name).group(1))
    region = "Macula" if "_macula_" in name else "Peripheral"
    enriched = "enriched" in name
    return {
        "file": name,
        "gsm": name.split("_")[0],
        "donor": donor,
        "region": region,
        "disease": "AMD" if donor in (3, 4) else "No_known_eye_disease",
        "experiment": 2 if enriched else 1,
        "enrichment": "CD31_enriched" if enriched else "Unselected",
        "rpe_cluster": 4 if enriched else 7,
        "rpe_annotation_source": (
            "Voigt_et_al_exp1_explicit_RPE_cluster7" if not enriched else
            "author_final_cluster4_marker_mapping_RPE65_BEST1_RLBP1_TTR"
        ),
    }


def analyse_set(means: pd.DataFrame, genes: list[str], threshold: int) -> pd.DataFrame:
    eligible = means.query("disease == 'No_known_eye_disease' and rpe_cells >= @threshold")
    paired = []
    for donor, block in eligible.groupby("donor"):
        if set(block.region) == {"Macula", "Peripheral"}:
            paired.append(donor)
    rows = []
    for gene in genes:
        diffs = []
        detect = 0
        exp_support = {}
        for donor in paired:
            b = eligible[eligible.donor == donor].set_index("region")
            m, p = float(b.loc["Macula", gene]), float(b.loc["Peripheral", gene])
            diffs.append(m - p)
            detect += int((m > 0) or (p > 0))
            exp = int(b.experiment.iloc[0])
            exp_support.setdefault(exp, []).append(m - p)
        arr = np.asarray(diffs, dtype=float)
        rows.append({
            "gene": gene,
            "cell_threshold": threshold,
            "paired_n": len(arr),
            "detectable_pair_n": detect,
            "mean_macula_minus_peripheral": np.mean(arr) if len(arr) else np.nan,
            "median_macula_minus_peripheral": np.median(arr) if len(arr) else np.nan,
            "positive_donor_n": int((arr > 0).sum()),
            "support_fraction": float((arr > 0).mean()) if len(arr) else np.nan,
            "direction": "Macula_high" if len(arr) and np.mean(arr) > 0 else ("Peripheral_high" if len(arr) else "NA"),
            "experiment1_positive_n": int(sum(x > 0 for x in exp_support.get(1, []))),
            "experiment1_pair_n": len(exp_support.get(1, [])),
            "experiment2_positive_n": int(sum(x > 0 for x in exp_support.get(2, []))),
            "experiment2_pair_n": len(exp_support.get(2, [])),
            "analysis_method": "donor-region mean paired analysis of author-normalized expression",
            "inference_scope": "descriptive_direction_only" if len(arr) < 5 else "small_n_paired_direction",
        })
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cand = pd.read_csv(CANDIDATES)
    genes = cand.loc[(cand.paired_t_fdr < 0.05) & (cand.mean_mac_minus_per > 0), "gene"].astype(str).tolist()
    if len(genes) != 82:
        raise RuntimeError(f"Locked candidate set has {len(genes)} genes, expected 82")
    wanted = set(genes + ["final_cluster_labels", "barcode", "library"])
    sample_rows, mean_rows = [], []
    for path in sorted(DATA.glob("*.tsv.gz")):
        meta = metadata(path)
        df = pd.read_csv(path, sep=r"\s+", usecols=lambda x: x in wanted)
        rpe = df[df["final_cluster_labels"] == meta["rpe_cluster"]]
        row = dict(meta)
        row.update(total_cells=len(df), rpe_cells=len(rpe), usable_paired_rpe_cells=len(rpe))
        sample_rows.append(row)
        means = rpe.reindex(columns=genes).mean(axis=0).fillna(0).to_dict()
        mean_rows.append({**meta, "total_cells": len(df), "rpe_cells": len(rpe), **means})

    audit = pd.DataFrame(sample_rows)
    means = pd.DataFrame(mean_rows)
    audit["eligible_ge20"] = audit.rpe_cells >= 20
    audit["eligible_ge10"] = audit.rpe_cells >= 10
    audit.to_csv(OUT / "GSE135922_FEASIBILITY_AUDIT.tsv", sep="\t", index=False)
    means.to_csv(OUT / "GSE135922_RPE_DONOR_REGION_MEANS.tsv", sep="\t", index=False)

    core = pd.concat([analyse_set(means, CORE5, 20), analyse_set(means, CORE5, 10)], ignore_index=True)
    allc = pd.concat([analyse_set(means, genes, 20), analyse_set(means, genes, 10)], ignore_index=True)
    core.to_csv(OUT / "GSE135922_CORE5.tsv", sep="\t", index=False)
    allc.to_csv(OUT / "GSE135922_CANDIDATE82.tsv", sep="\t", index=False)

    healthy = audit[audit.disease == "No_known_eye_disease"]
    def paired_n(th):
        q = healthy[healthy.rpe_cells >= th]
        return sum(set(b.region) == {"Macula", "Peripheral"} for _, b in q.groupby("donor"))
    c20 = core[core.cell_threshold == 20]
    c10 = core[core.cell_threshold == 10]
    text = f"""# GSE135922 feasibility decision

## Decision

**LIMITED directional replication.** The public processed matrices include author-normalized single-cell expression and `final_cluster_labels`, but not raw UMI counts. Voigt et al. explicitly identify cluster 7 as RPE in experiment 1. In experiment 2, cluster 4 was mapped to residual RPE by its aggregate RPE65/BEST1/RLBP1/TTR signal; this mapping is marker-supported rather than an explicit cell-type field in the downloaded matrices.

Among the five donors without known eye disease, only **{paired_n(20)}** donors have at least 20 RPE cells in both regions; **{paired_n(10)}** donors meet the 10-cell sensitivity gate. Consequently, the 20-cell result is descriptive only, and the 10-cell result is a small-n direction check rather than a powered independent replication test.

## Core-gene direction

At the primary 20-cell threshold, {int((c20.direction == 'Macula_high').sum())}/5 core genes are macula-high. At the 10-cell sensitivity threshold, {int((c10.direction == 'Macula_high').sum())}/5 are macula-high. Gene-level values and donor support are in `GSE135922_CORE5.tsv`.

## Boundaries

- Donors 3 and 4 are excluded from primary validation because they had neovascular AMD.
- Experiments used different dissociation/enrichment protocols; experiment-specific direction columns are retained.
- The matrices are normalized expression, so the analysis is explicitly named **donor-region mean paired analysis**.
- Low and uneven RPE recovery, especially after CD31 enrichment, precludes edgeR count pseudobulk and strong null-hypothesis claims.
"""
    (OUT / "GSE135922_FEASIBILITY_DECISION.md").write_text(text)


if __name__ == "__main__":
    main()
