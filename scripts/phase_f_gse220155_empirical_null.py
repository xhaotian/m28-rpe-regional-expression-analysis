#!/usr/bin/env python3
"""Build an expression/detection-matched empirical null for GSE220155 concordance."""

from __future__ import annotations

import gzip
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260811
N_ITER = 10_000
PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("M28_INPUT_DIR", PROJECT_ROOT / "inputs"))
RESULTS_ROOT = Path(os.environ.get("M28_RESULTS_DIR", PROJECT_ROOT / "results"))
GSE = Path(os.environ.get("M28_GSE220155_DIR", INPUT_ROOT / "GSE220155"))
DISCOVERY_RESULTS = Path(os.environ.get("M28_DISCOVERY_RESULTS", INPUT_ROOT / "discovery_paired_results.csv"))
OUT = RESULTS_ROOT / "03_GSE220155_EMPIRICAL_BACKGROUND"
META = {"barcode", "region", "new_cluster_number", "library"}


def header(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as h:
        return h.readline().rstrip("\n").split(",")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pattern = re.compile(r"GSM\d+_donor_(\d+)_(macula|peripheral)_normalized\.csv\.gz$")
    files = []
    for path in sorted(GSE.glob("*_normalized.csv.gz")):
        if "geneactivity" in path.name:
            continue
        m = pattern.match(path.name)
        if m:
            files.append((path, f"donor_{m.group(1)}", m.group(2)))
    if len(files) != 8:
        raise RuntimeError(f"Expected 8 GEX files, found {len(files)}")

    common = None
    ordered = None
    for path, _, _ in files:
        genes = [x for x in header(path) if x not in META]
        common = set(genes) if common is None else common.intersection(genes)
        if ordered is None:
            ordered = genes
    comparable = [g for g in ordered if g in common]

    sums = pd.Series(0.0, index=comparable)
    detected = pd.Series(0, index=comparable, dtype=np.int64)
    total_cells = 0
    donor_means = {}
    file_rows = []
    for path, donor, region in files:
        df = pd.read_csv(path, usecols=comparable, dtype=np.float32)
        arr_sum = df.sum(axis=0).astype(float)
        arr_detect = (df > 0).sum(axis=0).astype(np.int64)
        sums += arr_sum
        detected += arr_detect
        total_cells += len(df)
        donor_means[(donor, region)] = arr_sum / len(df)
        file_rows.append({"file": path.name, "donor": donor, "region": region, "n_cells": len(df), "gene_n": len(comparable)})
        del df
    pd.DataFrame(file_rows).to_csv(OUT / "GSE220155_GEX_FILE_AUDIT.tsv", sep="\t", index=False)

    bg = pd.DataFrame({"gene": comparable, "external_mean_expression": sums.to_numpy() / total_cells,
                       "external_detection_rate": detected.to_numpy() / total_cells})
    for donor in sorted({d for _, d, _ in files}):
        bg[f"{donor}_difference"] = (donor_means[(donor, "macula")] - donor_means[(donor, "peripheral")]).reindex(comparable).to_numpy()
    diff_cols = [c for c in bg.columns if c.endswith("_difference")]
    bg["mean_macula_minus_peripheral"] = bg[diff_cols].mean(axis=1)
    bg["positive_donor_number"] = (bg[diff_cols] > 0).sum(axis=1)
    bg["mean_direction_positive"] = bg.mean_macula_minus_peripheral > 0
    bg["donor_support_ge3"] = bg.positive_donor_number >= 3
    bg["joint_positive_and_ge3"] = bg.mean_direction_positive & bg.donor_support_ge3
    bg.to_csv(OUT / "GSE220155_ALL_GENE_BACKGROUND.tsv", sep="\t", index=False)

    old = pd.read_csv(DISCOVERY_RESULTS)
    candidates = old.loc[(old.paired_t_fdr < .05) & (old.mean_mac_minus_per > 0), "gene"].tolist()
    actual = bg[bg.gene.isin(candidates)].copy()
    if len(actual) != 74:
        raise RuntimeError(f"Expected 74 detectable candidate genes, found {len(actual)}")

    pool = bg[(bg.external_detection_rate > 0) & ~bg.gene.isin(candidates)].copy()
    combined = pd.concat([pool, actual], ignore_index=True)
    combined["expression_bin"] = pd.qcut(combined.external_mean_expression.rank(method="first"), 5, labels=False)
    combined["detection_bin"] = pd.qcut(combined.external_detection_rate.rank(method="first"), 5, labels=False)
    actual_bins = combined[combined.gene.isin(actual.gene)].groupby(["expression_bin", "detection_bin"]).size()
    pool = combined[~combined.gene.isin(candidates)].copy()
    rng = np.random.default_rng(SEED)
    null_rows = []
    replacement_events = 0
    for iteration in range(N_ITER):
        selected = []
        for key, n in actual_bins.items():
            candidates_bin = pool[(pool.expression_bin == key[0]) & (pool.detection_bin == key[1])]
            replace = len(candidates_bin) < n
            replacement_events += int(replace)
            choose = rng.choice(candidates_bin.index.to_numpy(), size=int(n), replace=replace)
            selected.extend(choose.tolist())
        draw = pool.loc[selected]
        null_rows.append({
            "iteration": iteration + 1,
            "mean_direction_positive_n": int(draw.mean_direction_positive.sum()),
            "donor_support_ge3_n": int(draw.donor_support_ge3.sum()),
            "joint_positive_and_ge3_n": int(draw.joint_positive_and_ge3.sum()),
        })
    null = pd.DataFrame(null_rows)
    null.to_csv(OUT / "GSE220155_MATCHED_NULL_DISTRIBUTION.tsv", sep="\t", index=False)

    actual_stats = {
        "mean_direction_positive_n": int(actual.mean_direction_positive.sum()),
        "donor_support_ge3_n": int(actual.donor_support_ge3.sum()),
        "joint_positive_and_ge3_n": int(actual.joint_positive_and_ge3.sum()),
    }
    rows = []
    for metric, observed in actual_stats.items():
        values = null[metric].to_numpy()
        rows.append({
            "metric": metric, "candidate_detected_n": len(actual), "observed_count": observed,
            "null_mean": values.mean(), "null_sd": values.std(ddof=1),
            "empirical_percentile": 100 * (np.mean(values < observed) + 0.5 * np.mean(values == observed)),
            "empirical_p_upper": (1 + int((values >= observed).sum())) / (N_ITER + 1),
            "iterations": N_ITER, "matching_bins": "5 expression quantiles x 5 detection quantiles",
            "replacement_events_across_bin_draws": replacement_events,
        })
    enrich = pd.DataFrame(rows)
    enrich.to_csv(OUT / "GSE220155_EMPIRICAL_ENRICHMENT.tsv", sep="\t", index=False)
    actual[["gene", "external_mean_expression", "external_detection_rate", *diff_cols,
            "mean_macula_minus_peripheral", "positive_donor_number", "mean_direction_positive",
            "donor_support_ge3", "joint_positive_and_ge3"]].to_csv(
                OUT / "GSE220155_CANDIDATE82_DETECTED_RESULTS.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
