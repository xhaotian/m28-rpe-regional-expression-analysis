#!/usr/bin/env python3
"""Compare edgeR with historical P3 and calculate exact paired robust inference/sensitivities."""

from __future__ import annotations

from itertools import product
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


SEED = 20260811
BOOT_N = 10_000
PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("M28_INPUT_DIR", PROJECT_ROOT / "inputs"))
RESULTS_ROOT = Path(os.environ.get("M28_RESULTS_DIR", PROJECT_ROOT / "results"))
OUT = RESULTS_ROOT / "01_DONOR_ROBUSTNESS"
CORE5 = ["WFDC1", "PDE3B", "SULF1", "COL8A1", "ABCA1"]
H5AD = Path(os.environ.get("M28_H5AD", INPUT_ROOT / "rpe_choroid.h5ad"))
DISCOVERY_RESULTS = Path(os.environ.get("M28_DISCOVERY_RESULTS", INPUT_ROOT / "discovery_paired_results.csv"))


def bh(values: pd.Series) -> np.ndarray:
    arr = values.to_numpy(float)
    valid = np.isfinite(arr)
    out = np.full(len(arr), np.nan)
    out[valid] = multipletests(arr[valid], method="fdr_bh")[1]
    return out


def signflip_p(diff: np.ndarray) -> float:
    diff = np.asarray(diff, float)
    signs = np.array(list(product((-1.0, 1.0), repeat=len(diff))), dtype=np.float64)
    perm = signs @ diff / len(diff)
    return float(np.mean(np.abs(perm) >= abs(diff.mean()) - 1e-14))


def robust_row(gene: str, diff: np.ndarray, rng: np.random.Generator) -> dict:
    d = np.asarray(diff, float)
    n = len(d)
    mean = d.mean()
    sd = d.std(ddof=1) if n > 1 else np.nan
    se = sd / np.sqrt(n) if n > 1 else np.nan
    ci_half = stats.t.ppf(0.975, n - 1) * se if n > 1 else np.nan
    boot_idx = rng.integers(0, n, size=(BOOT_N, n))
    boot_mean = d[boot_idx].mean(axis=1)
    nonzero = d[d != 0]
    wilcox = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox", method="auto").pvalue if nonzero.size else 1.0
    sign_p = stats.binomtest(int((nonzero > 0).sum()), n=int(nonzero.size), p=0.5, alternative="two-sided").pvalue if nonzero.size else 1.0
    return {
        "gene": gene, "n_donors": n, "mean_difference": mean, "sd": sd, "se": se,
        "classical_ci_low": mean - ci_half, "classical_ci_high": mean + ci_half,
        "median_difference": np.median(d), "iqr": np.quantile(d, 0.75) - np.quantile(d, 0.25),
        "positive_donor_count": int((d > 0).sum()), "support_fraction": float((d > 0).mean()),
        "signflip_p": signflip_p(d), "wilcoxon_p": float(wilcox), "sign_test_p": float(sign_p),
        "bootstrap_ci_low": float(np.quantile(boot_mean, 0.025)),
        "bootstrap_ci_high": float(np.quantile(boot_mean, 0.975)),
    }


def inference_table(diff_wide: pd.DataFrame, genes: list[str], seed_offset: int = 0) -> pd.DataFrame:
    rows = []
    for i, gene in enumerate(genes):
        if gene not in diff_wide.index:
            continue
        rng = np.random.default_rng(SEED + seed_offset + i)
        rows.append(robust_row(gene, diff_wide.loc[gene].dropna().to_numpy(float), rng))
    out = pd.DataFrame(rows)
    for col in ["signflip_p", "wilcoxon_p", "sign_test_p"]:
        out[col.replace("_p", "_fdr_bh")] = bh(out[col])
    return out


def historical_diffs_from_h5ad(genes: list[str], donors: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rebuild donor-region log-normalized means for every locked candidate.

    The historical P2 means file covered only 103 preselected genes and must not
    be used as the universe for the locked 82-gene inference table.
    """
    a = ad.read_h5ad(H5AD)
    name_to_idx = {str(name): i for i, name in enumerate(a.var["feature_name"].astype(str))}
    missing = [g for g in genes if g not in name_to_idx]
    if missing:
        raise RuntimeError(f"Candidate genes absent from H5AD: {missing}")
    idx = [name_to_idx[g] for g in genes]
    rows = []
    donor_values = a.obs["donor_id"].astype(str).to_numpy()
    region_values = a.obs["Region"].astype(str).to_numpy()
    for donor in donors:
        for region in ("Macular", "Peripheral"):
            mask = (donor_values == donor) & (region_values == region)
            values = np.asarray(a.X[mask][:, idx].mean(axis=0)).ravel()
            rows.extend({"donor_id": donor, "Region": region, "gene": gene, "mean_expr": value}
                        for gene, value in zip(genes, values))
    means = pd.DataFrame(rows)
    paired = means.pivot(index="gene", columns=["donor_id", "Region"], values="mean_expr")
    diffs = pd.DataFrame(index=genes, columns=donors, dtype=float)
    for donor in donors:
        diffs[donor] = paired[(donor, "Macular")] - paired[(donor, "Peripheral")]
    return diffs, means


def age_numeric(value: str) -> float:
    text = str(value).strip().lower()
    if "day" in text:
        return float(text.split()[0].replace(">", "")) / 365.25
    return float(text.split()[0].replace(">", ""))


def sensitivity_rows(diff_wide: pd.DataFrame, donors: list[str], label: str, cells: pd.DataFrame) -> list[dict]:
    rows = []
    for gene in CORE5:
        d = diff_wide.loc[gene, donors].dropna().to_numpy(float)
        base = robust_row(gene, d, np.random.default_rng(SEED + len(rows) + 5000))
        rows.append({
            "analysis_group": label, "gene": gene, "n_donors": len(d),
            "mean_difference": base["mean_difference"], "ci_low": base["classical_ci_low"],
            "ci_high": base["classical_ci_high"], "direction": "Macula_high" if base["mean_difference"] > 0 else "Peripheral_high",
            "positive_donor_count": base["positive_donor_count"], "support_fraction": base["support_fraction"],
            "signflip_p": base["signflip_p"],
            "minimum_macular_cells": int(cells.loc[donors, "Macular"].min()) if donors else np.nan,
            "minimum_peripheral_cells": int(cells.loc[donors, "Peripheral"].min()) if donors else np.nan,
            "formal_inference": "YES" if len(d) >= 6 else "NO_SAMPLE_STRUCTURE_ONLY",
        })
    return rows


def main() -> None:
    old = pd.read_csv(DISCOVERY_RESULTS)
    edger = pd.read_csv(OUT / "RAWCOUNT_PSEUDOBULK_ALL_GENES.tsv", sep="\t")
    candidates82 = old.loc[(old.paired_t_fdr < 0.05) & (old.mean_mac_minus_per > 0), "gene"].tolist()
    if len(candidates82) != 82:
        raise RuntimeError(f"Historical candidate set is {len(candidates82)}, expected 82")

    logcpm = pd.read_csv(OUT / "RAWCOUNT_TMM_LOGCPM.tsv.gz", sep="\t")
    value_cols = [c for c in logcpm.columns if "__" in c]
    lc = logcpm.set_index("gene")[value_cols]
    donors = sorted({c.split("__")[0] for c in value_cols})
    diffs = pd.DataFrame(index=lc.index, columns=donors, dtype=float)
    for donor in donors:
        diffs[donor] = lc[f"{donor}__Macular"] - lc[f"{donor}__Peripheral"]
    diffs.to_csv(OUT / "RAWCOUNT_TMM_DONOR_PAIRED_DIFFERENCES.tsv.gz", sep="\t", compression="gzip")

    support = pd.DataFrame({
        "gene": diffs.index,
        "rawcount_mean_logcpm_difference": diffs.mean(axis=1).to_numpy(),
        "rawcount_positive_donor_count": (diffs > 0).sum(axis=1).to_numpy(),
        "rawcount_support_fraction": (diffs > 0).mean(axis=1).to_numpy(),
    })
    comp = old[["gene_id", "gene", "mean_mac_minus_per", "paired_t_fdr", "support_fraction"]].merge(
        edger, on=["gene_id", "gene"], how="outer"
    ).merge(support, on="gene", how="left")
    comp["old_candidate82"] = (comp.paired_t_fdr < 0.05) & (comp.mean_mac_minus_per > 0)
    comp["new_macula_high_fdr05"] = (comp.FDR < 0.05) & (comp.logFC > 0)
    comp["candidate_change"] = np.select(
        [comp.old_candidate82 & comp.new_macula_high_fdr05, comp.old_candidate82 & ~comp.new_macula_high_fdr05,
         ~comp.old_candidate82 & comp.new_macula_high_fdr05],
        ["retained", "lost", "newly_gained"], default="neither"
    )
    comp["direction_concordant"] = np.sign(comp.mean_mac_minus_per) == np.sign(comp.logFC)
    comp.to_csv(OUT / "RAWCOUNT_VS_EXISTING_P3_COMPARISON.tsv", sep="\t", index=False)

    old_detected = comp.loc[comp.old_candidate82 & comp.logFC.notna()].copy()
    summary = pd.DataFrame({
        "metric": [
            "historical_candidate_n", "historical_candidates_tested_by_edger", "historical_candidates_retained_fdr05",
            "historical_candidate_direction_concordance", "all_tested_gene_effect_spearman",
            "new_macula_high_fdr05_n", "newly_gained_n", "lost_n"
        ],
        "value": [
            82, len(old_detected), int((old_detected.FDR < 0.05).sum()), old_detected.direction_concordant.mean(),
            comp[["mean_mac_minus_per", "logFC"]].dropna().corr(method="spearman").iloc[0, 1],
            int(comp.new_macula_high_fdr05.sum()), int((comp.candidate_change == "newly_gained").sum()),
            int((comp.candidate_change == "lost").sum()),
        ],
    })
    summary.to_csv(OUT / "RAWCOUNT_VS_EXISTING_P3_SUMMARY.tsv", sep="\t", index=False)

    # Robust inference using the historical donor-region mean-expression differences.
    old_diffs, locked_means = historical_diffs_from_h5ad(candidates82, donors)
    locked_means.to_csv(OUT / "CANDIDATE82_DONOR_REGION_LOGNORMALIZED_MEANS.tsv.gz", sep="\t", index=False, compression="gzip")
    core = inference_table(old_diffs, CORE5)
    cand = inference_table(old_diffs, candidates82, 100)
    core.to_csv(OUT / "CORE5_ROBUST_INFERENCE.tsv", sep="\t", index=False)
    cand.to_csv(OUT / "CANDIDATE82_ROBUST_INFERENCE.tsv", sep="\t", index=False)

    new_genes = comp.loc[comp.new_macula_high_fdr05, "gene"].dropna().tolist()
    raw_robust = inference_table(diffs, new_genes, 1000)
    raw_robust.to_csv(OUT / "RAWCOUNT_ROBUST_CANDIDATE_INFERENCE.tsv", sep="\t", index=False)

    # Cell-count, age and source sensitivities use TMM logCPM paired differences.
    cell_long = pd.read_csv(RESULTS_ROOT / "00_DATA_AUDIT" / "DONOR_REGION_CELL_COUNTS.tsv", sep="\t")
    cell_wide = cell_long.pivot(index="donor_id", columns="Region", values="n_cells").fillna(0)
    min_rows = []
    for threshold in (10, 20, 50, 100):
        selected = sorted(cell_wide.index[(cell_wide.Macular >= threshold) & (cell_wide.Peripheral >= threshold)].intersection(donors))
        min_rows.extend(sensitivity_rows(diffs, selected, f">={threshold}_cells_per_region", cell_wide))
    pd.DataFrame(min_rows).to_csv(OUT / "MIN_CELL_THRESHOLD_SENSITIVITY.tsv", sep="\t", index=False)

    meta = pd.read_csv(OUT / "RAWCOUNT_PSEUDOBULK_SAMPLE_METADATA.tsv", sep="\t")
    donor_meta = meta.drop_duplicates("donor_id")[["donor_id", "donor_age", "Study", "sex"]].copy()
    donor_meta["age_years"] = donor_meta.donor_age.map(age_numeric)
    donor_meta["paired_analysis"] = True
    donor_meta.to_csv(OUT / "PAIRED_DONOR_AGE_DISTRIBUTION.tsv", sep="\t", index=False)
    age_rows = []
    for label, selected in [
        ("all_paired", donor_meta.donor_id.tolist()),
        ("adult_ge18", donor_meta.loc[donor_meta.age_years >= 18, "donor_id"].tolist()),
        ("age_ge50", donor_meta.loc[donor_meta.age_years >= 50, "donor_id"].tolist()),
    ]:
        if label != "all_paired" and len(selected) < 8:
            for gene in CORE5:
                age_rows.append({"analysis_group": label, "gene": gene, "n_donors": len(selected), "formal_inference": "NO_N_LT_8"})
        else:
            age_rows.extend(sensitivity_rows(diffs, sorted(selected), label, cell_wide))
    pd.DataFrame(age_rows).to_csv(OUT / "AGE_SENSITIVITY.tsv", sep="\t", index=False)

    source_rows = []
    for label, selected in [
        ("all_paired", donor_meta.donor_id.tolist()),
        ("Chen_lab_only", donor_meta.loc[donor_meta.Study == "Chen_lab", "donor_id"].tolist()),
    ]:
        source_rows.extend(sensitivity_rows(diffs, sorted(selected), label, cell_wide))
    pd.DataFrame(source_rows).to_csv(OUT / "SOURCE_RESTRICTED_SENSITIVITY.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
