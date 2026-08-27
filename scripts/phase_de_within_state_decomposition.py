#!/usr/bin/env python3
"""Within-state paired inference and exact regional effect decomposition."""

from itertools import product
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("M28_INPUT_DIR", PROJECT_ROOT / "inputs"))
RESULTS_ROOT = Path(os.environ.get("M28_RESULTS_DIR", PROJECT_ROOT / "results"))
DISCOVERY_RESULTS = Path(os.environ.get("M28_DISCOVERY_RESULTS", INPUT_ROOT / "discovery_paired_results.csv"))
OUT = RESULTS_ROOT / "02_RPE_STATE_DECOMPOSITION"
CORE5 = ["WFDC1", "PDE3B", "SULF1", "COL8A1", "ABCA1"]


def signflip(diff):
    d = np.asarray(diff, float)
    signs = np.array(list(product((-1.0, 1.0), repeat=len(d))))
    return float(np.mean(np.abs(signs @ d / len(d)) >= abs(d.mean()) - 1e-14))


def infer(d):
    d = np.asarray(d, float)
    n = len(d)
    sd = d.std(ddof=1) if n > 1 else np.nan
    se = sd / np.sqrt(n) if n > 1 else np.nan
    half = stats.t.ppf(.975, n - 1) * se if n > 1 else np.nan
    return {
        "n_paired_donors": n, "mean_paired_difference": d.mean() if n else np.nan,
        "ci_low": d.mean() - half if n else np.nan, "ci_high": d.mean() + half if n else np.nan,
        "positive_donor_count": int((d > 0).sum()), "support_fraction": float((d > 0).mean()) if n else np.nan,
        "signflip_p": signflip(d) if n else np.nan,
        "wilcoxon_p": float(stats.wilcoxon(d, method="auto").pvalue) if n and np.any(d != 0) else 1.0,
    }


def main():
    old = pd.read_csv(DISCOVERY_RESULTS)
    genes = old.loc[(old.paired_t_fdr < .05) & (old.mean_mac_minus_per > 0), "gene"].tolist()
    means = pd.read_csv(OUT / "STATE_CANDIDATE82_DONOR_REGION_MEANS.tsv.gz", sep="\t")
    meta = pd.read_csv(OUT / "STATE_PSEUDOBULK_SAMPLE_METADATA.tsv", sep="\t")
    edger = pd.read_csv(OUT / "STATE_RAWCOUNT_EDGER_ALL_GENES.tsv.gz", sep="\t")
    states = sorted(means.state.unique())
    rows = []
    for threshold in (20, 10):
        for state in states:
            sm = meta[meta.state == state].pivot(index="donor_id", columns="Region", values="n_cells").fillna(0)
            donors = sorted(sm.index[(sm.Macular >= threshold) & (sm.Peripheral >= threshold)])
            for gene in genes:
                row = {"state": state, "threshold_cells": threshold, "gene": gene,
                       "formal_testing_eligible": len(donors) >= 8}
                if len(donors) >= 8:
                    sub = means[(means.state == state) & means.donor_id.isin(donors)]
                    wide = sub.pivot(index="donor_id", columns="Region", values=gene).dropna()
                    row.update(infer((wide.Macular - wide.Peripheral).to_numpy()))
                else:
                    row.update({"n_paired_donors": len(donors), "mean_paired_difference": np.nan,
                                "ci_low": np.nan, "ci_high": np.nan, "positive_donor_count": np.nan,
                                "support_fraction": np.nan, "signflip_p": np.nan, "wilcoxon_p": np.nan})
                rows.append(row)
    results = pd.DataFrame(rows)
    results["signflip_fdr_bh"] = np.nan
    for (state, threshold), idx in results.groupby(["state", "threshold_cells"]).groups.items():
        valid = results.loc[idx, "signflip_p"].notna()
        valid_idx = results.loc[idx].index[valid]
        if len(valid_idx):
            results.loc[valid_idx, "signflip_fdr_bh"] = multipletests(results.loc[valid_idx, "signflip_p"], method="fdr_bh")[1]
    edger_small = edger[edger.gene.isin(genes)].rename(columns={
        "logFC": "rawcount_logFC", "logCPM": "rawcount_logCPM", "F": "rawcount_F",
        "PValue": "rawcount_PValue", "FDR": "rawcount_FDR"
    })
    results["paired_donor_n"] = results["n_paired_donors"]
    results = results.merge(edger_small.drop(columns="gene_id"), on=["state", "threshold_cells", "gene", "paired_donor_n"], how="left")
    results[results.gene.isin(CORE5)].to_csv(OUT / "CORE5_WITHIN_STATE_RESULTS.tsv", sep="\t", index=False)
    results.to_csv(OUT / "CANDIDATE82_WITHIN_STATE_RESULTS.tsv", sep="\t", index=False)

    # Exact decomposition using normalized .X means. Primary coarse scheme is fully observed;
    # strict five-state decomposition is retained only for complete donor pairs.
    state_counts = meta.pivot_table(index=["donor_id", "Region"], columns="state", values="n_cells", fill_value=0)
    mean_idx = means.set_index(["donor_id", "Region", "state"])
    donors = sorted(set(means.donor_id))
    decomp_rows = []
    for scheme in ("atlas5_complete_case", "LGI1_vs_other"):
        for donor in donors:
            for gene in genes:
                if scheme == "atlas5_complete_case":
                    scheme_states = states
                    if any(state_counts.loc[(donor, r), s] <= 0 for r in ("Macular", "Peripheral") for s in scheme_states):
                        continue
                    p = {r: state_counts.loc[(donor, r), scheme_states] / state_counts.loc[(donor, r), scheme_states].sum() for r in ("Macular", "Peripheral")}
                    mu = {r: pd.Series({s: mean_idx.loc[(donor, r, s), gene] for s in scheme_states}) for r in ("Macular", "Peripheral")}
                else:
                    scheme_states = ["RPE_LGI1+", "non_LGI1"]
                    p = {}
                    mu = {}
                    valid = True
                    for r in ("Macular", "Peripheral"):
                        c_lgi = state_counts.loc[(donor, r), "RPE_LGI1+"]
                        c_all = state_counts.loc[(donor, r)].sum()
                        c_other = c_all - c_lgi
                        if c_lgi <= 0 or c_other <= 0:
                            valid = False
                            break
                        lgi_mu = mean_idx.loc[(donor, r, "RPE_LGI1+"), gene]
                        other_states = [s for s in states if s != "RPE_LGI1+" and state_counts.loc[(donor, r), s] > 0]
                        other_mu = sum(state_counts.loc[(donor, r), s] * mean_idx.loc[(donor, r, s), gene] for s in other_states) / c_other
                        p[r] = pd.Series({"RPE_LGI1+": c_lgi / c_all, "non_LGI1": c_other / c_all})
                        mu[r] = pd.Series({"RPE_LGI1+": lgi_mu, "non_LGI1": other_mu})
                    if not valid:
                        continue
                within = float((((p["Macular"] + p["Peripheral"]) / 2) * (mu["Macular"] - mu["Peripheral"])).sum())
                composition = float(((p["Macular"] - p["Peripheral"]) * ((mu["Macular"] + mu["Peripheral"]) / 2)).sum())
                observed_region = {}
                for r in ("Macular", "Peripheral"):
                    region_total = state_counts.loc[(donor, r)].sum()
                    observed_region[r] = sum(
                        state_counts.loc[(donor, r), s] * mean_idx.loc[(donor, r, s), gene]
                        for s in states if state_counts.loc[(donor, r), s] > 0
                    ) / region_total
                observed = float(observed_region["Macular"] - observed_region["Peripheral"])
                reconstructed = within + composition
                decomp_rows.append({
                    "decomposition_scheme": scheme, "donor": donor, "gene": gene,
                    "observed_total_difference": observed, "within_state_component": within,
                    "composition_component": composition, "reconstructed_difference": reconstructed,
                    "reconstruction_error": observed - reconstructed,
                    "within_state_fraction": within / reconstructed if abs(reconstructed) > 1e-12 else np.nan,
                })
    decomp = pd.DataFrame(decomp_rows)
    decomp[decomp.gene.isin(CORE5)].to_csv(OUT / "REGIONAL_EFFECT_DECOMPOSITION_CORE5.tsv", sep="\t", index=False)
    decomp.to_csv(OUT / "REGIONAL_EFFECT_DECOMPOSITION_CANDIDATE82.tsv", sep="\t", index=False)

    summary = decomp.groupby(["decomposition_scheme", "gene"]).agg(
        donor_n=("donor", "nunique"), mean_total=("observed_total_difference", "mean"),
        mean_within=("within_state_component", "mean"), mean_composition=("composition_component", "mean"),
        max_abs_reconstruction_error=("reconstruction_error", lambda x: np.abs(x).max())
    ).reset_index()
    summary["within_share_of_absolute_components"] = np.abs(summary.mean_within) / (np.abs(summary.mean_within) + np.abs(summary.mean_composition))
    summary["driver_class"] = np.select(
        [summary.within_share_of_absolute_components >= .67, summary.within_share_of_absolute_components <= .33],
        ["within_state_dominant", "composition_dominant"], default="mixed"
    )
    summary.to_csv(OUT / "REGIONAL_EFFECT_DECOMPOSITION_SUMMARY.tsv", sep="\t", index=False)
    core_summary = summary[(summary.gene.isin(CORE5)) & (summary.decomposition_scheme == "LGI1_vs_other")]
    strict_n = int(summary.loc[summary.decomposition_scheme == "atlas5_complete_case", "donor_n"].max()) if (summary.decomposition_scheme == "atlas5_complete_case").any() else 0
    lines = [f"- {r.gene}: {r.driver_class}; mean within-state={r.mean_within:.4f}, mean composition={r.mean_composition:.4f}, absolute-component within share={r.within_share_of_absolute_components:.1%}." for r in core_summary.itertuples()]
    report = f"""# Interpretation of regional effect decomposition

## Decomposition feasibility

The five original atlas states are not jointly observed in both regions for most donors. A strict five-state decomposition is therefore available for only {strict_n} donor pairs and is descriptive only. Missing conditional state means were not set to zero and were not imputed.

To retain an exact, auditable donor-level decomposition across a substantially larger paired set, the primary sensitivity collapses the same atlas labels into `RPE_LGI1+` versus `non_LGI1`. This two-state partition covers all RPE cells and preserves exact reconstruction wherever both partitions are observed in both regions.

## Core-gene results under the LGI1-versus-other partition

{chr(10).join(lines)}

`within_state_dominant`, `composition_dominant`, and `mixed` are descriptive classifications based on the mean absolute component shares (≥67%, ≤33%, or intermediate). Opposing component signs can produce fractions outside 0–1 at individual-donor level; the classification therefore uses absolute aggregated components rather than raw individual fractions.

## Scientific boundary

Within-LGI1+ paired testing can establish whether a regional expression difference persists within the only state that satisfies the primary ≥20-cell/≥8-donor gate. The coarse decomposition can separate a dominant-state composition contribution from within-partition expression, but it cannot identify a unique causal RPE transition or substitute for a well-powered five-state design.
"""
    (OUT / "STATE_DECOMPOSITION_INTERPRETATION.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
