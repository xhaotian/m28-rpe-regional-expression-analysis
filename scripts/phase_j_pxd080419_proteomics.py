#!/usr/bin/env python3
"""Apply the PXD080419 processed-table feasibility gate."""
from pathlib import Path
import json
import os
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("M28_INPUT_DIR", PROJECT_ROOT / "inputs"))
RESULTS_ROOT = Path(os.environ.get("M28_RESULTS_DIR", PROJECT_ROOT / "results"))
OUT = RESULTS_ROOT / "07_PXD080419_PROTEOMICS"
API = Path(os.environ.get("M28_PXD080419_METADATA", INPUT_ROOT / "PXD080419" / "metadata.json"))
CAND = Path(os.environ.get("M28_DISCOVERY_RESULTS", INPUT_ROOT / "discovery_paired_results.csv"))
CORE5 = {"WFDC1", "PDE3B", "SULF1", "COL8A1", "ABCA1"}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    j = json.loads(API.read_text())
    urls = [x["value"] for x in j["datasetFiles"]]
    rows = [
        {"donor": "NOT_EXPOSED_IN_PUBLIC_INDEX", "region": "Macula_and_Peripheral_paired_design_stated", "fraction": "EV", "replicate": "UNKNOWN", "protein_ID": "NOT_AVAILABLE", "quantification_type": "NOT_AVAILABLE", "public_archive": next(u for u in urls if "EV Proteome Raw" in u), "archive_bytes": 6493684314, "file_role": "raw_data"},
        {"donor": "NOT_EXPOSED_IN_PUBLIC_INDEX", "region": "Macula_and_Peripheral_paired_design_stated", "fraction": "Soluble_secretome", "replicate": "UNKNOWN", "protein_ID": "NOT_AVAILABLE", "quantification_type": "NOT_AVAILABLE", "public_archive": next(u for u in urls if "EV-depleted" in u), "archive_bytes": 14214543266, "file_role": "raw_data"},
        {"donor": "NOT_EXPOSED_IN_PUBLIC_INDEX", "region": "UNKNOWN", "fraction": "Search_project", "replicate": "UNKNOWN", "protein_ID": "NOT_AVAILABLE", "quantification_type": "Spectronaut_sne_proprietary_project", "public_archive": next(u for u in urls if "Spectronaut" in u), "archive_bytes": "UNKNOWN", "file_role": "search_engine_output_not_processed_table"},
    ]
    pd.DataFrame(rows).to_csv(OUT / "PXD080419_SAMPLE_MAP.tsv", sep="\t", index=False)

    locked = pd.read_csv(CAND)
    genes = locked.loc[(locked.paired_t_fdr < 0.05) & (locked.mean_mac_minus_per > 0), "gene"].astype(str).tolist()
    if len(genes) != 82:
        raise RuntimeError(f"Locked candidate set has {len(genes)} genes, expected 82")
    det = []
    for fraction in ["EV", "Soluble_secretome"]:
        for gene in genes:
            det.append({"gene": gene, "core5": gene in CORE5, "fraction": fraction, "detectable": "NOT_ASSESSABLE", "reason": "no_public_processed_quantitative_protein_table_or_sample_map"})
    pd.DataFrame(det).to_csv(OUT / "PXD080419_CANDIDATE_DETECTABILITY.tsv", sep="\t", index=False)

    concord = []
    enrich = []
    for fraction in ["EV", "Soluble_secretome"]:
        concord.append({"fraction": fraction, "candidate_n": 82, "quantifiable_candidate_n": "NA", "macula_high_n": "NA", "complete_pair_n": "NA", "status": "NOT_ASSESSABLE", "reason": "processed abundance table and donor/sample map unavailable"})
        enrich.append({"fraction": fraction, "method": "Fisher_exact_and_matched_random_set", "effect": "NA", "p_value": "NA", "status": "NOT_ASSESSABLE", "reason": "protein-level background universe and regional quantitative results unavailable"})
    pd.DataFrame(concord).to_csv(OUT / "PXD080419_TRANSCRIPT_PROTEIN_CONCORDANCE.tsv", sep="\t", index=False)
    pd.DataFrame(enrich).to_csv(OUT / "PXD080419_PROTEIN_ENRICHMENT.tsv", sep="\t", index=False)

    report = """# PXD080419 proteomics feasibility decision

## Decision

**FAIL for quantitative protein concordance with the currently public files.** ProteomeCentral/iProX confirms the paired human macular/peripheral RPE-choroid explant EV and EV-depleted soluble-secretome study, but the public index exposes two large raw-data archives and a proprietary Spectronaut `.sne` search-project archive. It does not expose a processed protein-abundance table, protein identifiers, or a donor/region/replicate sample map.

The two raw archives are approximately 6.49 GB and 14.21 GB. Reprocessing roughly 20.7 GB of raw mass-spectrometry data, reverse-engineering the design, and reproducing the proprietary search workflow is outside the preregistered “prioritize processed quantitative tables” branch and would introduce a new analysis route. Those archives were therefore not downloaded.

## Interpretation

- Candidate detectability, paired protein direction, and set-level enrichment are **not assessable**, not negative.
- No protein result is counted in the integrated evidence tiers.
- This branch can be reopened if the authors or repository release a processed quantitative export plus sample annotation.
"""
    (OUT / "PXD080419_FEASIBILITY_DECISION.md").write_text(report)


if __name__ == "__main__":
    main()
