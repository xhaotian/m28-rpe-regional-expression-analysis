#!/usr/bin/env python3
"""Cross-check M28 genes against the published GSE135092 RPE regional DE table."""
from __future__ import annotations

import gzip
import io
import os
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("M28_INPUT_DIR", PROJECT_ROOT / "inputs"))
RESULTS_ROOT = Path(os.environ.get("M28_RESULTS_DIR", PROJECT_ROOT / "results"))
OUT = RESULTS_ROOT / "06_GSE135092_BULK"
DISCOVERY_RESULTS = Path(os.environ.get("M28_DISCOVERY_RESULTS", INPUT_ROOT / "discovery_paired_results.csv"))
META = Path(os.environ.get("M28_GSE135092_METADATA", INPUT_ROOT / "GSE135092" / "metadata.tsv"))
TAR = Path(os.environ.get("M28_GSE135092_RAW", INPUT_ROOT / "GSE135092" / "GSE135092_RAW.tar"))
DATA_S1 = Path(os.environ.get("M28_GSE135092_SUPPLEMENT", INPUT_ROOT / "GSE135092" / "mmc2.xlsx"))
CORE5 = ["WFDC1", "PDE3B", "SULF1", "COL8A1", "ABCA1"]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    locked = pd.read_csv(DISCOVERY_RESULTS)
    genes = locked.loc[(locked.paired_t_fdr < 0.05) & (locked.mean_mac_minus_per > 0), "gene"].astype(str).tolist()
    if len(genes) != 82:
        raise RuntimeError(f"Locked candidate set has {len(genes)} genes, expected 82")
    p3 = pd.read_csv(DISCOVERY_RESULTS)
    gene_to_id = dict(zip(p3.gene.astype(str), p3.gene_id.astype(str)))
    targets = {gene_to_id[g]: g for g in genes if g in gene_to_id}

    meta = pd.read_csv(META, sep="\t")
    outmeta = meta[["sample_id", "samid", "tissue", "region", "disease_status", "age", "supplementary_file_1"]].copy()
    outmeta["donor_identifier"] = "NOT_PROVIDED_IN_GEO_METADATA"
    outmeta["paired_status"] = "UNRESOLVED"
    outmeta["pairing_evidence"] = "SAM IDs are sample IDs; age/order/BioSample adjacency were not used to infer donors"
    outmeta.to_csv(OUT / "GSE135092_SAMPLE_METADATA.tsv", sep="\t", index=False)

    eligible = meta[(meta.tissue.str.upper() == "RPE") & (meta.disease_status.str.lower() == "control")]
    gsm_set = set(eligible.sample_id)
    observations = {g: [] for g in genes}
    with tarfile.open(TAR, "r") as archive:
        for member in archive.getmembers():
            gsm = member.name.split("_")[0]
            if gsm not in gsm_set or not member.isfile():
                continue
            raw = archive.extractfile(member)
            if raw is None:
                continue
            with gzip.GzipFile(fileobj=raw) as gz, io.TextIOWrapper(gz, encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("#") or line.startswith("ID_REF"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if parts[0] in targets:
                        observations[targets[parts[0]]].append((float(parts[1]), float(parts[2])))

    published = pd.read_excel(DATA_S1, sheet_name="DE_rpe_non_macula_vs_macula", header=1)
    published = published.set_index("GeneSymbol", drop=False)
    rows = []
    for gene in genes:
        obs = observations.get(gene, [])
        if gene in published.index:
            r = published.loc[gene]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            effect = float(r.Log2FoldChange)
            fdr = float(r.AdjustedPvalue)
            listed = True
            direction = "Non_macula_high" if effect > 0 else "Macula_high"
        else:
            effect, fdr, listed, direction = np.nan, np.nan, False, "Not_listed"
        rows.append({
            "gene": gene,
            "gene_id": gene_to_id.get(gene, ""),
            "detectable": bool(obs and any(count > 0 for _, count in obs)),
            "control_rpe_samples_observed": len(obs),
            "control_rpe_samples_count_positive": sum(count > 0 for _, count in obs),
            "mean_control_rpe_nRPKM": np.mean([v for v, _ in obs]) if obs else np.nan,
            "published_region_DE_listed": listed,
            "direction": direction,
            "effect_statistic": effect,
            "effect_definition": "published log2 fold-change: RPE non-macula versus macula",
            "published_FDR": fdr,
            "published_or_regenerated": "published_Data_S1" if listed else "not_listed_in_published_Data_S1_DE_table",
            "support_interpretation": (
                "concordant_tissue_level" if direction == "Macula_high" else
                "discordant_tissue_level" if direction == "Non_macula_high" else
                "detectable_but_no_published_region_DE_statistic"
            ),
        })
    res = pd.DataFrame(rows)
    res[res.gene.isin(CORE5)].to_csv(OUT / "GSE135092_CORE5_CROSSCHECK.tsv", sep="\t", index=False)
    res.to_csv(OUT / "GSE135092_CANDIDATE82_CROSSCHECK.tsv", sep="\t", index=False)

    core = res[res.gene.isin(CORE5)]
    report = f"""# GSE135092 cross-platform decision

Donor identifiers are not present in official GEO sample metadata. Pairing was therefore not reconstructed from age, SAM number, BioSample adjacency, or sample order. The preregistered fallback was used: the Orozco et al. Data S1 sheet `DE_rpe_non_macula_vs_macula`.

The table's fold-change is interpreted according to its explicit contrast label: positive values are non-macula-high and negative values are macula-high. Of the core five, {int(core.published_region_DE_listed.sum())} appears in the published significant regional DE table: {', '.join(core.loc[core.published_region_DE_listed, 'gene']) or 'none'}. This is **RPE/choroid tissue-level cross-platform evidence**, not pure-RPE replication.

Absence from Data S1 means that no published regional statistic was available in that DE table; it is not treated as absence of expression. Detectability was independently checked in the GEO per-sample processed TSV files.
"""
    (OUT / "GSE135092_FEASIBILITY_DECISION.md").write_text(report)


if __name__ == "__main__":
    main()
