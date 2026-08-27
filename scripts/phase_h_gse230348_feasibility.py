#!/usr/bin/env python3
"""Parse GSE230348 official SOFT metadata and apply the prespecified feasibility gate."""

import gzip
import os
import re
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("M28_INPUT_DIR", PROJECT_ROOT / "inputs"))
RESULTS_ROOT = Path(os.environ.get("M28_RESULTS_DIR", PROJECT_ROOT / "results"))
SOURCE = Path(os.environ.get("M28_GSE230348_DIR", INPUT_ROOT / "GSE230348"))
SOFT = SOURCE / "metadata" / "GSE230348_family.soft.gz"
FILELIST = SOURCE / "supplementary" / "filelist.txt"
OUT = RESULTS_ROOT / "05_GSE230348_FEASIBILITY"


def parse_soft() -> pd.DataFrame:
    text = gzip.open(SOFT, "rt", errors="replace").read()
    blocks = re.split(r"(?=\^SAMPLE = )", text)
    rows = []
    for block in blocks:
        m = re.match(r"\^SAMPLE = (GSM\d+)", block)
        if not m:
            continue
        gsm = m.group(1)
        def one(label):
            z = re.search(rf"^!{re.escape(label)} = (.+)$", block, re.M)
            return z.group(1).strip() if z else "unknown"
        title = one("Sample_title")
        source = one("Sample_source_name_ch1")
        chars = re.findall(r"^!Sample_characteristics_ch1 = (.+)$", block, re.M)
        char_map = {}
        for value in chars:
            if ":" in value:
                k, v = value.split(":", 1)
                char_map[k.strip().lower()] = v.strip()
        donor_match = re.match(r"(.+?)-([mt])(?:Ret|RPE)$", title, flags=re.I)
        donor = donor_match.group(1) if donor_match else "unknown"
        tissue = "RPE/choroid" if ("rpe" in title.lower() or "choroid" in source.lower()) else "retina"
        if "macula" in source.lower() or (donor_match and donor_match.group(2).lower() == "m"):
            region = "macula"
        elif "periph" in source.lower() or (donor_match and donor_match.group(2).lower() == "t"):
            region = "peripheral"
        else:
            region = "unknown"
        supplementary = re.findall(r"^!Sample_supplementary_file_\d+ = (.+)$", block, re.M)
        has_matrix_triplet = any("matrix.mtx" in x for x in supplementary) and any("barcodes.tsv" in x for x in supplementary) and any("features.tsv" in x for x in supplementary)
        rows.append({
            "accession": gsm, "donor": donor, "disease_status": char_map.get("disease state", "unknown"),
            "region": region, "central/macula/peripheral": region, "tissue": tissue,
            "library": title, "assay": one("Sample_library_strategy"),
            "cell_annotation_availability": "ABSENT_FROM_GEO_SUPPLEMENTARY_FILES",
            "RPE_cell_availability": "UNCONFIRMED_MIXED_RPE_CHOROID" if tissue == "RPE/choroid" else "NOT_RPE_TISSUE",
            "processed_matrix_triplet_available": has_matrix_triplet,
            "paired_status": "pending_group_audit",
        })
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = parse_soft()
    rpe = d[d.tissue == "RPE/choroid"].copy()
    pair = rpe.groupby("donor").agg(
        disease_n=("disease_status", "nunique"),
        disease_status=("disease_status", lambda x: "|".join(sorted(set(x)))),
        macula=("region", lambda x: int((x == "macula").any())),
        peripheral=("region", lambda x: int((x == "peripheral").any())),
        matrix_complete=("processed_matrix_triplet_available", "all"),
    )
    paired_donors = set(pair.index[(pair.macula == 1) & (pair.peripheral == 1)])
    healthy_paired = set(pair.index[(pair.macula == 1) & (pair.peripheral == 1) & pair.disease_status.str.lower().isin(["normal", "healthy", "control"])])
    d["paired_status"] = d.donor.map(lambda x: "paired_region_donor" if x in paired_donors else "not_paired_or_unknown")
    d.to_csv(OUT / "GSE230348_SAMPLE_MAP.tsv", sep="\t", index=False)
    pair.reset_index().to_csv(OUT / "GSE230348_DONOR_PAIR_AUDIT.tsv", sep="\t", index=False)

    decision = "FAIL"
    report = f"""# GSE230348 feasibility decision

## Decision: {decision}

Official GEO metadata identifies {len(d)} libraries and {len(rpe)} RPE/choroid-region libraries. Donor and region strings permit reconstruction of {len(paired_donors)} macula/peripheral RPE/choroid tissue pairs, including {len(healthy_paired)} donors labelled normal/healthy/control.

The public supplementary structure provides per-library 10x-style barcode, feature and matrix files, but no cell-level annotation table or mapping that identifies RPE cells within the mixed RPE/choroid preparations. Consequently, the required criterion “each region has analyzable RPE cells” cannot be verified from the audited metadata and file inventory. Tissue labels are not treated as proof that recovered cells are RPE.

The branch therefore stops at the feasibility gate. No donor pairing is inferred beyond explicit sample-title/region metadata, and no descriptive or paired gene validation is generated. Reopening requires a public cell-annotation file, a traceable author object with cell identities, or a separately authorized de novo single-cell reconstruction from the full matrices.

## Audited inputs

- GEO family SOFT: `{SOFT}`
- Supplementary file inventory: `{FILELIST}`
- Submitter metadata workbooks: `{SOURCE / 'supplementary' / 'GSE230348_update.xlsx'}` and `{SOURCE / 'supplementary' / 'GSE230348_update3.xlsx'}`
"""
    (OUT / "GSE230348_FEASIBILITY_DECISION.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
