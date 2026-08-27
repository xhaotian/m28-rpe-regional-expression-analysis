#!/usr/bin/env python3
"""Build the M28 lineage audit and locked 82-gene evidence matrix."""
from pathlib import Path
import os
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("M28_INPUT_DIR", PROJECT_ROOT / "inputs"))
ROOT = Path(os.environ.get("M28_RESULTS_DIR", PROJECT_ROOT / "results"))
OUT = ROOT / "08_INTEGRATED_EVIDENCE"
GENE_ACTIVITY_RESULTS = Path(os.environ.get(
    "M28_GSE220155_GENE_ACTIVITY",
    INPUT_ROOT / "GSE220155" / "gene_activity_support.tsv",
))
CORE5 = {"WFDC1", "PDE3B", "SULF1", "COL8A1", "ABCA1"}


def build_lineage():
    resources = {
        "current_CELLxGENE_RPE_dataset": {"study": "2026_integrated_RPE_choroid_atlas", "accession": "CELLxGENE_collection_48c15b0c_dataset_f28376c8", "donors": "Chen_lab_plus_Sanes_GSE236566", "assay": "snRNA-seq"},
        "Chen_lab_source": {"study": "2026_integrated_RPE_choroid_atlas_primary_source", "accession": "within_CELLxGENE_48c15b0c", "donors": "Chen_lab", "assay": "snRNA-seq"},
        "Sanes_GSE236566": {"study": "Sanes_et_al", "accession": "GSE236566", "donors": "Sanes_GSE236566", "assay": "snRNA-seq"},
        "GSE220155": {"study": "Mullin_et_al_2023", "accession": "GSE220155", "donors": "four_paired_donors_RNA_and_gene_activity", "assay": "single-cell_multiome_RNA_ATAC"},
        "GSE135922": {"study": "Voigt_et_al_2019", "accession": "GSE135922", "donors": "seven_paired_donors", "assay": "scRNA-seq"},
        "GSE230348": {"study": "GSE230348_depositor_study", "accession": "GSE230348", "donors": "nineteen_donors_six_healthy_region_pairs", "assay": "scRNA-seq"},
        "GSE135092": {"study": "Orozco_et_al_2020", "accession": "GSE135092", "donors": "129_donors_bulk_cohort", "assay": "bulk_RNA-seq"},
        "PXD080419": {"study": "regional_RPE_choroid_EV_secretome_study", "accession": "PXD080419_IPX0018068000", "donors": "four_nondiseased_paired_donors", "assay": "mass_spectrometry_proteomics"},
        "2026_integrated_RPE_choroid_atlas": {"study": "2026_integrated_RPE_choroid_atlas", "accession": "doi_10.64898_2026.07.23.740352_CELLxGENE_48c15b0c", "donors": "integrated_Chen_and_published_sources", "assay": "sc_snRNA-seq_and_snATAC-seq"},
    }
    names = list(resources)
    same_groups = [
        {"current_CELLxGENE_RPE_dataset", "Chen_lab_source", "Sanes_GSE236566", "2026_integrated_RPE_choroid_atlas"},
    ]
    rows = []
    for a in names:
        for b in names:
            if a == b:
                status, overlap, independent, same_study, same_acc = "SAME_DONORS_OR_LAYER", "YES", "NO", "YES", "YES"
                evidence = "self-comparison"
            elif any(a in g and b in g for g in same_groups):
                status, overlap, independent, same_study = "SAME_DONORS_OR_LAYER", "YES_OR_SUBSET", "NO", "YES_OR_SOURCE_WITHIN_ATLAS"
                same_acc = "YES_OR_NESTED_SOURCE"
                evidence = "The current H5AD is the RPE dataset from CELLxGENE collection 48c15b0c for the 2026 integrated atlas; obs Study contains Chen_lab and Sanes_GSE236566."
            else:
                status, overlap, independent, same_study, same_acc = "INDEPENDENT_DONORS", "NO_KNOWN_OVERLAP", "YES", "NO", "NO"
                evidence = "Distinct accession/study and no donor overlap identified in official metadata; independence applies biologically, subject to each branch's feasibility gate."
            rows.append({
                "dataset_row": a, "dataset_column": b,
                "row_accession": resources[a]["accession"], "column_accession": resources[b]["accession"],
                "same_study": same_study, "same_accession": same_acc,
                "known_donor_overlap": overlap, "possible_overlap": "NO_IDENTIFIED" if status == "INDEPENDENT_DONORS" else "YES",
                "independent_donor_set": independent,
                "independent_assay_only": "YES" if a == b and a == "GSE220155" else "NO",
                "independence_status": status, "evidence": evidence,
            })
    pd.DataFrame(rows).to_csv(OUT / "DATA_LINEAGE_MATRIX.tsv", sep="\t", index=False)
    pd.DataFrame([{"resource": k, **v} for k, v in resources.items()]).to_csv(OUT / "DATA_LINEAGE_RESOURCE_SUMMARY.tsv", sep="\t", index=False)


def index(df, gene="gene"):
    return df.drop_duplicates(gene).set_index(gene)


def build_evidence():
    robust = index(pd.read_csv(ROOT / "01_DONOR_ROBUSTNESS/CANDIDATE82_ROBUST_INFERENCE.tsv", sep="\t"))
    raw = index(pd.read_csv(ROOT / "01_DONOR_ROBUSTNESS/RAWCOUNT_VS_EXISTING_P3_COMPARISON.tsv", sep="\t").query("old_candidate82 == True"))
    locked_means = pd.read_csv(ROOT / "01_DONOR_ROBUSTNESS/CANDIDATE82_DONOR_REGION_LOGNORMALIZED_MEANS.tsv.gz", sep="\t")
    wide = locked_means.pivot(index="gene", columns=["donor_id", "Region"], values="mean_expr")
    donor_names = sorted(set(wide.columns.get_level_values(0)))
    diff = pd.DataFrame({d: wide[(d, "Macular")] - wide[(d, "Peripheral")] for d in donor_names})
    lodo_rows = []
    for gene, vals in diff.iterrows():
        leaveout = [vals.drop(d).mean() for d in donor_names]
        lodo_rows.append({"gene": gene, "lodo_min_mean_difference": min(leaveout), "lodo_max_mean_difference": max(leaveout), "lodo_stable_macula_up": all(x > 0 for x in leaveout), "lodo_n_direction_flips": sum(x <= 0 for x in leaveout)})
    lodo_df = pd.DataFrame(lodo_rows)
    lodo_df.to_csv(ROOT / "01_DONOR_ROBUSTNESS/CANDIDATE82_LODO_SENSITIVITY.tsv", sep="\t", index=False)
    lodo = index(lodo_df)
    within_all = pd.read_csv(ROOT / "02_RPE_STATE_DECOMPOSITION/CANDIDATE82_WITHIN_STATE_RESULTS.tsv", sep="\t")
    within = index(within_all[(within_all.state == "RPE_LGI1+") & (within_all.threshold_cells == 20)])
    decomp = pd.read_csv(ROOT / "02_RPE_STATE_DECOMPOSITION/REGIONAL_EFFECT_DECOMPOSITION_SUMMARY.tsv", sep="\t")
    decomp = index(decomp[decomp.decomposition_scheme == "LGI1_vs_other"])
    g220 = index(pd.read_csv(ROOT / "03_GSE220155_EMPIRICAL_BACKGROUND/GSE220155_CANDIDATE82_DETECTED_RESULTS.tsv", sep="\t"))
    ga = index(pd.read_csv(GENE_ACTIVITY_RESULTS, sep="\t"))
    g135 = index(pd.read_csv(ROOT / "04_GSE135922_EXTERNAL/GSE135922_CANDIDATE82.tsv", sep="\t").query("cell_threshold == 10"))
    bulk = index(pd.read_csv(ROOT / "06_GSE135092_BULK/GSE135092_CANDIDATE82_CROSSCHECK.tsv", sep="\t"))
    minc = pd.read_csv(ROOT / "01_DONOR_ROBUSTNESS/MIN_CELL_THRESHOLD_SENSITIVITY.tsv", sep="\t")
    age = pd.read_csv(ROOT / "01_DONOR_ROBUSTNESS/AGE_SENSITIVITY.tsv", sep="\t")
    source = pd.read_csv(ROOT / "01_DONOR_ROBUSTNESS/SOURCE_RESTRICTED_SENSITIVITY.tsv", sep="\t")

    rows = []
    for gene in robust.index:
        r, e, lo, wi, dc = robust.loc[gene], raw.loc[gene], lodo.loc[gene], within.loc[gene], decomp.loc[gene]
        x220 = g220.loc[gene] if gene in g220.index else None
        xga = ga.loc[gene] if gene in ga.index else None
        x135 = g135.loc[gene] if gene in g135.index else None
        xb = bulk.loc[gene] if gene in bulk.index else None
        if gene in CORE5:
            mc = minc[minc.gene == gene]
            ag = age[age.gene == gene]
            sr = source[(source.gene == gene) & (source.analysis_group == "Chen_lab_only")]
            mincell = "all_thresholds_macula_high" if (mc.direction == "Macula_high").all() else "direction_instability"
            adult = ";".join(f"{z.analysis_group}:{z.direction}:n={z.n_donors}" for _, z in ag.iterrows())
            chen = f"{sr.iloc[0].direction};support={sr.iloc[0].support_fraction:.3g};p={sr.iloc[0].signflip_p:.3g}"
        else:
            mincell = adult = chen = "NOT_TESTED_CORE5_ONLY"

        g220_joint = bool(x220.joint_positive_and_ge3) if x220 is not None else False
        g135_limited = bool(x135 is not None and x135.mean_macula_minus_peripheral > 0 and x135.positive_donor_n >= 3)
        bulk_support = bool(xb is not None and xb.direction == "Macula_high")
        # GSE135922 does not qualify as a full replication source because only
        # n=2 passes the primary cell gate and n=4 the sensitivity gate.
        if g220_joint:
            tier = "Tier B"
            rationale = "stable discovery plus independent GSE220155 RNA support"
        else:
            tier = "Tier C"
            rationale = "stable discovery; qualifying independent replication insufficient"
        rows.append({
            "gene": gene,
            "current_donor_aware_effect_logFC": e.logFC,
            "current_FDR": e.FDR,
            "donor_support": e.rawcount_support_fraction,
            "robust_permutation_p": r.signflip_p,
            "robust_permutation_FDR": r.signflip_fdr_bh,
            "LODO": bool(lo.lodo_stable_macula_up),
            "min_cell_sensitivity": mincell,
            "adult_only_sensitivity": adult,
            "Chen_only_sensitivity": chen,
            "within_state_persistence": f"LGI1_mean={wi.mean_paired_difference:.6g};support={wi.support_fraction:.3g};permFDR={wi.signflip_fdr_bh:.3g};rawFDR={wi.rawcount_FDR:.3g}",
            "composition_component": f"within_share={dc.within_share_of_absolute_components:.3g};class={dc.driver_class}",
            "GSE220155_RNA": "NOT_DETECTED" if x220 is None else f"delta={x220.mean_macula_minus_peripheral:.6g};positive={int(x220.positive_donor_number)}/4;joint={bool(x220.joint_positive_and_ge3)}",
            "GSE220155_empirical_null_support": "candidate_set_joint_empirical_P=9.999e-05",
            "GSE220155_gene_activity": "NOT_DETECTED" if xga is None or not bool(xga.is_detected) else f"delta={xga.ga_mean_delta:.6g};positive={int(xga.ga_n_positive)}/4",
            "GSE135922": "NOT_DETECTED" if x135 is None else f"LIMITED_ge10_delta={x135.mean_macula_minus_peripheral:.6g};positive={int(x135.positive_donor_n)}/{int(x135.paired_n)};qualifying={g135_limited}",
            "GSE230348": "FAIL_NO_PUBLIC_RPE_CELL_ANNOTATION",
            "GSE135092_bulk": "NOT_AVAILABLE_IN_PUBLISHED_DE_TABLE" if xb is None or not bool(xb.published_region_DE_listed) else f"{xb.direction};log2FC_nonmac_vs_mac={xb.effect_statistic:.6g};FDR={xb.published_FDR:.3g};support={bulk_support}",
            "PXD_EV": "NOT_ASSESSABLE_NO_PROCESSED_TABLE",
            "PXD_soluble_secretome": "NOT_ASSESSABLE_NO_PROCESSED_TABLE",
            "independence_status": "GSE220155/GSE135922/GSE135092 independent donors; GSE220155 RNA and gene activity same 4 donors; current H5AD is part of 2026 atlas",
            "evidence_tier": tier,
            "tier_rationale": rationale,
            "limited_additional_support_GSE135922": g135_limited,
            "tissue_level_support_GSE135092": bulk_support,
        })
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "INTEGRATED_GENE_EVIDENCE_MATRIX.tsv", sep="\t", index=False)
    summary = result.evidence_tier.value_counts().rename_axis("evidence_tier").reset_index(name="gene_n")
    summary.to_csv(OUT / "EVIDENCE_TIER_SUMMARY.tsv", sep="\t", index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    build_lineage()
    build_evidence()


if __name__ == "__main__":
    main()
