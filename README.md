# Human RPE Regional Expression Analysis

This repository contains the analysis code used to study gene-expression differences between macular and peripheral human retinal pigment epithelium (RPE). The workflow combines donor-paired discovery analysis, RPE-state analysis, robustness checks, and bounded comparisons with independent public datasets.

## Analysis overview

The discovery workflow audits a CELLxGENE H5AD object, aggregates raw UMI counts by donor and region, fits donor-adjusted edgeR quasi-likelihood models, and evaluates paired effect robustness. A second branch tests whether regional effects persist within annotated RPE states and decomposes the observed regional difference into within-state and composition-associated components.

Independent analyses use the following public resources: GSE220155 (GEO; paired RNA and same-donor gene-activity support), GSE135922 (GEO; underpowered directional check), GSE230348 (GEO; metadata/manual feasibility review only), GSE135092 (GEO; mixed RPE/choroid tissue-level cross-reference), and PXD080419 (ProteomeXchange/PRIDE; proteomics feasibility assessment). These datasets differ in assay, sample size, annotation depth, and tissue composition. Their results are therefore retained as dataset-specific support or feasibility assessments rather than pooled as interchangeable replications.

The canonical discovery resource is the CELLxGENE RPE dataset (collection `48c15b0c-6039-4e0f-9668-b6b3c0b830ad`; dataset `d21aedc6-f3b9-4021-ab31-69fbd6d1fae5`; DOI `10.64898/2026.07.23.740352`). Its internal source annotations include `Chen_lab` and `Sanes_GSE236566`; the latter is nested within the CELLxGENE resource and is not an independent external cohort. Public source records: [CELLxGENE](https://cellxgene.cziscience.com), [GSE220155](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220155), [GSE135922](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135922), [GSE230348](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230348), [GSE135092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135092), and [PXD080419](https://www.ebi.ac.uk/pride/archive/projects/PXD080419).

## Repository structure

- `scripts/`: Python and R analysis scripts.
- `environment/`: software and package versions recorded during the analysis.

The scripts write structured results to `results/` by default. Input data are read from `inputs/` unless an environment variable overrides a path.

## Script sequence

| Script | Analysis role |
|---|---|
| `phase_a_h5ad_audit.py` | Audits H5AD matrices, annotations, and donor-by-region cell counts. |
| `phase_b_aggregate_raw_counts.py` | Aggregates raw UMI counts into donor-by-region pseudobulk samples. |
| `phase_b_edger_ql.R` | Fits the paired donor-adjusted edgeR quasi-likelihood model. |
| `phase_bc_robustness.py` | Computes paired robust inference, bootstrap intervals, leave-one-donor-out summaries, and threshold, age, and source sensitivities. |
| `phase_d_state_aggregate.py` | Aggregates counts and normalized expression within annotated RPE states. |
| `phase_d_state_edger.R` | Fits state-specific paired edgeR models for eligible states. |
| `phase_de_within_state_decomposition.py` | Tests within-state persistence and decomposes regional effects. |
| `phase_f_gse220155_empirical_null.py` | Builds the expression- and detection-matched empirical null for GSE220155. |
| `phase_g_gse135922_external.py` | Performs the bounded GSE135922 RPE direction analysis. |
| `phase_h_gse230348_feasibility.py` | Audits whether GSE230348 supports an RPE-specific paired analysis. |
| `phase_i_gse135092_bulk.py` | Cross-checks the predefined gene set against GSE135092 bulk RPE results. |
| `phase_j_pxd080419_proteomics.py` | Assesses whether PXD080419 supports quantitative protein concordance. |
| `phase_k_l_lineage_integrated_evidence.py` | Builds data-lineage records and the integrated gene-level evidence table. |
| `phase_n_manifest_session.py` | Records input checksums, result checksums, and software-session information. |

## Input configuration

Set `M28_INPUT_DIR` to a directory containing the local inputs, or set individual variables when inputs are stored separately:

| Variable | Expected input |
|---|---|
| `M28_H5AD` | Discovery H5AD object. |
| `M28_DISCOVERY_RESULTS` | Discovery paired-expression results with gene identifiers, paired-test FDR, and regional mean difference. |
| `M28_GSE220155_DIR` | GSE220155 normalized expression files. |
| `M28_GSE135922_DIR` | GSE135922 processed expression tables. |
| `M28_GSE230348_DIR` | GSE230348 metadata and supplementary inventory. |
| `M28_GSE135092_METADATA` | Standardized GSE135092 sample metadata. |
| `M28_GSE135092_RAW` | GSE135092 processed-expression archive. |
| `M28_GSE135092_SUPPLEMENT` | Published GSE135092 supplementary workbook. |
| `M28_PXD080419_METADATA` | PXD080419 repository metadata in JSON format. |
| `M28_GSE220155_GENE_ACTIVITY` | Precomputed GSE220155 gene-activity support table. |
| `M28_RESULTS_DIR` | Output directory; defaults to `results/`. |

Default file names under `inputs/` are encoded near the top of each script. The input files themselves are not required to be stored inside the repository.

## Running the workflow

Python dependencies recorded for the analysis include Python 3.13, anndata, NumPy, pandas, SciPy, statsmodels, openpyxl, and tabulate. The count models use R 4.3 with edgeR and limma. Exact recorded versions are provided in `environment/`.

After configuring the inputs, run scripts in the order shown above. For example:

```bash
python scripts/phase_a_h5ad_audit.py
python scripts/phase_b_aggregate_raw_counts.py
Rscript scripts/phase_b_edger_ql.R
python scripts/phase_bc_robustness.py
```

Continue with the state and independent-dataset scripts after their prerequisite result tables are available.

## Statistical scope

The primary discovery comparison is donor paired. State-specific inference is restricted to states meeting the prespecified donor and cell-count thresholds. The GSE220155 analysis uses an expression- and detection-matched empirical null with a fixed random seed and 10,000 draws. Small independent datasets and tissue-level datasets are interpreted within their sample-size and annotation constraints. The workflow supports regional expression associations; it does not establish causality, disease mechanism, or clinical utility. GSE230348 is not a quantitative validation dataset because of incomplete reviewed sample coverage, and PXD080419 is not an assessable negative protein result when a usable processed abundance table and sample map are unavailable.

## License

The analysis code is available under the MIT License.
