#!/usr/bin/env python3
"""Write a checksum manifest and software/session record for the upgrade."""
from pathlib import Path
import hashlib
import os
import platform
import subprocess
import sys

import anndata, numpy, pandas, scipy, statsmodels
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = Path(os.environ.get("M28_INPUT_DIR", PROJECT_ROOT / "inputs"))
ROOT = Path(os.environ.get("M28_RESULTS_DIR", PROJECT_ROOT / "results"))
H5AD = Path(os.environ.get("M28_H5AD", INPUT_ROOT / "rpe_choroid.h5ad"))
DISCOVERY_RESULTS = Path(os.environ.get("M28_DISCOVERY_RESULTS", INPUT_ROOT / "discovery_paired_results.csv"))
GSE135922_RAW = Path(os.environ.get("M28_GSE135922_RAW", INPUT_ROOT / "GSE135922" / "GSE135922_RAW.tar"))
GSE135092_RAW = Path(os.environ.get("M28_GSE135092_RAW", INPUT_ROOT / "GSE135092" / "GSE135092_RAW.tar"))
GSE135092_SUPPLEMENT = Path(os.environ.get("M28_GSE135092_SUPPLEMENT", INPUT_ROOT / "GSE135092" / "mmc2.xlsx"))
PXD080419_METADATA = Path(os.environ.get("M28_PXD080419_METADATA", INPUT_ROOT / "PXD080419" / "metadata.json"))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()


def script_for(path):
    rel=str(path.relative_to(ROOT))
    mapping={
        "00_DATA_AUDIT/":"scripts/phase_a_h5ad_audit.py",
        "01_DONOR_ROBUSTNESS/":"scripts/phase_b_aggregate_raw_counts.py; scripts/phase_b_edger_ql.R; scripts/phase_bc_robustness.py",
        "02_RPE_STATE_DECOMPOSITION/":"scripts/phase_d_state_aggregate.py; scripts/phase_d_state_edger.R; scripts/phase_de_within_state_decomposition.py",
        "03_GSE220155_EMPIRICAL_BACKGROUND/":"scripts/phase_f_gse220155_empirical_null.py",
        "04_GSE135922_EXTERNAL/":"scripts/phase_g_gse135922_external.py",
        "05_GSE230348_FEASIBILITY/":"scripts/phase_h_gse230348_feasibility.py",
        "06_GSE135092_BULK/":"scripts/phase_i_gse135092_bulk.py",
        "07_PXD080419_PROTEOMICS/":"scripts/phase_j_pxd080419_proteomics.py",
        "08_INTEGRATED_EVIDENCE/":"scripts/phase_k_l_lineage_integrated_evidence.py",
    }
    for prefix,script in mapping.items():
        if rel.startswith(prefix): return script
    return "manual synthesis from audited outputs" if path.suffix==".md" else "scripts/phase_n_manifest_session.py"


def main():
    inputs=[
        ("discovery_h5ad",H5AD,"CELLxGENE collection 48c15b0c"),
        ("discovery_results",DISCOVERY_RESULTS,"discovery analysis"),
        ("GSE135922_RAW",GSE135922_RAW,"GSE135922"),
        ("GSE135092_RAW",GSE135092_RAW,"GSE135092"),
        ("GSE135092_Data_S1",GSE135092_SUPPLEMENT,"DOI 10.1016/j.celrep.2019.12.082"),
        ("PXD080419_metadata",PXD080419_METADATA,"PXD080419"),
    ]
    output_files=[p for p in ROOT.rglob("*") if p.is_file() and p.name not in {"MANIFEST.tsv","SOFTWARE_SESSION_INFO.txt"} and "__pycache__" not in p.parts]
    rows=[]
    for name,p,acc in inputs:
        rows.append({"record_type":"input","input":name,"accession_or_path":acc if acc else str(p),"path":str(p),"size_bytes":p.stat().st_size if p.exists() else "MISSING","sha256":sha(p) if p.exists() else "MISSING","script":"","output":""})
    for p in sorted(output_files):
        rows.append({"record_type":"output","input":"see producing script","accession_or_path":"","path":str(p),"size_bytes":p.stat().st_size,"sha256":sha(p),"script":script_for(p),"output":str(p.relative_to(ROOT))})
    pd.DataFrame(rows).to_csv(ROOT/"MANIFEST.tsv",sep="\t",index=False)

    try:
        rinfo=subprocess.run(["Rscript","-e","sessionInfo()"],capture_output=True,text=True,timeout=30).stdout
    except Exception as exc: rinfo=f"R sessionInfo unavailable: {exc}"
    text=f"""M28 analysis software/session information
Platform: {platform.platform()}
Python: {sys.version}
anndata: {anndata.__version__}
numpy: {numpy.__version__}
pandas: {pandas.__version__}
scipy: {scipy.__version__}
statsmodels: {statsmodels.__version__}
Random seed: 20260811
Bootstrap replicates: 10000
GSE220155 matched-null draws: 10000

R sessionInfo
-------------
{rinfo}
"""
    (ROOT/"SOFTWARE_SESSION_INFO.txt").write_text(text)
    manifest = pd.read_csv(ROOT/"MANIFEST.tsv", sep="\t")
    session_path = ROOT/"SOFTWARE_SESSION_INFO.txt"
    session_row = {"record_type":"output","input":"runtime environments","accession_or_path":"","path":str(session_path),"size_bytes":session_path.stat().st_size,"sha256":sha(session_path),"script":"scripts/phase_n_manifest_session.py","output":"SOFTWARE_SESSION_INFO.txt"}
    manifest = pd.concat([manifest, pd.DataFrame([session_row])], ignore_index=True)
    manifest.to_csv(ROOT/"MANIFEST.tsv", sep="\t", index=False)


if __name__=="__main__":main()
