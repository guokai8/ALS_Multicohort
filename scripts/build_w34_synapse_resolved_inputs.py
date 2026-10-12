#!/usr/bin/env python3
"""Build downstream inputs from the Synapse/NYGC-resolved W33 results.

W01/W03 remain frozen GEO-label analyses. W34 is a separate input layer so
that downstream regulatory evidence can be traced to corrected diagnoses.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from p01_paths import project_root


ROOT = project_root()
W33 = ROOT / "results" / "W33_synapse_resolved_diagnosis_sensitivity"
W01 = ROOT / "results" / "W01_RNA_DE"
W03 = ROOT / "results" / "W03_ATAC_DA"
OUT = ROOT / "results" / "W34_synapse_resolved_regulatory_reanalysis"
RNA_SEEDS = OUT / "inputs" / "rna_seed_lists"
PEAK_SETS = OUT / "inputs" / "peak_set_lists"
TABLES = OUT / "tables"


def read_table(path: Path) -> pd.DataFrame:
    with path.open("rb") as handle:
        magic = handle.read(2)
    compression = "gzip" if magic == b"\x1f\x8b" else None
    return pd.read_csv(path, sep="\t", compression=compression)


def write_table(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def build_rna_inputs() -> dict[str, int]:
    w33 = read_table(W33 / "tables" / "MTC_Astrocyte_RNA.synapse_resolved.tsv.gz")
    old = read_table(W01 / "tables" / "GSE330130_MTC_RNA__Astrocyte.limma_voom.tsv.gz")
    annotation_cols = [
        "gene_id", "common_name", "chrom", "logFC", "AveExpr", "t", "P.Value",
        "adj.P.Val", "B", "method", "model_formula", "covariates_used",
        "n_samples_model", "n_als_donors_model", "n_control_donors_model", "n_genes_tested",
    ]
    annotation = old[[col for col in annotation_cols if col in old.columns]].copy()
    annotation = annotation.drop_duplicates("gene_id")
    full = w33.rename(columns={"feature_id": "gene_id"}).merge(
        annotation.drop(columns=["logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"], errors="ignore"),
        on="gene_id", how="left", suffixes=("_w33", "_frozen")
    )
    full["analysis_label_source"] = "Synapse_NYGC_resolved_W33"
    full["primary_feature_rule"] = "W33_adj.P.Val<0.05_no_abslogFC"
    write_table(full, TABLES / "GSE330130_MTC_RNA__Astrocyte.W33_synapse_resolved.tsv.gz", "gzip")

    primary = full[full["adj.P.Val"] < 0.05].copy()
    counts: dict[str, int] = {}
    for direction, mask in (("up", primary["logFC"] > 0), ("down", primary["logFC"] < 0)):
        subset = primary.loc[mask].copy()
        counts[f"rna_{direction}"] = len(subset)
        subset["dataset_id"] = "GSE330130_MTC_RNA"
        subset["cell_type"] = "Astrocyte"
        subset["gene_id"] = subset["gene_id"].astype(str)
        write_table(subset, RNA_SEEDS / f"GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_{direction}.tsv")
    return counts


def build_atac_inputs() -> dict[str, int]:
    w33 = read_table(W33 / "tables" / "MTC_Astrocyte_ATAC.synapse_resolved.tsv.gz")
    old = read_table(W03 / "tables" / "GSE330130_MTC_ATAC__Astrocyte.DA_limma_voom.tsv.gz")
    coord_cols = [
        "peak_id", "id", "genome", "interval", "feature_type", "chrom", "start", "end",
        "common_name", "n_counts",
    ]
    coords = old[[col for col in coord_cols if col in old.columns]].drop_duplicates("peak_id")
    # W33 is the statistical authority for Synapse/NYGC-resolved tested peaks.
    # W03 annotation is auxiliary only; its feature-filtering universe differs
    # slightly from W33 and must not silently drop W33-tested peaks.
    full = w33.rename(columns={"feature_id": "peak_id"}).merge(
        coords,
        on="peak_id",
        how="left",
        suffixes=("", "_w03"),
    )
    parsed = full["peak_id"].astype(str).str.extract(
        r"^(?P<parsed_chrom>chr[^:]+):(?P<parsed_start>[0-9]+)-(?P<parsed_end>[0-9]+)$"
    )
    full["chrom"] = full["chrom"].where(full["chrom"].notna(), parsed["parsed_chrom"])
    full["start"] = full["start"].where(full["start"].notna(), pd.to_numeric(parsed["parsed_start"], errors="coerce"))
    full["end"] = full["end"].where(full["end"].notna(), pd.to_numeric(parsed["parsed_end"], errors="coerce"))
    full["id"] = full["id"].where(full["id"].notna(), full["peak_id"])
    full["interval"] = full["interval"].where(full["interval"].notna(), full["peak_id"])
    full["common_name"] = full["common_name"].where(full["common_name"].notna(), full["peak_id"])
    full["genome"] = full["genome"].where(full["genome"].notna(), "GRCh38")
    full["feature_type"] = full["feature_type"].where(full["feature_type"].notna(), "Peaks")
    full["coordinate_annotation_source"] = "W03_frozen_annotation"
    full.loc[~full["peak_id"].isin(set(coords["peak_id"].astype(str))), "coordinate_annotation_source"] = "parsed_peak_id"
    # Keep feature_id for downstream scripts that expect the W33 result schema.
    full["feature_id"] = full["peak_id"]
    front_cols = [
        "peak_id", "id", "genome", "interval", "feature_type", "chrom", "start", "end",
        "common_name", "n_counts", "coordinate_annotation_source",
    ]
    full = full[[col for col in front_cols if col in full.columns] + [
        col for col in full.columns if col not in front_cols
    ]]
    full["analysis_label_source"] = "Synapse_NYGC_resolved_W33"
    full["primary_feature_rule"] = "W33_adj.P.Val<0.05_no_abslogFC"
    write_table(full, TABLES / "GSE330130_MTC_ATAC__Astrocyte.W33_synapse_resolved.tsv.gz", "gzip")

    w33_ids = set(w33["feature_id"].astype(str))
    old_ids = set(coords["peak_id"].astype(str))
    full_ids = set(full["peak_id"].astype(str))
    discrepancy_rows = [
        {
            "comparison": "W33_tested_not_in_W03_annotation",
            "n_peaks": len(w33_ids - old_ids),
            "n_fdr05": int(w33[w33["feature_id"].astype(str).isin(w33_ids - old_ids)]["adj.P.Val"].lt(0.05).sum()),
            "min_adj_p": w33[w33["feature_id"].astype(str).isin(w33_ids - old_ids)]["adj.P.Val"].min(),
            "interpretation": "W33-tested peaks retained in W34; coordinates parsed from peak_id; no FDR05 primary peaks affected.",
        },
        {
            "comparison": "W03_annotation_not_in_W33_tested",
            "n_peaks": len(old_ids - w33_ids),
            "n_fdr05": int(old[old["peak_id"].astype(str).isin(old_ids - w33_ids)]["adj.P.Val"].lt(0.05).sum()),
            "min_adj_p": old[old["peak_id"].astype(str).isin(old_ids - w33_ids)]["adj.P.Val"].min(),
            "interpretation": "Old W03-only peaks are not part of the Synapse/NYGC-resolved W33 tested universe.",
        },
        {
            "comparison": "W34_final_tested_universe",
            "n_peaks": len(full_ids),
            "n_fdr05": int(full["adj.P.Val"].lt(0.05).sum()),
            "min_adj_p": full["adj.P.Val"].min(),
            "interpretation": "Canonical W34 ATAC tested universe equals W33 full-model tested universe.",
        },
    ]
    write_table(pd.DataFrame(discrepancy_rows), OUT / "w34_atac_feature_universe_discrepancy_audit.tsv")

    primary = full[full["adj.P.Val"] < 0.05].copy()
    counts: dict[str, int] = {}
    for direction, mask, token in (("open", primary["logFC"] > 0, "open_in_ALS"),
                                   ("closed", primary["logFC"] < 0, "closed_in_ALS")):
        subset = primary.loc[mask].copy()
        counts[f"atac_{direction}"] = len(subset)
        subset["dataset_id"] = "GSE330130_MTC_ATAC"
        subset["cell_type"] = "Astrocyte"
        write_table(subset, PEAK_SETS / f"GSE330130_MTC_ATAC__Astrocyte.FDR05_no_abslogFC_{token}.tsv")
        write_table(subset, PEAK_SETS / f"GSE330130_MTC_ATAC__Astrocyte.FDR05_logFC025_{token}.tsv")
    return counts


def main() -> None:
    for directory in (RNA_SEEDS, PEAK_SETS, TABLES):
        directory.mkdir(parents=True, exist_ok=True)
    counts = {**build_rna_inputs(), **build_atac_inputs()}
    audit = pd.DataFrame([
        {"input_layer": "W34", "source": "W33_synapse_resolved_diagnosis_sensitivity", "metric": key, "value": value}
        for key, value in counts.items()
    ])
    audit.to_csv(OUT / "w34_input_build_audit.tsv", sep="\t", index=False)
    (OUT / "README.md").write_text(
        "# W34 Synapse/NYGC-resolved regulatory reanalysis\n\n"
        "W34 is a separate downstream input layer built from W33 resolved diagnosis results. "
        "W01/W03 remain frozen GEO-label analyses and are retained for sensitivity comparison.\n\n"
        "Primary features use W33 adjusted P<0.05 without an absolute logFC cutoff.\n",
        encoding="utf-8",
    )
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
