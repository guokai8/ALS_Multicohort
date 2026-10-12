#!/usr/bin/env python3
"""Run W08 Pineda ALS/FTLD cortical RNA module validation.

This workflow projects W01 RNA seed modules into the Pineda_snRNAseq processed
sample folders. It streams MatrixMarket files sample-by-sample and writes
sample/donor-level RNA module validation only. It does not create ATAC,
regulatory, or splicing claims.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from p01_paths import project_root, source_data_root


ROOT = project_root()
SOURCE_DATA_ROOT = source_data_root()
PINEDA_DIR = SOURCE_DATA_ROOT / "Other_Datasets/Pineda_snRNAseq/synapse_download/Processed Data"
MANIFEST = ROOT / "metadata/pineda_sample_manifest.tsv"
W01_SEEDS = ROOT / "results/W01_RNA_DE/module_seed_lists"
OUT = ROOT / "results/W08_Pineda_validation"

MIN_MODULE_GENES = 3
MIN_CELLS_PER_CELLTYPE_SAMPLE = 20
PRIMARY_REGIONS = {"BA4", "BA9"}


def safe_name(value: object) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value))
    return re.sub(r"_+", "_", text).strip("_") or "NA"


def parse_seed_file(path: Path) -> tuple[str, str, str]:
    name = path.name.replace(".tsv", "")
    prefix, direction = name.split(".FDR05_logFC025_")
    dataset_id, cell_type_safe = prefix.split("__", 1)
    return dataset_id, cell_type_safe.replace("_", " "), direction


def resolve_source_path(value: object) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else SOURCE_DATA_ROOT / path


def source_path_id(path: Path) -> str:
    try:
        return path.relative_to(SOURCE_DATA_ROOT).as_posix()
    except ValueError:
        return path.name


def broad_celltype(row: pd.Series) -> str:
    cell_type = str(row.get("CellType", "")).strip()
    cell_class = str(row.get("CellClass", "")).strip()
    if cell_type == "Astro":
        return "Astrocyte"
    if cell_type == "Micro":
        return "Microglia"
    if cell_type == "OPC":
        return "OPC"
    if cell_type == "Oligo":
        return "Oligodendrocyte"
    if cell_class == "Ex":
        return "Excitatory neuron"
    if cell_class == "In":
        return "Inhibitory neuron"
    if cell_type == "Endo":
        return "Endothelial"
    if cell_type in {"Mural", "Fibro"}:
        return "Vascular/Meningeal"
    if cell_type == "T_Cell":
        return "T cell"
    return cell_type or "Unknown"


def load_modules(gene_universe: set[str]) -> tuple[pd.DataFrame, list[dict[str, object]], pd.DataFrame]:
    modules: list[dict[str, object]] = []
    skipped: list[dict[str, object]] = []

    for path in sorted(W01_SEEDS.glob("*.FDR05_logFC025_*.tsv")):
        dataset_id, cell_type, direction = parse_seed_file(path)
        df = pd.read_csv(path, sep="\t")
        genes = sorted(set(df["gene_id"].dropna().astype(str))) if "gene_id" in df.columns else []
        genes_in_pineda = [gene for gene in genes if gene in gene_universe]
        module_id = f"{dataset_id}__{safe_name(cell_type)}__{direction}"
        record = {
            "module_id": module_id,
            "source_dataset_id": dataset_id,
            "source_cell_type": cell_type,
            "direction": direction,
            "n_seed_genes": len(genes),
            "n_genes_in_pineda": len(genes_in_pineda),
            "seed_path": str(path),
        }
        if len(genes_in_pineda) < MIN_MODULE_GENES:
            skipped.append({**record, "skip_reason": f"module_genes_below_min_{MIN_MODULE_GENES}"})
            continue
        modules.append({**record, "genes": genes_in_pineda})

    registry = pd.DataFrame([{k: v for k, v in module.items() if k != "genes"} for module in modules])
    skipped_df = pd.DataFrame(skipped)
    return registry, modules, skipped_df


def first_existing_row_metadata(manifest: pd.DataFrame) -> Path:
    for sample_path in manifest["path"].astype(str):
        candidate = resolve_source_path(sample_path) / "row_metadata.tsv"
        if candidate.exists():
            return candidate
    raise FileNotFoundError("No row_metadata.tsv found from Pineda manifest paths")


def read_matrix_market_selected(
    matrix_path: Path,
    selected_rows: dict[int, int],
    n_selected: int,
) -> tuple[np.ndarray, np.ndarray, tuple[int, int, int]]:
    """Read total cell counts and selected-gene counts from MatrixMarket.

    MatrixMarket row/column indices are one-based. Pineda matrices are genes x
    cells, so selected_rows maps one-based gene row index to a compact row index.
    """

    opener = gzip.open if matrix_path.suffix == ".gz" else open
    with opener(matrix_path, "rt") as handle:
        header = handle.readline()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"{matrix_path} is not a MatrixMarket file")
        dims_line = ""
        for line in handle:
            if not line.startswith("%"):
                dims_line = line.strip()
                break
        if not dims_line:
            raise ValueError(f"No MatrixMarket dimension line in {matrix_path}")
        n_rows, n_cols, n_nonzero = map(int, dims_line.split()[:3])
        totals = np.zeros(n_cols, dtype=np.float64)
        selected = np.zeros((n_selected, n_cols), dtype=np.float32)

        for line in handle:
            if not line.strip():
                continue
            row_s, col_s, value_s = line.split()[:3]
            row_i = int(row_s)
            col_i = int(col_s) - 1
            value = float(value_s)
            totals[col_i] += value
            compact_i = selected_rows.get(row_i)
            if compact_i is not None:
                selected[compact_i, col_i] += value

    return totals, selected, (n_rows, n_cols, n_nonzero)


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    mask = values.notna() & weights.notna() & (weights > 0)
    if not mask.any():
        return float("nan")
    return float(np.average(values[mask].astype(float), weights=weights[mask].astype(float)))


def project_sample(
    sample_path: Path,
    modules: list[dict[str, object]],
    all_module_genes: list[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    row_meta = pd.read_csv(sample_path / "row_metadata.tsv", sep="\t")
    col_meta = pd.read_csv(sample_path / "col_metadata.tsv", sep="\t")
    col_meta["broad_cell_type"] = col_meta.apply(broad_celltype, axis=1)

    gene_to_row = {str(gene): i + 1 for i, gene in enumerate(row_meta["ENSEMBL"].astype(str))}
    selected_gene_rows: dict[int, int] = {}
    selected_gene_names: list[str] = []
    for gene in all_module_genes:
        row_idx = gene_to_row.get(gene)
        if row_idx is None:
            continue
        selected_gene_rows[row_idx] = len(selected_gene_names)
        selected_gene_names.append(gene)

    totals, selected_counts, dims = read_matrix_market_selected(
        sample_path / "counts_fil.mtx",
        selected_gene_rows,
        len(selected_gene_names),
    )
    if dims[1] != col_meta.shape[0]:
        raise ValueError(f"Cell count mismatch in {sample_path}: mtx={dims[1]}, metadata={col_meta.shape[0]}")

    gene_to_selected = {gene: i for i, gene in enumerate(selected_gene_names)}
    rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []

    sample_meta_cols = ["Sample_ID", "Donor", "Region", "Sex", "Condition", "Group", "C9_pos"]
    sample_meta = {col: col_meta[col].iloc[0] if col in col_meta.columns else np.nan for col in sample_meta_cols}

    for target_cell_type, cell_idx in col_meta.groupby("broad_cell_type", observed=True).indices.items():
        cell_idx_arr = np.asarray(cell_idx, dtype=int)
        valid_cells = cell_idx_arr[totals[cell_idx_arr] > 0]
        n_cells = int(valid_cells.size)
        audit_rows.append(
            {
                "sample_folder": sample_path.name,
                **sample_meta,
                "target_cell_type": target_cell_type,
                "n_cells": n_cells,
                "total_counts": float(totals[valid_cells].sum()) if n_cells else 0.0,
            }
        )
        if n_cells < MIN_CELLS_PER_CELLTYPE_SAMPLE:
            continue

        cell_totals = totals[valid_cells]
        for module in modules:
            if target_cell_type != module["source_cell_type"]:
                continue
            module_gene_idx = [gene_to_selected[gene] for gene in module["genes"] if gene in gene_to_selected]
            if len(module_gene_idx) < MIN_MODULE_GENES:
                continue
            raw = selected_counts[np.asarray(module_gene_idx, dtype=int)[:, None], valid_cells[None, :]]
            cpm = raw.astype(np.float64) / cell_totals[None, :] * 1_000_000.0
            per_cell_score = np.nanmean(np.log2(cpm + 1.0), axis=0)
            rows.append(
                {
                    "sample_folder": sample_path.name,
                    "sample_path": source_path_id(sample_path),
                    **sample_meta,
                    "target_cell_type": target_cell_type,
                    "module_id": module["module_id"],
                    "source_dataset_id": module["source_dataset_id"],
                    "source_cell_type": module["source_cell_type"],
                    "direction": module["direction"],
                    "n_seed_genes": module["n_seed_genes"],
                    "n_genes_in_pineda": module["n_genes_in_pineda"],
                    "n_genes_scored_in_sample": len(module_gene_idx),
                    "n_cells": n_cells,
                    "mean_library_size": float(np.mean(cell_totals)),
                    "module_score": float(np.nanmean(per_cell_score)),
                    "module_score_sd_cells": float(np.nanstd(per_cell_score, ddof=1)) if n_cells > 1 else np.nan,
                }
            )

    return rows, audit_rows


def normalize_c9(value: object) -> object:
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    if text in {"", "NA", "N/A", "nan", "None"}:
        return np.nan
    try:
        return int(float(text))
    except ValueError:
        return np.nan


def aggregate_scores(sample_scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = sample_scores.copy()
    scores["C9_pos"] = scores["C9_pos"].map(normalize_c9)
    scores["Region"] = scores["Region"].astype(str)

    grouping = [
        "Donor",
        "Region",
        "target_cell_type",
        "module_id",
        "source_dataset_id",
        "source_cell_type",
        "direction",
        "n_seed_genes",
        "n_genes_in_pineda",
    ]
    first_cols = ["Sex", "Condition", "Group", "C9_pos"]
    donor_region = (
        scores.groupby(grouping, dropna=False, observed=True)
        .apply(
            lambda df: pd.Series(
                {
                    **{col: df[col].iloc[0] for col in first_cols},
                    "n_samples": int(df["sample_folder"].nunique()),
                    "n_cells": int(df["n_cells"].sum()),
                    "module_score": weighted_mean(df["module_score"], df["n_cells"]),
                    "mean_library_size": weighted_mean(df["mean_library_size"], df["n_cells"]),
                }
            )
        )
        .reset_index()
    )

    donor_mean_input = donor_region[donor_region["Region"].isin(PRIMARY_REGIONS)].copy()
    grouping_mean = [
        "Donor",
        "target_cell_type",
        "module_id",
        "source_dataset_id",
        "source_cell_type",
        "direction",
        "n_seed_genes",
        "n_genes_in_pineda",
    ]
    donor_mean = (
        donor_mean_input.groupby(grouping_mean, dropna=False, observed=True)
        .apply(
            lambda df: pd.Series(
                {
                    **{col: df[col].iloc[0] for col in first_cols},
                    "region_scope": "BA4_BA9_mean",
                    "regions_included": ";".join(sorted(df["Region"].astype(str).unique())),
                    "n_regions": int(df["Region"].nunique()),
                    "n_samples": int(df["n_samples"].sum()),
                    "n_cells": int(df["n_cells"].sum()),
                    "module_score": weighted_mean(df["module_score"], df["n_cells"]),
                    "mean_library_size": weighted_mean(df["mean_library_size"], df["n_cells"]),
                }
            )
        )
        .reset_index()
    )
    return donor_region, donor_mean


def add_binary_predictor(df: pd.DataFrame, contrast: str) -> tuple[pd.DataFrame, str]:
    model_df = df.copy()
    if contrast == "ALS_vs_PN":
        model_df = model_df[model_df["Condition"].isin(["ALS", "PN"])].copy()
        model_df["predictor_value"] = (model_df["Condition"] == "ALS").astype(float)
        return model_df, "Condition: ALS vs PN"
    if contrast == "FTLD_vs_PN":
        model_df = model_df[model_df["Condition"].isin(["FTLD", "PN"])].copy()
        model_df["predictor_value"] = (model_df["Condition"] == "FTLD").astype(float)
        return model_df, "Condition: FTLD vs PN"
    if contrast == "ALS_FTLD_vs_PN":
        model_df = model_df[model_df["Condition"].isin(["ALS", "FTLD", "PN"])].copy()
        model_df["predictor_value"] = model_df["Condition"].isin(["ALS", "FTLD"]).astype(float)
        return model_df, "Condition: ALS/FTLD vs PN"
    if contrast in {"SALS_vs_PN", "C9ALS_vs_PN", "SFTLD_vs_PN", "C9FTLD_vs_PN"}:
        group = contrast.replace("_vs_PN", "")
        model_df = model_df[model_df["Group"].isin([group, "PN"])].copy()
        model_df["predictor_value"] = (model_df["Group"] == group).astype(float)
        return model_df, f"Group: {group} vs PN"
    if contrast == "C9ALS_vs_SALS":
        model_df = model_df[model_df["Group"].isin(["C9ALS", "SALS"])].copy()
        model_df["predictor_value"] = (model_df["Group"] == "C9ALS").astype(float)
        return model_df, "Group: C9ALS vs SALS"
    if contrast == "C9FTLD_vs_SFTLD":
        model_df = model_df[model_df["Group"].isin(["C9FTLD", "SFTLD"])].copy()
        model_df["predictor_value"] = (model_df["Group"] == "C9FTLD").astype(float)
        return model_df, "Group: C9FTLD vs SFTLD"
    raise ValueError(f"Unknown contrast: {contrast}")


def direction_consistent(direction: str, beta: float) -> object:
    if math.isnan(beta):
        return np.nan
    if direction == "up":
        return beta > 0
    if direction == "down":
        return beta < 0
    return np.nan


def fit_contrast(
    df: pd.DataFrame,
    contrast: str,
    region_scope: str,
    min_group_donors: int = 5,
    min_total_donors: int = 10,
) -> dict[str, object]:
    model_df, contrast_label = add_binary_predictor(df, contrast)
    result = {
        "contrast": contrast,
        "contrast_label": contrast_label,
        "region_scope": region_scope,
        "status": "not_run",
        "n_observations": int(model_df.shape[0]),
        "n_donors": int(model_df["Donor"].nunique()) if not model_df.empty else 0,
        "n_cells_model_total": int(model_df["n_cells"].sum()) if "n_cells" in model_df.columns and not model_df.empty else 0,
        "n_samples_model_total": int(model_df["n_samples"].sum()) if "n_samples" in model_df.columns and not model_df.empty else 0,
        "predictor_min_group_donors": np.nan,
        "model_formula_registered": "module_score ~ predictor_value + sex_if_estimable",
        "covariates_used_registered": "",
        "covariates_dropped": "",
        "beta": np.nan,
        "se": np.nan,
        "p_value": np.nan,
        "message": "",
    }
    if model_df.empty:
        result["status"] = "skipped"
        result["message"] = "empty_after_contrast_filter"
        return result
    group_counts = model_df.groupby("predictor_value")["Donor"].nunique()
    result["predictor_min_group_donors"] = int(group_counts.min()) if not group_counts.empty else 0
    if len(group_counts) < 2 or group_counts.min() < min_group_donors:
        result["status"] = "skipped"
        result["message"] = "insufficient_predictor_group_donors"
        return result
    if model_df["Donor"].nunique() < min_total_donors:
        result["status"] = "skipped"
        result["message"] = "insufficient_total_donors"
        return result

    x_cols = ["predictor_value"]
    dropped: list[str] = []
    used: list[str] = []
    if "Sex" in model_df.columns and model_df["Sex"].notna().all() and model_df["Sex"].astype(str).nunique() >= 2:
        dummies = pd.get_dummies(model_df["Sex"].astype(str), prefix="Sex", drop_first=True, dtype=float)
        for col in dummies.columns:
            model_df[col] = dummies[col]
        x_cols.extend(dummies.columns.tolist())
        used.append("Sex")
    else:
        dropped.append("Sex:not_estimable")

    model_df = model_df[["module_score", "predictor_value", "Donor", *x_cols[1:]]].dropna().copy()
    x = sm.add_constant(model_df[x_cols], has_constant="add")
    rank = np.linalg.matrix_rank(x.to_numpy(dtype=float))
    if rank < x.shape[1]:
        for col in list(x.columns)[::-1]:
            if col in {"const", "predictor_value"}:
                continue
            x = x.drop(columns=[col])
            dropped.append(f"{col}:rank_deficient")
            if np.linalg.matrix_rank(x.to_numpy(dtype=float)) == x.shape[1]:
                break
    if np.linalg.matrix_rank(x.to_numpy(dtype=float)) < x.shape[1]:
        result["status"] = "skipped"
        result["message"] = "rank_deficient_after_covariate_drop"
        result["covariates_used_registered"] = ";".join(used)
        result["covariates_dropped"] = "|".join(dropped)
        return result

    try:
        fit = sm.OLS(model_df["module_score"].astype(float), x.astype(float)).fit()
        result.update(
            {
                "status": "completed",
                "beta": float(fit.params["predictor_value"]),
                "se": float(fit.bse["predictor_value"]),
                "p_value": float(fit.pvalues["predictor_value"]),
                "covariates_used_registered": ";".join(used),
                "covariates_dropped": "|".join(dropped),
                "message": "completed_OLS_donor_level",
            }
        )
    except Exception as exc:  # noqa: BLE001
        result["status"] = "error"
        result["message"] = str(exc)
    return result


def run_models(donor_mean: pd.DataFrame, donor_region: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    meta_cols = [
        "module_id",
        "source_dataset_id",
        "source_cell_type",
        "target_cell_type",
        "direction",
        "n_seed_genes",
        "n_genes_in_pineda",
    ]
    contrasts = [
        "ALS_vs_PN",
        "FTLD_vs_PN",
        "ALS_FTLD_vs_PN",
        "SALS_vs_PN",
        "C9ALS_vs_PN",
        "SFTLD_vs_PN",
        "C9FTLD_vs_PN",
        "C9ALS_vs_SALS",
        "C9FTLD_vs_SFTLD",
    ]

    for module_id, module_df in donor_mean.groupby("module_id", observed=True):
        meta = module_df.iloc[0][[col for col in meta_cols if col in module_df.columns]].to_dict()
        for contrast in contrasts:
            result = fit_contrast(module_df, contrast, "BA4_BA9_mean")
            rows.append({**meta, **result, "analysis_type": "donor_mean_primary"})

    region_input = donor_region[donor_region["Region"].isin(PRIMARY_REGIONS)].copy()
    for (module_id, region), module_df in region_input.groupby(["module_id", "Region"], observed=True):
        meta = module_df.iloc[0][[col for col in meta_cols if col in module_df.columns]].to_dict()
        meta["Region"] = region
        for contrast in ["ALS_vs_PN", "FTLD_vs_PN", "ALS_FTLD_vs_PN"]:
            result = fit_contrast(module_df, contrast, str(region), min_group_donors=5, min_total_donors=10)
            rows.append({**meta, **result, "analysis_type": "region_specific_sensitivity"})

    results = pd.DataFrame(rows)
    if results.empty:
        return results
    completed = results["status"].eq("completed")
    results["adj_p_global_BH"] = np.nan
    results["direction_consistent_with_source_ALS"] = pd.Series([pd.NA] * len(results), dtype="object")
    if completed.any():
        results.loc[completed, "adj_p_global_BH"] = multipletests(results.loc[completed, "p_value"], method="fdr_bh")[1]
        als_mask = completed & results["contrast"].isin(["ALS_vs_PN", "SALS_vs_PN", "C9ALS_vs_PN", "ALS_FTLD_vs_PN"])
        results.loc[als_mask, "direction_consistent_with_source_ALS"] = [
            direction_consistent(direction, beta)
            for direction, beta in zip(results.loc[als_mask, "direction"], results.loc[als_mask, "beta"])
        ]
    return results


def validation_tier(row: pd.Series) -> str:
    if row["status"] != "completed":
        return "not_tested"
    if row["contrast"] != "ALS_vs_PN" or row["region_scope"] != "BA4_BA9_mean":
        return "secondary_context"
    consistent = bool(row["direction_consistent_with_source_ALS"])
    if consistent and pd.notna(row["adj_p_global_BH"]) and row["adj_p_global_BH"] < 0.05:
        return "cross_cohort_validated_BH05_directional"
    if consistent and row["p_value"] < 0.05:
        return "nominal_directional_support"
    if consistent:
        return "direction_only_no_significance"
    return "not_supported_or_opposite_direction"


def write_harmonization() -> None:
    rows = [
        ("CellType=Astro", "Astrocyte"),
        ("CellType=Micro", "Microglia"),
        ("CellType=OPC", "OPC"),
        ("CellType=Oligo", "Oligodendrocyte"),
        ("CellClass=Ex", "Excitatory neuron"),
        ("CellClass=In", "Inhibitory neuron"),
        ("CellType=Endo", "Endothelial"),
        ("CellType=Mural/Fibro", "Vascular/Meningeal"),
        ("CellType=T_Cell", "T cell"),
    ]
    pd.DataFrame(rows, columns=["pineda_rule", "harmonized_broad_cell_type"]).to_csv(
        OUT / "w08_celltype_harmonization.tsv", sep="\t", index=False
    )


def iter_manifest_rows(manifest: pd.DataFrame) -> Iterable[tuple[int, pd.Series]]:
    for idx, row in manifest.reset_index(drop=True).iterrows():
        yield idx, row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reuse-sample-scores",
        action="store_true",
        help="Reuse existing w08_sample_module_scores.tsv.gz and recompute donor/model tables only.",
    )
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    write_harmonization()

    manifest = pd.read_csv(MANIFEST, sep="\t")
    manifest = manifest[manifest["counts_fil_mtx_exists"].astype(str).eq("True")].copy()
    row_meta_path = first_existing_row_metadata(manifest)
    gene_universe = set(pd.read_csv(row_meta_path, sep="\t")["ENSEMBL"].dropna().astype(str))

    registry, modules, skipped = load_modules(gene_universe)
    registry.to_csv(OUT / "w08_module_registry.tsv", sep="\t", index=False)
    skipped.to_csv(OUT / "w08_skipped_modules.tsv", sep="\t", index=False)
    if registry.empty:
        raise RuntimeError("No W01 modules with at least three genes in Pineda")

    all_module_genes = sorted({gene for module in modules for gene in module["genes"]})
    sample_score_path = OUT / "w08_sample_module_scores.tsv.gz"
    if args.reuse_sample_scores and sample_score_path.exists():
        sample_scores = pd.read_csv(sample_score_path, sep="\t", compression="gzip")
        print(f"[W08] reused existing sample scores: {sample_score_path}", flush=True)
    else:
        sample_rows: list[dict[str, object]] = []
        audit_rows: list[dict[str, object]] = []
        failures: list[dict[str, object]] = []

        for idx, row in iter_manifest_rows(manifest):
            sample_path = resolve_source_path(row["path"])
            try:
                rows, audits = project_sample(sample_path, modules, all_module_genes)
                sample_rows.extend(rows)
                audit_rows.extend(audits)
                print(f"[W08] processed {idx + 1}/{len(manifest)} {sample_path.name}: score_rows={len(rows)}", flush=True)
            except Exception as exc:  # noqa: BLE001
                failures.append(
                    {
                        "sample_folder": row.get("sample_folder", sample_path.name),
                        "path": source_path_id(sample_path),
                        "error": str(exc),
                    }
                )
                print(f"[W08] failed {idx + 1}/{len(manifest)} {sample_path.name}: {exc}", flush=True)

        sample_scores = pd.DataFrame(sample_rows)
        sample_scores.to_csv(sample_score_path, sep="\t", index=False, compression="gzip")
        pd.DataFrame(audit_rows).to_csv(OUT / "w08_celltype_sample_coverage.tsv", sep="\t", index=False)
        pd.DataFrame(failures, columns=["sample_folder", "path", "error"]).to_csv(
            OUT / "w08_failed_samples.tsv", sep="\t", index=False
        )

    if sample_scores.empty:
        raise RuntimeError("No Pineda module scores were generated")

    donor_region, donor_mean = aggregate_scores(sample_scores)
    donor_region.to_csv(OUT / "w08_donor_region_module_scores.tsv", sep="\t", index=False)
    donor_mean.to_csv(OUT / "w08_donor_mean_module_scores.tsv", sep="\t", index=False)

    results = run_models(donor_mean, donor_region)
    if not results.empty:
        results["validation_tier"] = results.apply(validation_tier, axis=1)
    front = [
        "analysis_type",
        "region_scope",
        "Region",
        "module_id",
        "source_dataset_id",
        "source_cell_type",
        "direction",
        "contrast",
        "contrast_label",
        "validation_tier",
        "status",
        "beta",
        "se",
        "p_value",
        "adj_p_global_BH",
        "direction_consistent_with_source_ALS",
        "n_observations",
        "n_donors",
        "n_cells_model_total",
        "n_samples_model_total",
        "predictor_min_group_donors",
        "model_formula_registered",
        "covariates_used_registered",
        "covariates_dropped",
        "message",
        "n_seed_genes",
        "n_genes_in_pineda",
    ]
    present = [col for col in front if col in results.columns]
    results = results[present + [col for col in results.columns if col not in present]]
    results.to_csv(OUT / "w08_validation_results.tsv", sep="\t", index=False)

    audit = (
        results.groupby(["analysis_type", "region_scope", "status"], dropna=False)
        .size()
        .reset_index(name="n_tests")
        .sort_values(["analysis_type", "region_scope", "status"])
    )
    audit.to_csv(OUT / "w08_model_audit_summary.tsv", sep="\t", index=False)

    primary = results[
        results["analysis_type"].eq("donor_mean_primary")
        & results["region_scope"].eq("BA4_BA9_mean")
        & results["contrast"].eq("ALS_vs_PN")
    ].copy()
    primary = primary.sort_values(["validation_tier", "p_value"], na_position="last")
    primary.to_csv(OUT / "w08_primary_als_validation_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
