#!/usr/bin/env python3
"""Run W50 Pineda motor-cortex validation of the final W34 astrocyte modules.

W08 projected the historical W01 seed lists into Pineda. W50 projects the final
W34 Synapse/NYGC-resolved MTC astrocyte modules (15 ALS-up, 14 ALS-down genes)
into Pineda astrocytes, using the same per-cell score definition as W08/W09:
mean log2(CPM + 1) over module genes per nucleus, averaged within sample and
donor-region with nucleus-count weights.

Analysis plan, fixed before results were inspected:
- Primary: BA4 (primary motor cortex), ALS vs PN, OLS donor-level
  module_score ~ ALS + sex; BH across the up and down modules (as in W09).
- Sensitivities: + age at death + PMI (complete case); + log10 mean nucleus
  library size; donor pseudobulk summarization; SALS and C9ALS subgroups;
  BA4/BA9 donor mean (W08 scope); BA9 only; FTLD vs PN as context;
  leave-one-gene-out and gene-level directions in BA4.
The planned "BA4 plus M1 controls" sensitivity was dropped: the five M1 PN
samples (220616_BL6_MCX_*, donors MS1-MS5) are C57BL/6 mouse motor cortex and
share no human ENSG module genes.
Pineda is RNA-only, so W50 validates the RNA module only.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from p01_paths import project_root, source_data_root
from run_w08_pineda_validation import broad_celltype, resolve_source_path

ROOT = project_root()
SOURCE_DATA_ROOT = source_data_root()
MANIFEST = ROOT / "metadata/pineda_sample_manifest.tsv"
SEEDS = ROOT / "results/W34_synapse_resolved_regulatory_reanalysis/inputs/rna_seed_lists"
CLINICAL = ROOT / "results/W22_Pineda_S01_clinical_pathology_context/w22_pineda_lead_module_donor_context.tsv"
OUT = ROOT / "results/W50_Pineda_final15_validation"

MODULE_FILES = {
    "up": SEEDS / "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv",
    "down": SEEDS / "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_down.tsv",
}
TARGET_CELL_TYPE = "Astrocyte"
MIN_CELLS_PER_SAMPLE = 20
CHUNK_ROWS = 20_000_000


def load_modules() -> dict[str, pd.DataFrame]:
    modules = {}
    for direction, path in MODULE_FILES.items():
        df = pd.read_csv(path, sep="\t")[["gene_id", "common_name", "logFC"]].drop_duplicates("gene_id")
        modules[direction] = df.reset_index(drop=True)
    return modules


def read_astrocyte_counts(
    matrix_path: Path, gene_rows: np.ndarray, astro_cols: np.ndarray, n_cols: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return per-cell totals and module-gene counts for astrocyte columns only.

    MatrixMarket indices are one-based; Pineda matrices are genes x cells.
    gene_rows holds one-based row indices of module genes in compact order.
    """
    col_to_compact = np.full(n_cols + 1, -1, dtype=np.int64)
    col_to_compact[astro_cols + 1] = np.arange(astro_cols.size)
    row_lookup = np.full(int(gene_rows.max()) + 2 if gene_rows.size else 2, -1, dtype=np.int64)
    row_lookup[gene_rows] = np.arange(gene_rows.size)

    n_header = 0
    with open(matrix_path) as handle:
        for line in handle:
            n_header += 1
            if not line.startswith("%"):
                if int(line.split()[1]) != n_cols:
                    raise ValueError(f"Cell count mismatch in {matrix_path}: mtx={line.split()[1]}, metadata={n_cols}")
                break

    totals = np.zeros(astro_cols.size, dtype=np.float64)
    selected = np.zeros((gene_rows.size, astro_cols.size), dtype=np.float64)
    reader = pd.read_csv(
        matrix_path, sep=" ", header=None, skiprows=n_header,
        names=["row", "col", "value"], dtype={"row": np.int64, "col": np.int64, "value": np.float64},
        chunksize=CHUNK_ROWS, engine="c",
    )
    for chunk in reader:
        cols = col_to_compact[chunk["col"].to_numpy()]
        keep = cols >= 0
        if not keep.any():
            continue
        rows = chunk["row"].to_numpy()[keep]
        vals = chunk["value"].to_numpy()[keep]
        cols = cols[keep]
        totals += np.bincount(cols, weights=vals, minlength=astro_cols.size)
        in_range = rows < row_lookup.size
        gene_idx = np.full(rows.size, -1, dtype=np.int64)
        gene_idx[in_range] = row_lookup[rows[in_range]]
        hit = gene_idx >= 0
        np.add.at(selected, (gene_idx[hit], cols[hit]), vals[hit])
    return totals, selected


def project_sample(args: tuple[str, dict[str, list[str]]]) -> dict[str, Any]:
    sample_dir, module_genes = args
    sample_path = Path(sample_dir)
    row_meta = pd.read_csv(sample_path / "row_metadata.tsv", sep="\t")
    col_meta = pd.read_csv(sample_path / "col_metadata.tsv", sep="\t")
    col_meta["broad_cell_type"] = col_meta.apply(broad_celltype, axis=1)
    meta_cols = ["Sample_ID", "Donor", "Region", "Sex", "Condition", "Group", "C9_pos"]
    meta = {c: col_meta[c].iloc[0] if c in col_meta.columns else np.nan for c in meta_cols}

    gene_to_row = {g: i + 1 for i, g in enumerate(row_meta["ENSEMBL"].astype(str))}
    all_genes = sorted({g for genes in module_genes.values() for g in genes if g in gene_to_row})
    gene_rows = np.array([gene_to_row[g] for g in all_genes], dtype=np.int64)
    astro_cols = np.flatnonzero(col_meta["broad_cell_type"].eq(TARGET_CELL_TYPE).to_numpy())

    record: dict[str, Any] = {"sample_folder": sample_path.name, **meta, "n_astro_barcodes": int(astro_cols.size)}
    if astro_cols.size == 0:
        return {**record, "n_cells": 0, "genes": all_genes}
    totals, selected = read_astrocyte_counts(sample_path / "counts_fil.mtx", gene_rows, astro_cols, col_meta.shape[0])
    valid = totals > 0
    record.update({
        "n_cells": int(valid.sum()),
        "total_counts": float(totals[valid].sum()),
        "mean_library_size": float(totals[valid].mean()) if valid.any() else np.nan,
        "genes": all_genes,
        # Gene sums support donor pseudobulk and gene-level summaries.
        "gene_sums": selected[:, valid].sum(axis=1),
        # Per-cell log2(CPM + 1) matrix kept only as per-gene cell means for leave-one-gene-out.
        "gene_cell_mean_logcpm": np.log2(selected[:, valid] / totals[None, valid] * 1e6 + 1.0).mean(axis=1)
        if valid.any() else np.full(len(all_genes), np.nan),
    })
    return record


def sample_table(records: list[dict[str, Any]]) -> pd.DataFrame:
    """Long table: one row per sample x gene with cell-mean log2CPM and count sums."""
    rows = []
    for rec in records:
        if rec.get("n_cells", 0) < MIN_CELLS_PER_SAMPLE:
            continue
        for i, gene in enumerate(rec["genes"]):
            rows.append({
                **{k: rec[k] for k in ["sample_folder", "Sample_ID", "Donor", "Region", "Sex", "Condition", "Group",
                                       "C9_pos", "n_cells", "total_counts", "mean_library_size"]},
                "gene_id": gene,
                "cell_mean_logcpm": float(rec["gene_cell_mean_logcpm"][i]),
                "gene_count_sum": float(rec["gene_sums"][i]),
            })
    return pd.DataFrame(rows)


def module_scores(long: pd.DataFrame, genes: list[str], region_map: dict[str, str]) -> pd.DataFrame:
    """Donor x scope module scores (cell-weighted mean score and pseudobulk score)."""
    sub = long[long["gene_id"].isin(genes)].copy()
    sub["scope"] = sub["Region"].map(region_map)
    sub = sub.dropna(subset=["scope"])
    # Mean over genes of per-cell log2CPM equals per-cell module score averaged over cells.
    per_sample = sub.groupby(["sample_folder", "Donor", "scope"], observed=True).agg(
        score=("cell_mean_logcpm", "mean"), n_cells=("n_cells", "first"),
        total_counts=("total_counts", "first"), mean_library_size=("mean_library_size", "first"),
    ).reset_index()
    per_sample["w_score"] = per_sample["score"] * per_sample["n_cells"]
    per_sample["w_lib"] = per_sample["mean_library_size"] * per_sample["n_cells"]
    donor = per_sample.groupby(["Donor", "scope"], observed=True).agg(
        w_score=("w_score", "sum"), w_lib=("w_lib", "sum"), n_cells=("n_cells", "sum"),
        total_counts=("total_counts", "sum"), n_samples=("sample_folder", "nunique"),
    ).reset_index()
    donor["module_score"] = donor["w_score"] / donor["n_cells"]
    donor["mean_library_size"] = donor["w_lib"] / donor["n_cells"]

    gene_sums = sub.groupby(["Donor", "scope", "gene_id"], observed=True)["gene_count_sum"].sum().reset_index()
    gene_sums = gene_sums.merge(donor[["Donor", "scope", "total_counts"]], on=["Donor", "scope"])
    gene_sums["pb_logcpm"] = np.log2(gene_sums["gene_count_sum"] / gene_sums["total_counts"] * 1e6 + 1.0)
    pb = gene_sums.groupby(["Donor", "scope"], observed=True)["pb_logcpm"].mean().rename("pseudobulk_score")
    return donor.drop(columns=["w_score", "w_lib"]).merge(pb.reset_index(), on=["Donor", "scope"])


def contrast_frame(df: pd.DataFrame, contrast: str) -> pd.DataFrame:
    case_col, case = ("Group", contrast.split("_vs_")[0]) if contrast.startswith(("SALS", "C9ALS")) else ("Condition", contrast.split("_vs_")[0])
    out = df[df[case_col].eq(case) | df["Condition"].eq("PN")].copy()
    out["case"] = out[case_col].eq(case).astype(float)
    return out


def fit(df: pd.DataFrame, response: str, covariates: list[str]) -> dict[str, Any]:
    data = df.dropna(subset=[response, "case", *covariates]).copy()
    x = pd.DataFrame({"case": data["case"]}, index=data.index)
    for cov in covariates:
        if cov == "Sex":
            x["sex_M"] = data["Sex"].astype(str).eq("M").astype(float)
        else:
            x[cov] = data[cov].astype(float)
    x = sm.add_constant(x, has_constant="add")
    n_case, n_ctrl = int(data["case"].sum()), int((1 - data["case"]).sum())
    base = {"n_case_donors": n_case, "n_control_donors": n_ctrl, "covariates": "+".join(covariates),
            "beta": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan}
    if n_case < 5 or n_ctrl < 5 or np.linalg.matrix_rank(x.to_numpy()) < x.shape[1]:
        return {**base, "status": "skipped"}
    res = sm.OLS(data[response].astype(float), x).fit()
    lo, hi = res.conf_int().loc["case"]
    return {**base, "status": "completed", "beta": res.params["case"], "ci_low": lo, "ci_high": hi,
            "p_value": res.pvalues["case"]}


def run_models(scores: dict[str, pd.DataFrame], clinical: pd.DataFrame) -> pd.DataFrame:
    specs = [
        # analysis_id, scope, contrast, response, covariates, tier
        ("primary_BA4", "BA4", "ALS_vs_PN", "module_score", ["Sex"], "primary"),
        ("BA4_age_PMI", "BA4", "ALS_vs_PN", "module_score", ["Sex", "age_death", "PMI_h"], "sensitivity"),
        ("BA4_library_depth", "BA4", "ALS_vs_PN", "module_score", ["Sex", "log10_library_size"], "sensitivity"),
        ("BA4_pseudobulk", "BA4", "ALS_vs_PN", "pseudobulk_score", ["Sex"], "sensitivity"),
        ("BA4_SALS", "BA4", "SALS_vs_PN", "module_score", ["Sex"], "subgroup"),
        ("BA4_C9ALS", "BA4", "C9ALS_vs_PN", "module_score", ["Sex"], "subgroup"),
        ("BA4_BA9_mean", "BA4_BA9", "ALS_vs_PN", "module_score", ["Sex"], "sensitivity"),
        ("BA9_only", "BA9", "ALS_vs_PN", "module_score", ["Sex"], "sensitivity"),
        ("BA4_FTLD_context", "BA4", "FTLD_vs_PN", "module_score", ["Sex"], "context"),
    ]
    rows = []
    for direction, df in scores.items():
        data = df.merge(clinical, on="Donor", how="left")
        data["log10_library_size"] = np.log10(data["mean_library_size"])
        for analysis_id, scope, contrast, response, covs, tier in specs:
            sub = contrast_frame(data[data["scope"].eq(scope)], contrast)
            rows.append({"analysis_id": analysis_id, "tier": tier, "module": direction, "scope": scope,
                         "contrast": contrast, "response": response, **fit(sub, response, covs)})
    res = pd.DataFrame(rows)
    res["BH_within_analysis"] = np.nan
    for _, idx in res[res["status"].eq("completed")].groupby("analysis_id").groups.items():
        res.loc[idx, "BH_within_analysis"] = multipletests(res.loc[idx, "p_value"], method="fdr_bh")[1]
    res["direction_consistent"] = np.where(res["module"].eq("up"), res["beta"] > 0, res["beta"] < 0)
    return res


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=None, help="Process only the first N samples (smoke test).")
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--reuse", action="store_true", help="Reuse w50_sample_gene_table.tsv.gz.")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    modules = load_modules()
    module_genes = {d: m["gene_id"].tolist() for d, m in modules.items()}
    pd.concat([m.assign(direction=d) for d, m in modules.items()]).to_csv(OUT / "w50_module_registry.tsv", sep="\t", index=False)

    long_path = OUT / "w50_sample_gene_table.tsv.gz"
    if args.reuse and long_path.exists():
        long = pd.read_csv(long_path, sep="\t")
    else:
        manifest = pd.read_csv(MANIFEST, sep="\t")
        manifest = manifest[manifest["counts_fil_mtx_exists"].astype(str).eq("True")]
        if args.limit:
            manifest = manifest.head(args.limit)
        jobs = [(str(resolve_source_path(p)), module_genes) for p in manifest["path"]]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            records = []
            for i, rec in enumerate(pool.map(project_sample, jobs), 1):
                records.append(rec)
                print(f"[W50] {i}/{len(jobs)} {rec['sample_folder']} astro_cells={rec.get('n_cells', 0)}", flush=True)
        pd.DataFrame([{k: v for k, v in r.items() if k not in {"genes", "gene_sums", "gene_cell_mean_logcpm"}}
                      for r in records]).to_csv(OUT / "w50_sample_coverage.tsv", sep="\t", index=False)
        long = sample_table(records)
        long.to_csv(long_path, sep="\t", index=False)

    meta = long.drop_duplicates("Donor").set_index("Donor")[["Sex", "Condition", "Group", "C9_pos"]]
    region_maps = {
        "BA4": {"BA4": "BA4"}, "BA9": {"BA9": "BA9"},
        "BA4_BA9": {"BA4": "BA4_BA9", "BA9": "BA4_BA9"},
    }
    scores = {}
    for direction, genes in module_genes.items():
        parts = [module_scores(long, genes, rmap) for rmap in region_maps.values()]
        scores[direction] = pd.concat(parts).merge(meta, left_on="Donor", right_index=True)
    pd.concat([s.assign(module=d) for d, s in scores.items()]).to_csv(OUT / "w50_donor_module_scores.tsv", sep="\t", index=False)

    clinical = pd.read_csv(CLINICAL, sep="\t").drop_duplicates("Donor")[["Donor", "PMI_h", "Age of death (Y)"]]
    clinical = clinical.rename(columns={"Age of death (Y)": "age_death"})
    results = run_models(scores, clinical)
    results.to_csv(OUT / "w50_validation_results.tsv", sep="\t", index=False)

    # Leave-one-gene-out and gene-level directions for the up module in BA4.
    loo_rows, gene_rows = [], []
    ba4 = long[long["Region"].eq("BA4")]
    for gene in module_genes["up"]:
        s = module_scores(ba4, [g for g in module_genes["up"] if g != gene], {"BA4": "BA4"}).merge(meta, left_on="Donor", right_index=True)
        loo_rows.append({"left_out_gene": gene, **fit(contrast_frame(s, "ALS_vs_PN"), "module_score", ["Sex"])})
        g = module_scores(ba4, [gene], {"BA4": "BA4"}).merge(meta, left_on="Donor", right_index=True)
        gene_rows.append({"gene_id": gene, **fit(contrast_frame(g, "ALS_vs_PN"), "module_score", ["Sex"])})
    names = modules["up"].set_index("gene_id")["common_name"]
    loo = pd.DataFrame(loo_rows)
    loo.insert(1, "gene", loo["left_out_gene"].map(names))
    loo.to_csv(OUT / "w50_up_module_leave_one_gene_out_BA4.tsv", sep="\t", index=False)
    genes_df = pd.DataFrame(gene_rows)
    genes_df.insert(1, "gene", genes_df["gene_id"].map(names))
    genes_df["BH"] = multipletests(genes_df["p_value"].fillna(1), method="fdr_bh")[1]
    genes_df.to_csv(OUT / "w50_up_gene_level_BA4.tsv", sep="\t", index=False)
    print(results.to_string(), flush=True)


if __name__ == "__main__":
    main()
