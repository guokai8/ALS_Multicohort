#!/usr/bin/env python3
"""Audit Pineda S01 sample information for P01 metadata use.

The Pineda S01 workbook can enrich the Pineda ALS/FTLD cortical context with
donor-level clinical/pathology variables. It must not be used to fill
GSE330130/Takeuchi metadata because the ID systems are independent; numeric
Pineda donor IDs such as 101 are not GSE330130 MO101 libraries.

Inputs may be supplied by environment variable:
- PINEDA_S01_SAMPLE_INFORMATION
- GSE330130_SERIES_MATRIX
"""

from __future__ import annotations

import hashlib
import os
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

from p01_paths import plan_root, project_root, source_data_root


ROOT = project_root()
SOURCE_DATA_ROOT = source_data_root()
PLAN_ROOT = plan_root()
OUT = ROOT / "results" / "Supplemental_QC_covariate_audit"
DATE_TAG = "2026-08-12"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_input(env_name: str, candidates: Iterable[Path]) -> tuple[Path, str]:
    value = os.environ.get(env_name)
    if value:
        return Path(value).expanduser().resolve(), f"{{{env_name}}}"
    for candidate in candidates:
        if candidate.exists():
            resolved = candidate.resolve()
            if resolved.is_relative_to(ROOT):
                return resolved, str(resolved.relative_to(ROOT))
            if resolved.is_relative_to(SOURCE_DATA_ROOT):
                return resolved, "{P01_SOURCE_DATA_ROOT}/" + str(resolved.relative_to(SOURCE_DATA_ROOT))
            return resolved, f"{{{env_name}}}"
    raise FileNotFoundError(
        f"Set {env_name}; checked: " + "; ".join(str(candidate) for candidate in candidates)
    )


def normalize_id(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    try:
        number = float(text)
        if number.is_integer():
            return str(int(number))
    except ValueError:
        pass
    return text


def normalize_numeric(value: object) -> float | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"na", "nan", "unknown", "not recorded", "not applicable"}:
        return None
    match = re.search(r"-?\d+(?:\.\d+)?", text)
    if not match:
        return None
    return float(match.group(0))


def parse_gse330130_ids(path: Path) -> pd.DataFrame:
    fields: dict[str, list[list[str]]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("!Sample_"):
            continue
        parts = line.rstrip("\n").split("\t")
        fields.setdefault(parts[0], []).append([value.strip().strip('"') for value in parts[1:]])
    titles = fields.get("!Sample_title", [[]])[0]
    gsms = fields.get("!Sample_geo_accession", [[]])[0]
    chars = fields.get("!Sample_characteristics_ch1", [])
    rows: list[dict[str, object]] = []
    for i, title in enumerate(titles):
        rec: dict[str, object] = {
            "gsm": gsms[i] if i < len(gsms) else "",
            "sample_title": title,
            "mo_id": "",
            "mo_numeric": "",
            "cgnd_id": "",
            "u54_id": "",
        }
        mo = re.search(r"(MO\d+)", title)
        cgnd = re.search(r"(CGND-HRA-\d+(?:-\d+)?)", title)
        u54 = re.search(r"(U54-HRA-\d+(?:-A\d+)?)", title)
        if mo:
            rec["mo_id"] = mo.group(1)
            rec["mo_numeric"] = mo.group(1).replace("MO", "").lstrip("0")
        if cgnd:
            rec["cgnd_id"] = cgnd.group(1)
        if u54:
            rec["u54_id"] = u54.group(1)
        for values in chars:
            if i >= len(values) or ":" not in values[i]:
                continue
            key, raw = values[i].split(":", 1)
            rec[key.strip().lower().replace(" ", "_")] = raw.strip()
        rows.append(rec)
    return pd.DataFrame(rows)


def field_coverage(donor: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for col in donor.columns:
        values = donor[col]
        rows.append(
            {
                "field": col,
                "non_missing_rows": int(values.notna().sum()),
                "total_rows": len(values),
                "coverage_fraction": round(float(values.notna().mean()), 4),
                "n_unique_non_missing": int(values.nunique(dropna=True)),
                "p01_use": field_use(col),
            }
        )
    return pd.DataFrame(rows)


def field_use(field: str) -> str:
    if field in {"Donor", "Group", "Sex"}:
        return "pineda_context_core_metadata"
    if field in {"PMI_h", "Race", "Age of death (Y)"}:
        return "pineda_context_covariate_or_demographic"
    if field in {"Braak_stage", "Thal_phase"}:
        return "pineda_context_pathology_covariate_only_not_gse330130_pathology"
    if field in {"Age of onset cognitive (Y)", "Age of onset motor (Y)", "Disease duration (Y)"}:
        return "pineda_context_clinical_variable_only_not_gse330130_clinical_endpoint"
    return "pineda_context_other"


def summarize_s01(
    s01_path: Path,
    s01_label: str,
    donor: pd.DataFrame,
    sample_mapping: pd.DataFrame,
    manifest: pd.DataFrame,
    gse330130: pd.DataFrame,
) -> pd.DataFrame:
    sample_mapping = sample_mapping.copy()
    manifest = manifest.copy()
    donor = donor.copy()
    sample_mapping["Sample_ID_norm"] = sample_mapping["Sample_ID"].map(normalize_id)
    sample_mapping["Donor_norm"] = sample_mapping["Donor"].map(normalize_id)
    donor["Donor_norm"] = donor["Donor"].map(normalize_id)
    manifest["sample_folder_norm"] = manifest["sample_folder"].map(normalize_id)
    manifest["Donor_norm"] = manifest["Donor"].map(normalize_id)

    s01_samples = set(sample_mapping["Sample_ID_norm"])
    manifest_samples = set(manifest["sample_folder_norm"])
    s01_donors = {x for x in sample_mapping["Donor_norm"] if x}
    manifest_donors = {x for x in manifest["Donor_norm"] if x}

    gse_mo_raw = {x for x in gse330130["mo_id"].astype(str) if x}
    gse_mo_numeric = {x for x in gse330130["mo_numeric"].astype(str) if x}
    gse_cgnd = {x for x in gse330130["cgnd_id"].astype(str) if x}
    gse_u54 = {x for x in gse330130["u54_id"].astype(str) if x}
    gse_titles = set(gse330130["sample_title"].astype(str))

    rows = [
        ("s01_source_label", s01_label, "Release-safe source label for S01 workbook."),
        ("s01_source_sha256", sha256(s01_path), "Input checksum."),
        ("s01_donor_information_rows", len(donor), "Rows in S01 donor_information sheet, including note rows."),
        ("s01_sample_mapping_rows", len(sample_mapping), "Rows in S01 sample_mapping sheet."),
        ("s01_donor_ids_nonmissing", int(donor["Donor_norm"].astype(bool).sum()), "Non-empty donor IDs in donor_information."),
        ("s01_sample_mapping_unique_donors", sample_mapping["Donor_norm"].nunique(), "Unique donors represented in sample_mapping."),
        ("s01_sample_mapping_region_counts", counts(sample_mapping["Region"]), "BA4/BA9/M1 sample coverage in S01."),
        ("s01_group_counts", counts(donor["Group"]), "Disease/group rows in donor_information."),
        ("manifest_rows", len(manifest), "P01 Pineda manifest rows."),
        ("manifest_unique_donors", manifest["Donor_norm"].nunique(), "P01 Pineda manifest unique donors."),
        (
            "s01_samples_to_p01_pineda_manifest",
            f"{len(s01_samples & manifest_samples)}/{len(s01_samples)}/{len(manifest_samples)}",
            "S01 Sample_ID overlap with P01 pineda_sample_manifest sample_folder.",
        ),
        (
            "s01_donors_to_p01_pineda_manifest",
            f"{len(s01_donors & manifest_donors)}/{len(s01_donors)}/{len(manifest_donors)}",
            "S01 donor overlap with P01 pineda_sample_manifest Donor.",
        ),
        (
            "s01_donor_ids_to_gse330130_mo_raw",
            len(s01_donors & gse_mo_raw),
            "Exact Pineda donor ID to GSE330130 MOxxx ID overlap; expected zero.",
        ),
        (
            "s01_sample_ids_to_gse330130_titles",
            len(s01_samples & gse_titles),
            "Exact Pineda Sample_ID to GSE330130 sample-title overlap; expected zero.",
        ),
        (
            "s01_donor_ids_to_gse330130_cgnd",
            len(s01_donors & gse_cgnd),
            "Exact Pineda donor ID to GSE330130 CGND-HRA ID overlap; expected zero.",
        ),
        (
            "s01_donor_ids_to_gse330130_u54",
            len(s01_donors & gse_u54),
            "Exact Pineda donor ID to GSE330130 U54-HRA ID overlap; expected zero.",
        ),
        (
            "false_positive_numeric_overlap_with_gse330130_mo_numbers",
            len(s01_donors & gse_mo_numeric),
            "Numeric-only overlap caused by Pineda donor numbers and GSE330130 MO numbers; do not interpret as donor identity.",
        ),
        (
            "p01_use_decision",
            "pineda_context_metadata_only_not_gse330130_metadata",
            "S01 can enrich Pineda context but cannot fill GSE330130/Takeuchi metadata.",
        ),
    ]
    return pd.DataFrame(rows, columns=["metric", "value", "interpretation"])


def counts(values: pd.Series) -> str:
    clean = values.fillna("NA").astype(str)
    return ";".join(f"{key}:{value}" for key, value in clean.value_counts(dropna=False).items())


def build_manifest_join(donor: pd.DataFrame, sample_mapping: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    donor = donor.copy()
    sample_mapping = sample_mapping.copy()
    manifest = manifest.copy()
    donor["Donor_norm"] = donor["Donor"].map(normalize_id)
    sample_mapping["Sample_ID_norm"] = sample_mapping["Sample_ID"].map(normalize_id)
    sample_mapping["Donor_norm"] = sample_mapping["Donor"].map(normalize_id)
    manifest["sample_folder_norm"] = manifest["sample_folder"].map(normalize_id)
    manifest["Donor_norm"] = manifest["Donor"].map(normalize_id)

    sample_mapping_for_join = sample_mapping[["Sample_ID_norm", "Donor_norm", "Region"]].rename(
        columns={"Donor_norm": "Donor_norm_s01_mapping"}
    )
    joined = manifest.merge(
        sample_mapping_for_join,
        left_on="sample_folder_norm",
        right_on="Sample_ID_norm",
        how="left",
        suffixes=("_manifest", "_s01_mapping"),
    ).merge(
        donor[
            [
                "Donor_norm",
                "Group",
                "Sex",
                "Race",
                "PMI_h",
                "Braak_stage",
                "Thal_phase",
                "Age of onset cognitive (Y)",
                "Age of onset motor (Y)",
                "Age of death (Y)",
                "Disease duration (Y)",
            ]
        ],
        on="Donor_norm",
        how="left",
        suffixes=("", "_s01_donor"),
    )
    joined["sample_mapping_match"] = joined["Sample_ID_norm"].notna()
    joined["donor_match"] = joined["Group"].notna()
    joined["donor_mismatch"] = joined["Donor_norm"] != joined["Donor_norm_s01_mapping"]
    joined["region_mismatch"] = joined["Region_manifest"].astype(str) != joined["Region_s01_mapping"].astype(str)
    keep = [
        "sample_folder",
        "Donor",
        "Donor_norm",
        "Region_manifest",
        "Region_s01_mapping",
        "Condition",
        "Group_s01_donor" if "Group_s01_donor" in joined.columns else "Group",
        "Sex",
        "Race",
        "PMI_h",
        "Braak_stage",
        "Thal_phase",
        "Age of onset cognitive (Y)",
        "Age of onset motor (Y)",
        "Age of death (Y)",
        "Disease duration (Y)",
        "n_cells",
        "sample_mapping_match",
        "donor_match",
        "donor_mismatch",
        "region_mismatch",
    ]
    keep = [col for col in keep if col in joined.columns]
    return joined[keep].copy()


def build_gse330130_boundary(donor: pd.DataFrame, sample_mapping: pd.DataFrame, gse330130: pd.DataFrame) -> pd.DataFrame:
    donor_ids = {normalize_id(x) for x in donor["Donor"] if normalize_id(x)}
    sample_ids = {normalize_id(x) for x in sample_mapping["Sample_ID"] if normalize_id(x)}
    comparisons = [
        ("pineda_donor_vs_gse330130_mo_raw", donor_ids, set(gse330130["mo_id"].astype(str))),
        ("pineda_donor_vs_gse330130_mo_numeric", donor_ids, set(gse330130["mo_numeric"].astype(str))),
        ("pineda_donor_vs_gse330130_cgnd_id", donor_ids, set(gse330130["cgnd_id"].astype(str))),
        ("pineda_donor_vs_gse330130_u54_id", donor_ids, set(gse330130["u54_id"].astype(str))),
        ("pineda_sample_vs_gse330130_title", sample_ids, set(gse330130["sample_title"].astype(str))),
    ]
    rows: list[dict[str, object]] = []
    for comparison, left, right in comparisons:
        overlap = sorted({x for x in left & right if x})
        rows.append(
            {
                "comparison": comparison,
                "n_overlap": len(overlap),
                "overlap_examples": ";".join(overlap[:30]) if overlap else "none",
                "interpretation": "false_positive_numeric_overlap_do_not_use"
                if comparison == "pineda_donor_vs_gse330130_mo_numeric" and overlap
                else ("no_exact_id_overlap" if not overlap else "manual_review_required"),
            }
        )
    return pd.DataFrame(rows)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.copy().where(pd.notna(df), "NA").to_csv(path, sep="\t", index=False)


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "No rows."
    safe = df.copy().where(pd.notna(df), "NA").astype(str)
    lines = [
        "| " + " | ".join(safe.columns) + " |",
        "| " + " | ".join(["---"] * len(safe.columns)) + " |",
    ]
    for _, row in safe.iterrows():
        lines.append("| " + " | ".join(row[col].replace("|", "/") for col in safe.columns) + " |")
    return "\n".join(lines)


def write_report(summary: pd.DataFrame, coverage: pd.DataFrame, boundary: pd.DataFrame) -> Path:
    path = PLAN_ROOT / f"pineda_s01_metadata_audit_{DATE_TAG}.md"
    key_metrics = summary[
        summary["metric"].isin(
            [
                "s01_donor_information_rows",
                "s01_sample_mapping_rows",
                "s01_samples_to_p01_pineda_manifest",
                "s01_donors_to_p01_pineda_manifest",
                "s01_donor_ids_to_gse330130_mo_raw",
                "s01_donor_ids_to_gse330130_cgnd",
                "false_positive_numeric_overlap_with_gse330130_mo_numbers",
                "p01_use_decision",
            ]
        )
    ].copy()
    coverage_preview = coverage[
        coverage["field"].isin(
            [
                "Group",
                "Sex",
                "PMI_h",
                "Braak_stage",
                "Thal_phase",
                "Age of onset cognitive (Y)",
                "Age of onset motor (Y)",
                "Age of death (Y)",
                "Disease duration (Y)",
            ]
        )
    ].copy()
    text = f"""# Pineda S01 sample-information metadata audit ({DATE_TAG})

## Bottom Line

`S01_sample_information.xlsx` is useful for Pineda-specific donor/sample metadata. It fully maps to the P01 Pineda sample manifest by `Sample_ID` and donor ID, and adds donor-level clinical/pathology context such as PMI, Braak stage, Thal phase, cognitive/motor onset age, age at death and disease duration.

It cannot supplement GSE330130/Takeuchi metadata. Exact overlap with GSE330130 `MOxxx`, `CGND-HRA-*`, `U54-HRA-*` and sample-title IDs is zero. Numeric overlap between Pineda donor IDs such as `101` and GSE330130 library numbers such as `MO101` is a false-positive identifier collision and must not be interpreted as donor identity.

## Key Metrics

{markdown_table(key_metrics)}

## Clinically Useful Pineda Fields

{markdown_table(coverage_preview)}

## GSE330130 Boundary Check

{markdown_table(boundary)}

## P01 Use Decision

Allowed now:

- Use S01 to strengthen Pineda ALS/FTLD cortical-context tables and Methods.
- Use donor-level S01 fields as Pineda context covariates or descriptive metadata, with clear labels that they belong to the Pineda cohort.
- Use Braak/Thal only as Pineda context/pathology covariates, not as GSE330130 pathology.

Forbidden now:

- Do not fill GSE330130 age, onset, duration, pathology, cognitive or genotype fields from Pineda S01.
- Do not equate Pineda donor `101` with GSE330130 `MO101`; the matching numeric strings are not shared biological IDs.
- Do not upgrade Brain-route clinical/cognitive/pathology claims from this audit alone.

## Machine-Readable Outputs

- `results/Supplemental_QC_covariate_audit/pineda_s01_metadata_summary_{DATE_TAG}.tsv`
- `results/Supplemental_QC_covariate_audit/pineda_s01_field_coverage_{DATE_TAG}.tsv`
- `results/Supplemental_QC_covariate_audit/pineda_s01_manifest_join_{DATE_TAG}.tsv`
- `results/Supplemental_QC_covariate_audit/pineda_s01_gse330130_overlap_boundary_{DATE_TAG}.tsv`
"""
    path.write_text(text, encoding="utf-8")
    return path


def main() -> None:
    s01_path, s01_label = resolve_input(
        "PINEDA_S01_SAMPLE_INFORMATION",
        [
            SOURCE_DATA_ROOT
            / "Other_Datasets"
            / "Pineda_snRNAseq"
            / "files"
            / "pineda"
            / "S01_sample_information.xlsx",
        ],
    )
    gse_path, _ = resolve_input(
        "GSE330130_SERIES_MATRIX",
        [
            SOURCE_DATA_ROOT / "GSE330130" / "GSE330130_series_matrix",
            ROOT / "metadata" / "GSE330130_series_matrix",
        ],
    )
    manifest_path = ROOT / "metadata" / "pineda_sample_manifest.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError(manifest_path)

    donor = pd.read_excel(s01_path, sheet_name="donor_information")
    sample_mapping = pd.read_excel(s01_path, sheet_name="sample_mapping")
    manifest = pd.read_csv(manifest_path, sep="\t")
    gse330130 = parse_gse330130_ids(gse_path)

    summary = summarize_s01(s01_path, s01_label, donor, sample_mapping, manifest, gse330130)
    coverage = field_coverage(donor)
    manifest_join = build_manifest_join(donor, sample_mapping, manifest)
    boundary = build_gse330130_boundary(donor, sample_mapping, gse330130)

    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv(summary, OUT / f"pineda_s01_metadata_summary_{DATE_TAG}.tsv")
    write_tsv(coverage, OUT / f"pineda_s01_field_coverage_{DATE_TAG}.tsv")
    write_tsv(manifest_join, OUT / f"pineda_s01_manifest_join_{DATE_TAG}.tsv")
    write_tsv(boundary, OUT / f"pineda_s01_gse330130_overlap_boundary_{DATE_TAG}.tsv")
    report_path = write_report(summary, coverage, boundary)

    print(f"summary_rows={len(summary)}")
    print(f"field_coverage_rows={len(coverage)}")
    print(f"manifest_join_rows={len(manifest_join)}")
    print(f"gse330130_boundary_rows={len(boundary)}")
    print(f"report={report_path}")


if __name__ == "__main__":
    main()
