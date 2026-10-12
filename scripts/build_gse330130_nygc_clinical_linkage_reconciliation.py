#!/usr/bin/env python3
"""Reconcile GSE330130 GEO metadata with a NYGC clinical linkage table.

This is an older compatibility audit for the NYGC clinical linkage CSV. Current
P01 clinical/diagnosis interpretation follows the latest Synapse/NYGC bridge:
for direct matches, Synapse/NYGC curated metadata is the resolved standard and
GEO fields are retained as raw sequencing provenance plus conflict notes.

Inputs are intentionally provided by environment variable so local absolute
paths do not enter the release package:
- GSE330130_SERIES_MATRIX
- NYGC_CLINICAL_LINKAGE_CSV
"""

from __future__ import annotations

import hashlib
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd

from p01_paths import plan_root, project_root, source_data_root


ROOT = project_root()
PLAN_ROOT = plan_root()
OUT = ROOT / "results" / "Supplemental_QC_covariate_audit"
DATE_TAG = "2026-08-12"

SERIES_ENV = "GSE330130_SERIES_MATRIX"
NYGC_ENV = "NYGC_CLINICAL_LINKAGE_CSV"
SERIES_DURATION_RAW = "series_disease_duration_(onset_to_tracheostomy_or_death)_in_months"
SERIES_DURATION_ALIAS = "series_disease_duration_onset_to_tracheostomy_or_death_months"


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
            return candidate.resolve(), str(candidate.relative_to(ROOT)) if candidate.is_relative_to(ROOT) else f"{{{env_name}}}"
    raise FileNotFoundError(
        f"Set {env_name} to the local source file. Checked candidates: "
        + "; ".join(str(x) for x in candidates)
    )


def read_series_sample_fields(path: Path) -> dict[str, list[list[str]]]:
    fields: dict[str, list[list[str]]] = defaultdict(list)
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("!Sample_"):
            continue
        parts = line.rstrip("\n").split("\t")
        fields[parts[0]].append([value.strip().strip('"') for value in parts[1:]])
    return dict(fields)


def parse_series_libraries(path: Path, source_label: str) -> pd.DataFrame:
    fields = read_series_sample_fields(path)
    titles = fields["!Sample_title"][0]
    gsms = fields["!Sample_geo_accession"][0]
    sources = fields["!Sample_source_name_ch1"][0]
    characteristics = fields.get("!Sample_characteristics_ch1", [])

    sample_rows: list[dict[str, object]] = []
    for i, title in enumerate(titles):
        rec: dict[str, object] = {
            "series_source_label": source_label,
            "series_sha256": sha256(path),
            "gsm": gsms[i],
            "sample_title": title,
            "sample_source_name": sources[i],
        }
        mo = re.search(r"(MO\d+)", title)
        if mo:
            rec["mo_id"] = mo.group(1)
        elif "U54-HRA-135" in title:
            rec["mo_id"] = "MO023"
        else:
            rec["mo_id"] = ""

        rec["assay_from_title"] = "ATAC" if title.endswith("-ATAC") else ("GEX" if title.endswith("-GEX") else "OTHER")

        cgnd = re.search(r"(CGND-HRA-\d+)", title)
        u54 = re.search(r"(U54-HRA-\d+(?:-A\d+)?)", title)
        rec["series_cgnd_id"] = cgnd.group(1) if cgnd else ""
        rec["series_u54_id"] = u54.group(1) if u54 else ""

        for values in characteristics:
            value = values[i].strip().strip('"') if i < len(values) else ""
            if not value or ":" not in value:
                continue
            key, raw = value.split(":", 1)
            rec[f"series_{key.strip().lower().replace(' ', '_')}"] = raw.strip()
        sample_rows.append(rec)

    collapsed: list[dict[str, object]] = []
    for mo_id, sub in pd.DataFrame(sample_rows).groupby("mo_id", dropna=False):
        row: dict[str, object] = {
            "mo_id": mo_id,
            "n_gsm_columns": len(sub),
            "series_assays": ";".join(sorted(set(sub["assay_from_title"].astype(str)))),
            "series_gsms": ";".join(sub["gsm"].astype(str)),
            "series_titles": ";".join(sub["sample_title"].astype(str)),
        }
        for col in [
            "series_source_label",
            "series_sha256",
            "series_cgnd_id",
            "series_u54_id",
            "series_tissue",
            "series_genotype",
            "series_sex",
            "series_pmi",
            "series_age_at_death",
            "series_age_at_symptom_onset",
            "series_site_of_motor_onset",
            SERIES_DURATION_RAW,
        ]:
            if col not in sub.columns:
                row[col] = ""
                continue
            values = sorted({str(v).strip() for v in sub[col].dropna() if str(v).strip()})
            row[col] = ";".join(values)
        row[SERIES_DURATION_ALIAS] = row.pop(SERIES_DURATION_RAW, "")
        row["series_broad_region"] = broad_region(row.get("series_tissue", ""))
        row["series_duration_years_from_age"] = duration_years(row.get("series_age_at_death"), row.get("series_age_at_symptom_onset"))
        collapsed.append(row)
    return pd.DataFrame(collapsed).sort_values("mo_id").reset_index(drop=True)


def duration_years(age_death: object, age_onset: object) -> object:
    try:
        return float(age_death) - float(age_onset)
    except Exception:
        return ""


def broad_region(value: object) -> str:
    text = str(value)
    if "Spinal" in text or "spinal" in text:
        return "LSC" if "Lumbar" in text or "lumbar" in text else "spinal_cord"
    if "Cortex_Motor" in text or "Motor_Cortex" in text or "Motor" in text:
        return "MTC"
    if "Cortex" in text:
        return "cortex_other"
    return "unknown"


def normalize_disease(value: object) -> str:
    text = str(value).strip().lower().replace(" ", "_")
    if text == "als":
        return "ALS"
    if text == "control":
        return "Control"
    if text in {"als-ftd", "als_ftd"}:
        return "ALS_FTD"
    if text:
        return text
    return "missing"


def numeric_equal(a: object, b: object) -> bool:
    try:
        return abs(float(a) - float(b)) < 1e-6
    except Exception:
        return str(a).strip() == str(b).strip() or not str(a).strip() or not str(b).strip()


def read_nygc(path: Path, source_label: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["nygc_source_label"] = source_label
    df["nygc_sha256"] = sha256(path)
    return df


def attach_nygc(series_df: pd.DataFrame, nygc: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "rna_id",
        "donor_id",
        "disease_group_rna",
        "tissue_clean",
        "region",
        "rin",
        "pmi",
        "used_in_analysis",
        "sra_run",
        "sex_call",
        "age_at_symptom_onset",
        "age_at_death",
        "disease_duration_in_months",
        "site_of_motor_onset",
        "family_history_of_als_ftd",
        "c9orf72_genotype",
        "atxn2_genotype",
        "in_nygc_count_matrix",
        "in_nygc_splice_matrix",
        "usable_als_spinal_duration_sample",
        "usable_als_spinal_duration_primary",
        "nygc_source_label",
        "nygc_sha256",
    ]
    missing = [col for col in keep if col not in nygc.columns]
    if missing:
        raise ValueError(f"NYGC linkage table is missing required columns: {missing}")
    nygc_small = nygc[keep].copy().add_prefix("nygc_")
    merged = series_df.merge(
        nygc_small,
        left_on="series_cgnd_id",
        right_on="nygc_rna_id",
        how="left",
    )
    merged["nygc_direct_match"] = merged["nygc_rna_id"].notna().map({True: "yes", False: "no"})
    return merged


def classify_rows(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, raw in df.iterrows():
        row = raw.to_dict()
        conflict_fields: list[str] = []
        warning_fields: list[str] = []

        if row["nygc_direct_match"] == "yes":
            if normalize_disease(row.get("series_genotype")) != normalize_disease(row.get("nygc_disease_group_rna")):
                conflict_fields.append("disease_label")

            series_region = broad_region(row.get("series_tissue"))
            nygc_region = broad_region(row.get("nygc_tissue_clean"))
            if series_region != nygc_region:
                conflict_fields.append("broad_region")
            elif str(row.get("series_tissue", "")).strip() != str(row.get("nygc_tissue_clean", "")).strip():
                warning_fields.append("fine_tissue_label")

            for series_col, nygc_col, label in [
                ("series_sex", "nygc_sex_call", "sex"),
                ("series_pmi", "nygc_pmi", "pmi"),
                ("series_age_at_death", "nygc_age_at_death", "age_at_death"),
                ("series_age_at_symptom_onset", "nygc_age_at_symptom_onset", "age_at_symptom_onset"),
            ]:
                if not numeric_equal(row.get(series_col, ""), row.get(nygc_col, "")):
                    conflict_fields.append(label)

            if str(row.get("series_site_of_motor_onset", "")).strip() and str(row.get("nygc_site_of_motor_onset", "")).strip():
                if str(row.get("series_site_of_motor_onset")).strip() != str(row.get("nygc_site_of_motor_onset")).strip():
                    conflict_fields.append("site_of_motor_onset")

        row["reconciliation_conflict_fields"] = ";".join(conflict_fields) if conflict_fields else "none"
        row["reconciliation_warning_fields"] = ";".join(warning_fields) if warning_fields else "none"
        row["resolved_source_policy"] = (
            "geo_series_matrix_no_direct_synapse_nygc_match"
            if row["nygc_direct_match"] == "no"
            else "synapse_nygc_curated_priority_for_direct_match"
        )
        note_parts = []
        if row["nygc_direct_match"] == "no":
            note_parts.append("No direct Synapse/NYGC RNA-ID match; GEO series metadata remains the only resolved local context.")
        else:
            note_parts.append("Direct Synapse/NYGC match; use curated NYGC values for clinical/diagnosis context and keep GEO as provenance.")
            if conflict_fields:
                note_parts.append("GEO-vs-Synapse differences noted: " + ";".join(conflict_fields) + ".")
            if warning_fields:
                note_parts.append("Additional label warnings: " + ";".join(warning_fields) + ".")
        row["resolved_conflict_note"] = " ".join(note_parts)

        if row["nygc_direct_match"] == "no":
            row["nygc_auxiliary_use_decision"] = "series_primary_metadata_only_no_direct_nygc_match"
            row["allowed_primary_clinical_fields"] = (
                "series_age_at_symptom_onset;series_site_of_motor_onset;"
                "series_disease_duration_onset_to_tracheostomy_or_death_months;"
                "series_duration_years_from_age"
            )
            row["forbidden_primary_use"] = "do_not_impute_nygc_fields_without_manual_id_mapping"
        elif conflict_fields:
            row["nygc_auxiliary_use_decision"] = "direct_match_synapse_nygc_standard_with_conflict_note"
            row["allowed_primary_clinical_fields"] = (
                "series_age_at_symptom_onset;series_site_of_motor_onset;"
                "series_disease_duration_onset_to_tracheostomy_or_death_months;"
                "series_duration_years_from_age;"
                "nygc_disease_group_rna;nygc_disease_duration_in_months;"
                "nygc_family_history_of_als_ftd;nygc_c9orf72_genotype;nygc_atxn2_genotype"
            )
            row["forbidden_primary_use"] = "do_not_claim_brain_ready_clinical_severity_pathology_or_cognition"
        else:
            row["nygc_auxiliary_use_decision"] = "direct_match_synapse_nygc_standard_for_context_sensitivity"
            row["allowed_primary_clinical_fields"] = (
                "series_age_at_symptom_onset;series_site_of_motor_onset;"
                "series_disease_duration_onset_to_tracheostomy_or_death_months;"
                "series_duration_years_from_age;"
                "nygc_disease_group_rna;nygc_disease_duration_in_months;"
                "nygc_family_history_of_als_ftd;nygc_c9orf72_genotype;nygc_atxn2_genotype"
            )
            row["forbidden_primary_use"] = "do_not_claim_brain_ready_clinical_severity_pathology_or_cognition"
        rows.append(row)
    return pd.DataFrame(rows)


def conflict_table(reconciled: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in reconciled.iterrows():
        fields = str(row.get("reconciliation_conflict_fields", "none"))
        if fields == "none":
            continue
        for field in fields.split(";"):
            rows.append(
                {
                    "mo_id": row["mo_id"],
                    "series_cgnd_id": row.get("series_cgnd_id", ""),
                    "field": field,
                    "series_value": row.get(series_field_for_conflict(field), ""),
                    "nygc_value": row.get(nygc_field_for_conflict(field), ""),
                    "decision": "use_synapse_nygc_value_with_note",
                    "resolved_conflict_note": row.get("resolved_conflict_note", ""),
                }
            )
    return pd.DataFrame(rows)


def series_field_for_conflict(field: str) -> str:
    return {
        "disease_label": "series_genotype",
        "broad_region": "series_tissue",
        "sex": "series_sex",
        "pmi": "series_pmi",
        "age_at_death": "series_age_at_death",
        "age_at_symptom_onset": "series_age_at_symptom_onset",
        "site_of_motor_onset": "series_site_of_motor_onset",
    }.get(field, "")


def nygc_field_for_conflict(field: str) -> str:
    return {
        "disease_label": "nygc_disease_group_rna",
        "broad_region": "nygc_tissue_clean",
        "sex": "nygc_sex_call",
        "pmi": "nygc_pmi",
        "age_at_death": "nygc_age_at_death",
        "age_at_symptom_onset": "nygc_age_at_symptom_onset",
        "site_of_motor_onset": "nygc_site_of_motor_onset",
    }.get(field, "")


def summary_table(reconciled: pd.DataFrame, nygc: pd.DataFrame) -> pd.DataFrame:
    als = reconciled[reconciled["series_genotype"].eq("ALS")]
    rows = [
        ("gse330130_mo_libraries", len(reconciled), "Collapsed from GEO series matrix sample columns."),
        ("gse330130_gsm_columns", int(pd.to_numeric(reconciled["n_gsm_columns"], errors="coerce").sum()), "ATAC and GEX are separate GSM columns; MO105 has GEX only."),
        ("gse330130_libraries_with_nygc_direct_match", int(reconciled["nygc_direct_match"].eq("yes").sum()), "Direct match uses series CGND-HRA id to NYGC rna_id."),
        ("gse330130_libraries_without_nygc_direct_match", int(reconciled["nygc_direct_match"].eq("no").sum()), "No Synapse/NYGC resolved fields are assigned without a direct RNA-ID match."),
        ("gse330130_libraries_with_reconciliation_conflict", int((~reconciled["reconciliation_conflict_fields"].eq("none")).sum()), "Any disease, broad-region, sex, PMI, age or onset-site difference; direct matches use Synapse/NYGC values with notes."),
        ("gse330130_als_libraries_with_series_onset_age", int(als["series_age_at_symptom_onset"].astype(str).str.len().gt(0).sum()), "GEO series matrix supports onset-age context for ALS libraries."),
        ("gse330130_als_libraries_with_series_onset_site", int(als["series_site_of_motor_onset"].astype(str).str.len().gt(0).sum()), "GEO series matrix supports onset-site context for ALS libraries."),
        (
            "gse330130_als_libraries_with_series_direct_duration_months",
            int(als[SERIES_DURATION_ALIAS].astype(str).str.len().gt(0).sum()) if SERIES_DURATION_ALIAS in als.columns else 0,
            "GEO series matrix has direct disease-duration months for a low-coverage ALS subset; use as descriptive/sensitivity context, not progression rate.",
        ),
        (
            "gse330130_als_libraries_with_series_age_derived_duration_years",
            int(als["series_duration_years_from_age"].astype(str).str.len().gt(0).sum()),
            "Age-at-death minus age-at-onset is coarse context and not a longitudinal progression endpoint.",
        ),
        ("nygc_rows_total", len(nygc), "Rows in NYGC clinical linkage CSV."),
        ("nygc_unique_donors", nygc["donor_id"].nunique(), "Broader NYGC/ST1-ST3 universe; not GSE330130-specific."),
        ("nygc_als_rows_with_duration_months", int(nygc[nygc["disease_group_rna"].eq("ALS")]["disease_duration_in_months"].astype(str).str.len().gt(0).sum()), "Useful clinical-context source for direct Synapse/NYGC matches."),
    ]
    return pd.DataFrame(rows, columns=["metric", "value", "interpretation"])


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


def write_report(reconciled: pd.DataFrame, conflicts: pd.DataFrame, summary: pd.DataFrame) -> Path:
    path = PLAN_ROOT / f"gse330130_nygc_clinical_linkage_reconciliation_{DATE_TAG}.md"
    summary_preview = markdown_table(summary)
    conflict_preview = markdown_table(conflicts.head(20))
    decision_counts = (
        reconciled.groupby("nygc_auxiliary_use_decision", dropna=False)
        .size()
        .reset_index(name="n_libraries")
        .sort_values("nygc_auxiliary_use_decision")
    )
    text = f"""# GSE330130 / NYGC clinical linkage reconciliation ({DATE_TAG})

## Bottom Line

This early reconciliation note is retained for provenance. Current P01 clinical/diagnosis interpretation uses the latest Synapse/NYGC-prioritized bridge when an exact RNA-ID match exists; GEO fields remain raw sequencing provenance and conflict notes.

GSE330130 GEO series matrix still confirms native onset-age/site context and limited direct duration fields. For direct NYGC matches, conflicts are no longer blockers: use Synapse/NYGC curated values as resolved clinical/diagnosis context, and keep `resolved_conflict_note` visible.

## Summary

{summary_preview}

## Current Use Decision

{markdown_table(decision_counts)}

Allowed now:

- Use Synapse/NYGC values as the resolved clinical/diagnosis standard for direct matched libraries, with GEO retained as provenance.
- Use GSE330130 series-matrix fields for libraries without a direct Synapse/NYGC match.
- Use direct series-matrix disease-duration months where present only as descriptive/sensitivity context.
- Use `age_at_death - age_at_symptom_onset` only as coarse disease-duration context, not as a precise progression-rate endpoint.
- Use NYGC duration/genotype/family-history/splice flags in direct matched libraries only as descriptive/sensitivity context.

Forbidden now:

- Do not use NYGC-only fields for unmatched MO libraries.
- Do not claim ALSFRS, cognitive, pathology, TDP-43, STMN2/UNC13A or Brain-ready clinical severity evidence from this reconciliation.

## Conflict Preview

{conflict_preview}

## Consequence For Target Journal Strategy

This improves the Brain-route metadata situation from \"no onset metadata\" to \"onset-age/site and partial duration context exists in GSE330130.\" It still does not make Brain first-ready because the current P01 package lacks ALSFRS(-R), progression rate, survival, cognitive endpoint in the main discovery cohort, neuropathology burden, TDP-43/splicing endpoint and an independent phenotype association with BH/FDR support.

## Machine-Readable Outputs

- `results/Supplemental_QC_covariate_audit/gse330130_nygc_clinical_reconciliation_{DATE_TAG}.tsv`
- `results/Supplemental_QC_covariate_audit/gse330130_nygc_clinical_reconciliation_conflicts_{DATE_TAG}.tsv`
- `results/Supplemental_QC_covariate_audit/gse330130_nygc_clinical_reconciliation_summary_{DATE_TAG}.tsv`
"""
    path.write_text(text, encoding="utf-8")
    return path


def main() -> None:
    series_path, series_label = resolve_input(
        SERIES_ENV,
        [
            source_data_root() / "GSE330130" / "GSE330130_series_matrix",
            ROOT / "metadata" / "GSE330130_series_matrix",
        ],
    )
    nygc_path, nygc_label = resolve_input(
        NYGC_ENV,
        [
            source_data_root() / "NYGC" / "nygc_downloaded_st3_st1_sample_clinical_linkage.csv",
            ROOT / "metadata" / "nygc_downloaded_st3_st1_sample_clinical_linkage.csv",
        ],
    )

    series = parse_series_libraries(series_path, series_label)
    nygc = read_nygc(nygc_path, nygc_label)
    reconciled = classify_rows(attach_nygc(series, nygc))
    conflicts = conflict_table(reconciled)
    summary = summary_table(reconciled, nygc)

    OUT.mkdir(parents=True, exist_ok=True)
    recon_path = OUT / f"gse330130_nygc_clinical_reconciliation_{DATE_TAG}.tsv"
    conflict_path = OUT / f"gse330130_nygc_clinical_reconciliation_conflicts_{DATE_TAG}.tsv"
    summary_path = OUT / f"gse330130_nygc_clinical_reconciliation_summary_{DATE_TAG}.tsv"
    write_tsv(reconciled, recon_path)
    write_tsv(conflicts, conflict_path)
    write_tsv(summary, summary_path)
    report_path = write_report(reconciled, conflicts, summary)

    print(f"reconciliation={recon_path} rows={len(reconciled)}")
    print(f"conflicts={conflict_path} rows={len(conflicts)}")
    print(f"summary={summary_path} rows={len(summary)}")
    print(f"report={report_path}")


if __name__ == "__main__":
    main()
