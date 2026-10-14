#!/usr/bin/env python3
"""Build Supplementary Tables S1-S4 for the v2 multi-cohort manuscript.

S1: frozen analysis plan and amendments, with commit hashes, times and SHA-256.
S2: registry of every test of the frozen modules and the consensus signature
    (W08 superseded module, pre-plan checks W50-W54, governed analyses W55-W62, W64, W66-W71).
S3: deviations from the plan and amendments.
S4: exploratory biological context (W63 pathways, W64 overlap with NYGC DE, W65 ligand-receptor pairs, W66 spatial, W67 motor-system selectivity, W68 independent spatial replication, W69 dose-response and cross-cohort combination, W70 microglia-adjusted territory test, W71 Maniatis spatial time course, W72 disease-state cell source).

All values are read from committed result tables and from git; nothing is re-computed.
"""

from __future__ import annotations

import hashlib
import logging
import subprocess
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
PLANS = ROOT / "docs" / "analysis_plans"
OUT = RES / "Supplementary_Tables_v2"
PLAN_STEM = "2026-10-11_P01_multicohort_astrocyte_consensus_analysis_plan"

S2_COLUMNS = [
    "analysis", "test_id", "phase", "cohort", "region", "contrast", "design", "gene_set",
    "n_case", "n_control", "n_genes_tested", "effect_type", "effect", "effect_ci_low", "effect_ci_high",
    "fry_p", "camera_p", "other_p_type", "other_p", "fdr", "fdr_family", "source_file",
]


def git(*args: str) -> str:
    """Run git in the repository and return stdout."""
    return subprocess.run(["git", "-C", str(ROOT), *args], check=True, capture_output=True, text=True).stdout.strip()


def commit_info(rev: str) -> tuple[str, str]:
    """Return (short hash, ISO commit time) for a revision."""
    h, t = git("show", "-s", "--format=%h|%ad", "--date=format:%Y-%m-%d %H:%M:%S %z", rev).split("|")
    return h, t


def first_commit(path: Path) -> str:
    """Hash of the commit that added a file; the file must be unchanged since."""
    revs = git("log", "--format=%h", "--", str(path.relative_to(ROOT))).split()
    if len(revs) != 1:
        raise RuntimeError(f"{path.name} changed after it was committed: {revs}")
    return revs[0]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def status_line(path: Path) -> str:
    """The 'Status:' line of a plan document, without markdown emphasis."""
    for line in path.read_text().splitlines():
        if line.startswith("Status:"):
            return line.removeprefix("Status:").replace("**", "").strip()
    return ""


# ---------------------------------------------------------------- S1
S1_ROWS = [
    # document suffix, version, scope, governed analysis, governed result commit
    (None, "Pre-plan", "External single-nucleus checks W50-W54, observed before the plan was frozen", "W50-W54", "267949d"),
    ("", "Plan v1.0", "Study question, hypotheses H1-H5, cohorts, confirmatory hold-out (Section 4), meta-analysis (Section 5), exploratory analyses (Section 6)", "W55-W59", None),
    ("_AMENDMENT_1", "Amendment 1", "NYGC hold-out implementation: discovery-donor exclusion, covariates, BH families", "W55 NYGC bulk hold-out (confirmatory)", "2861aef"),
    ("_AMENDMENT_2", "Amendment 2", "Exploratory cell-source checks (E1 specificity, E2 cell-type shift, E3 myeloid/vascular adjustment)", "W56 cell-source check (exploratory)", "9b1fd37"),
    ("_AMENDMENT_3", "Amendment 3", "Exploratory ambient-RNA / contamination test", "W57 ambient-RNA test (exploratory)", "8ce592d"),
    ("_AMENDMENT_4", "Amendment 4", "Section 5 implementation: per-cohort designs, standardised effect, Ulm age correction, consensus runs, reproducibility null", "W58 meta-analysis (consolidation)", "957a06a"),
    ("_AMENDMENT_5", "Amendment 5", "Second confirmatory test: frozen discovery-agnostic consensus signature in the NYGC hold-out", "W59 consensus-signature hold-out (confirmatory)", "fedddd5"),
    ("_AMENDMENT_6", "Amendment 6", "Restores the plan 5.3 sensitivity: beta_age estimated in donors aged 70 or older for the Ulm age correction", "W60 age >= 70 sensitivity (sensitivity of W58)", "7f16fd3"),
    ("_AMENDMENT_7", "Amendment 7", "Review-prompted checks: HKSJ small-k inference (already observed) and discovery confounder audit with tissue-programme sensitivity", "W61 HKSJ inference; W62 discovery audit (post hoc)", "9342a49"),
    ("_AMENDMENT_8", "Amendment 8", "Exploratory biological context: themed pathway GSEA, bulk hold-out against published deconvolution, ligand-receptor context", "W63 pathways; W64 published composition; W65 ligand-receptor (exploratory)", "73864d2"),
    ("_AMENDMENT_9", "Amendment 9", "Exploratory spatial context in GSE330224 Visium: donor pseudobulk disease tests, lineage niches, co-localisation, ligand-receptor co-detection", "W66 spatial context (exploratory)", "89c7c6a"),
    ("_AMENDMENT_10", "Amendment 10", "Review-prompted W66 sensitivity (empirical co-detection null, depth-partial co-localisation; already observed) and W67 motor- vs sensory-territory test (fixed before running)", "W66 post hoc; W67 motor-system selectivity (exploratory)", "96335a5"),
    ("_AMENDMENT_11", "Amendment 11", "Replication of motor-territory selectivity in an independent spatial cohort (GSE288365 Visium CytAssist cervical cord), with motor-neuron proximity and protein co-localisation", "W68 independent spatial replication (exploratory)", "7f9c562"),
    ("_AMENDMENT_12", "Amendment 12", "Dose-response of the motor-territory response with motor-neuron loss (fixed before running) and Stouffer combination of the W67 and W68 territory tests (post hoc summary of known results)", "W69 (exploratory; B post hoc)", "1c819a7"),
    ("_AMENDMENT_13", "Amendment 13", "Territory x disease interaction adjusted for a unit-level microglial score (all, homeostatic and activated marker sets), with prespecified combination and attenuation reading", "W70 microglia-adjusted territory test (exploratory)", "954b138"),
    ("_AMENDMENT_14", "Amendment 14", "SOD1-G93A spatial time course and descriptive human co-localisation in Maniatis et al. 2019 (W71); disease-state lineage share of module genes in GSE287257 (W72)", "W71 Maniatis spatial time course; W72 Zelic disease-state cell source (exploratory)", "833bb13"),
    ("_AMENDMENT_15", "Amendment 15", "Generic inflammation, TYROBP-network microglial and ageing gene sets in the territory test, and the ALS-up module adjusted for each generic score", "W73 generic signature comparison (exploratory)", "174391a"),
]


def build_s1() -> pd.DataFrame:
    rows = []
    for suffix, version, scope, governs, result_rev in S1_ROWS:
        row = {"document": "", "version": version, "scope": scope, "governs": governs, "commit": "", "commit_time": "",
               "sha256": "", "status_at_commit": "", "governed_result_commit": "", "governed_result_time": ""}
        if suffix is not None:
            path = PLANS / f"{PLAN_STEM}{suffix}.md"
            row["commit"], row["commit_time"] = commit_info(first_commit(path))
            row["document"] = f"docs/analysis_plans/{path.name}"
            row["sha256"] = sha256(path)
            row["status_at_commit"] = status_line(path)
        if result_rev:
            row["governed_result_commit"], row["governed_result_time"] = commit_info(result_rev)
            if suffix is None:
                row["commit"], row["commit_time"] = row["governed_result_commit"], row["governed_result_time"]
        rows.append(row)
    s1 = pd.DataFrame(rows)
    gov = s1[(s1["governed_result_time"] != "") & (s1["document"] != "")]
    if (gov["commit_time"] >= gov["governed_result_time"]).any():
        raise RuntimeError("An amendment is not committed before its governed analysis")
    return s1


# ---------------------------------------------------------------- S2 helpers
def read(rel: str) -> pd.DataFrame:
    return pd.read_csv(RES / rel, sep="\t")


def rows(df: pd.DataFrame, rel: str, **fixed: object) -> pd.DataFrame:
    """Attach constant columns and the source file, keep the S2 column order."""
    out = df.copy()
    for key, value in fixed.items():
        out[key] = value
    out["source_file"] = f"results/{rel}"
    return out.reindex(columns=S2_COLUMNS)


def gene_set_label(module: pd.Series) -> pd.Series:
    labels = {"up": "ALS-up module (15 genes)", "down": "ALS-down module (14 genes)",
              "up_consensus": "Consensus up (85 genes)", "down_consensus": "Consensus down (73 genes)",
              "up_astro_enriched": "ALS-up astrocyte-enriched subset"}
    return module.map(labels).fillna(module)


def w08() -> pd.DataFrame:
    rel = "W08_Pineda_validation/w08_primary_als_validation_summary.tsv"
    d = read(rel)
    d = d[d["module_id"].str.startswith("GSE330130_MTC_RNA__Astrocyte__")].copy()
    d["gene_set"] = d["direction"].map({"up": "Earlier ALS-up module", "down": "Earlier ALS-down module"}) + " (" + d["n_seed_genes"].astype(str) + " genes)"
    d = d.rename(columns={"model_formula_registered": "design", "n_genes_in_pineda": "n_genes_tested", "beta": "effect",
                          "p_value": "other_p", "adj_p_global_BH": "fdr"})
    d["test_id"] = d["analysis_type"]
    d["region"], d["contrast"] = "BA4 + BA9 donor mean", "ALS vs PN"
    d["n_case"], d["n_control"] = pd.NA, pd.NA
    return rows(d, rel, analysis="W08", phase="Superseded (earlier module definition, pre-plan)", cohort="Pineda 2024",
                effect_type="beta (module score)", other_p_type="OLS P", fdr_family="W08 global BH")


def w50() -> pd.DataFrame:
    rel = "W50_Pineda_final15_validation/w50_validation_results.tsv"
    d = read(rel)
    d["phase"] = "Pre-plan external check: " + d["tier"]
    d = d.rename(columns={"analysis_id": "test_id", "scope": "region", "covariates": "design", "n_case_donors": "n_case",
                          "n_control_donors": "n_control", "beta": "effect", "ci_low": "effect_ci_low",
                          "ci_high": "effect_ci_high", "p_value": "other_p", "BH_within_analysis": "fdr"})
    d["gene_set"] = gene_set_label(d["module"])
    a = rows(d, rel, analysis="W50", cohort="Pineda 2024", effect_type="beta (module score)", other_p_type="OLS P",
             fdr_family="W50 within analysis")
    rel2 = "W50_Pineda_final15_validation/w50_post_hoc_depth_adjusted.tsv"
    p = read(rel2).rename(columns={"covariates": "design", "beta": "effect", "ci_low": "effect_ci_low",
                                   "ci_high": "effect_ci_high", "p_value": "other_p"})
    p["test_id"] = p["analysis"] + ":" + p["response"]
    p["gene_set"] = gene_set_label(p["module"])
    b = rows(p, rel2, analysis="W50", phase="Pre-plan: post hoc", cohort="Pineda 2024", region="BA4", contrast="ALS vs PN",
             effect_type="beta (" + p["response"] + ")", other_p_type="OLS P")
    return pd.concat([a, b])


def limma_sets(rel: str, analysis: str, cohort: str, region: str, phase_map: dict[str, str], contrast: str,
               fry_col: str = "fry_p", case: str = "n_case", control: str = "n_control") -> pd.DataFrame:
    d = read(rel)
    d["phase"] = d["analysis_id"].map(phase_map)
    if d["phase"].isna().any():
        raise KeyError(f"Unmapped analysis_id in {rel}: {sorted(d.loc[d['phase'].isna(), 'analysis_id'].unique())}")
    d = d.rename(columns={"analysis_id": "test_id", "formula": "design", "n_module_genes_tested": "n_genes_tested",
                          fry_col: "fry_p", "fry_BH": "fdr", case: "n_case", control: "n_control"})
    d["gene_set"] = gene_set_label(d["module"])
    if "contrast" not in d:
        d["contrast"] = contrast
    return rows(d, rel, analysis=analysis, cohort=cohort, region=region, effect_type="", fdr_family="fry BH within analysis")


def w51() -> pd.DataFrame:
    phases = {"primary_discovery_matched": "Pre-plan external check: primary", "sensitivity_no_PMI": "Pre-plan: sensitivity",
              "sensitivity_min20_cells": "Pre-plan: sensitivity", "campaign2019_primary_no_PMI": "Pre-plan: post hoc",
              "campaign2019_with_PMI": "Pre-plan: post hoc"}
    a = limma_sets("W51_Pineda_BA4_discovery_matched_limma/w51_module_set_tests.tsv", "W51", "Pineda 2024", "BA4",
                   phases, "ALS vs PN", fry_col="fry_p_directional", case="n_als")
    b = limma_sets("W51_Pineda_BA4_discovery_matched_limma/w51_campaign2019_module_set_tests.tsv", "W51", "Pineda 2024",
                   "BA4 (2019 campaign)", phases, "ALS vs PN", fry_col="fry_p_directional", case="n_als")
    return pd.concat([a, b])


def w52() -> pd.DataFrame:
    phases = {"primary_donor_chip_blocked": "Pre-plan external check: primary", "sensA_donor_main_chip": "Pre-plan: sensitivity",
              "sensB_primary_plus_QC": "Pre-plan: sensitivity", "sensC_donor_no_chip": "Pre-plan: sensitivity",
              "context_ALS_FTD_nonC9": "Pre-plan: context", "context_C9_ALS_FTD": "Pre-plan: context"}
    a = limma_sets("W52_Ulm_motor_cortex_validation/w52_module_set_tests.tsv", "W52", "Ulm/DZNE", "M1", phases, "",
                   case="n_case_donors", control="n_control_donors")
    rel = "W52_Ulm_motor_cortex_validation/w52_post_hoc_age_bias_bound_summary.tsv"
    d = read(rel).rename(columns={"median_ulm_age_corrected": "effect"})
    d["test_id"] = "post_hoc_age_bias_bound"
    d["gene_set"] = gene_set_label(d["module"])
    b = rows(d, rel, analysis="W52", phase="Pre-plan: post hoc (bound, no test)", cohort="Ulm/DZNE", region="M1",
             contrast="ALS vs HC", effect_type="median age-corrected logFC")
    return pd.concat([a, b])


def w53() -> pd.DataFrame:
    rel = "W53_Masrori2025_astrocyte_DE_lookup/w53_module_direction_summary.tsv"
    d = read(rel).rename(columns={"comparison": "test_id", "n_reported": "n_genes_tested", "median_als_logFC": "effect",
                                  "sign_test_p_one_sided": "other_p"})
    d["region"] = d["test_id"].str.split("_").str[0].map({"MC": "Motor cortex", "SC": "Spinal cord"})
    d["contrast"] = d["test_id"].str.split("_", n=1).str[1] + " vs HC"
    d["gene_set"] = gene_set_label(d["module"])
    return rows(d, rel, analysis="W53", phase="Pre-plan external check: direction only", cohort="Masrori 2025",
                design="Published cell-level DE (summary)", n_case=5, n_control=5, effect_type="median logFC",
                other_p_type="one-sided sign test P")


def w54() -> pd.DataFrame:
    phases = {"primary_C9ALS_rule": "Pre-plan external check: primary", "sens_C9ALS_sex_only": "Pre-plan: sensitivity",
              "context_C9FTD_rule": "Pre-plan: context"}
    a = limma_sets("W54_Li2023_motor_cortex_validation/w54_module_set_tests.tsv", "W54", "Li 2023", "BA4", phases, "")
    rel = "W54_Li2023_motor_cortex_validation/w54_post_hoc_depth_sensitivity.tsv"
    d = read(rel).rename(columns={"formula": "design"})
    d["test_id"] = ["post_hoc_" + ("depth" if "depth" in f else "primary_design") for f in d["design"]]
    d["gene_set"] = gene_set_label(d["module"])
    b = rows(d, rel, analysis="W54", phase="Pre-plan: post hoc", cohort="Li 2023", region="BA4", contrast="C9-ALS vs Control",
             n_case=6, n_control=6)
    return pd.concat([a, b])


def bulk_sets(rel: str, analysis: str, phase: str, test_id: str | None = None, fdr_family: str = "",
              subset: str | None = None) -> pd.DataFrame:
    """NYGC bulk set tests; `subset` keeps one value of the table's `analysis` column."""
    d = read(rel)
    if subset is not None:
        d = d[d["analysis"] == subset].copy()
    d["test_id"] = d["analysis"] if test_id is None else test_id
    d = d.rename(columns={"n_als_donors": "n_case", "n_control_donors": "n_control", "n_module_genes_tested": "n_genes_tested",
                          "median_logFC": "effect", "fry_BH": "fdr"})
    d["region"] = d["region"].str.replace("_", " ")
    d["gene_set"] = gene_set_label(d["module"])
    d["design"] = "donor-blocked limma-voom; " + ("PMI in model" if d["pmi_in_model"].all() else "no PMI")
    return rows(d, rel, analysis=analysis, phase=phase, cohort="NYGC bulk hold-out", contrast="ALS vs control",
                effect_type="median logFC", fdr_family=fdr_family)


def w55() -> pd.DataFrame:
    a = bulk_sets("W55_NYGC_bulk_holdout/w55_primary_module_tests.tsv", "W55", "Confirmatory", "primary",
                  "fry BH over 4 tests")
    b = bulk_sets("W55_NYGC_bulk_holdout/w55_sensitivity_no_pmi_module_tests.tsv", "W55", "Prespecified sensitivity",
                  "sensitivity_no_PMI", "fry BH over 4 tests")
    rel = "W55_NYGC_bulk_holdout/w55_clinical_coefficients_with_ci.tsv"
    c = read(rel).merge(read("W55_NYGC_bulk_holdout/w55_clinical_association.tsv"), on=["region", "module"])
    c = c.rename(columns={"n": "n_case", "beta_per_sd": "effect", "ci_low": "effect_ci_low", "ci_high": "effect_ci_high",
                          "p": "other_p", "lm_BH": "fdr"})
    c["region"] = c["region"].str.replace("_", " ")
    c["gene_set"] = gene_set_label(c["module"])
    c = rows(c, rel, analysis="W55", test_id="clinical_duration_lm", phase="Prespecified secondary (H4)",
             cohort="NYGC bulk hold-out", contrast="Within ALS: disease duration",
             design="log(duration) ~ score + onset age + sex + onset site + C9orf72",
             effect_type="log-months per SD of score", other_p_type="linear model P", fdr_family="BH over 4 tests")
    return pd.concat([a, b, c])


def w56() -> pd.DataFrame:
    rel = "W56_cell_source_check/w56_e2_celltype_disease_shift.tsv"
    d = read(rel).rename(columns={"analysis": "test_id", "n_als": "n_case", "median_logFC": "effect"})
    pineda = d["test_id"].str.startswith("Pineda")
    d["cohort"] = pineda.map({True: "Pineda 2024 BA9 batch 210526", False: "Ulm/DZNE"})
    d["region"] = pineda.map({True: "BA9", False: "M1"})
    d["design"] = "per-lineage pseudobulk, " + d["test_id"].str.split("_").str[-1]
    d["gene_set"] = gene_set_label(d["module"])
    a = rows(d, rel, analysis="W56", phase="Exploratory (Amendment 2, E2)", contrast="ALS vs control",
             effect_type="median logFC")
    b = bulk_sets("W56_cell_source_check/w56_e3_nygc_composition_retest.tsv", "W56", "Exploratory (Amendment 2, E3)")
    return pd.concat([a, b])


def w57() -> pd.DataFrame:
    rel = "W57_ambient_contamination_test/w57_module_tests.tsv"
    d = read(rel).rename(columns={"analysis": "test_id", "median_logFC": "effect"})
    cohort = d["test_id"].str.split(":").str[0]
    d["cohort"] = cohort.map({"discovery": "NYGC GSE330130 (discovery)", "PinedaBA9_210526": "Pineda 2024 BA9 batch 210526",
                              "Ulm": "Ulm/DZNE (estimator assumption violated)"})
    d["region"] = cohort.map({"discovery": "MTC", "PinedaBA9_210526": "BA9", "Ulm": "M1"})
    d["design"] = d["test_id"].str.split(":").str[1]
    d["gene_set"] = gene_set_label(d["module"])
    return rows(d, rel, analysis="W57", phase="Exploratory (Amendment 3)", contrast="ALS vs control",
                effect_type="median logFC")


def w58() -> pd.DataFrame:
    rel = "W58_multicohort_meta/w58_module_effects_by_cohort.tsv"
    d = read(rel).rename(columns={"cohort": "test_id", "n_als": "n_case", "n_genes": "n_genes_tested", "d": "effect"})
    d["effect_ci_low"] = d["effect"] - 1.96 * d["se_d"]
    d["effect_ci_high"] = d["effect"] + 1.96 * d["se_d"]
    d["cohort"] = d["test_id"]
    d["region"] = d["region"].str.replace("_", " ")
    d["gene_set"] = gene_set_label(d["module"])
    a = rows(d, rel, analysis="W58", phase="Consolidation: per-cohort effect (Amendment 4)", contrast="ALS vs control",
             design="donor-level score ~ disease + cohort covariates", effect_type="standardised effect d")
    rel2 = "W58_multicohort_meta/w58_module_meta.tsv"
    m = read(rel2).rename(columns={"analysis": "test_id", "cohorts": "cohort", "pooled_d": "effect", "ci_low": "effect_ci_low",
                                   "ci_high": "effect_ci_high", "p": "other_p"})
    m["design"] = "REML random effects; k = " + m["k"].astype(str) + "; I2 = " + m["I2"].round(1).astype(str) + "%; criterion met: " + m["criterion_met"].astype(str)
    m["gene_set"] = gene_set_label(m["module"])
    b = rows(m, rel2, analysis="W58", phase="Consolidation: pooled (Amendment 4)", contrast="ALS vs control",
             effect_type="pooled d", other_p_type="meta-analysis P")
    return pd.concat([a, b])


def w59() -> pd.DataFrame:
    phase = {"1_primary_W55_design": "Confirmatory (Amendment 5)", "2_secondary_myeloid_vascular": "Prespecified secondary (Amendment 5)",
             "3_sensitivity_signature_cohort_overlap_excluded": "Prespecified sensitivity (Amendment 5)"}
    rel = "W59_NYGC_consensus_signature_holdout/w59_signature_tests.tsv"
    return pd.concat([bulk_sets(rel, "W59", label, key, "fry BH over 4 tests", subset=key) for key, label in phase.items()])


def w60() -> pd.DataFrame:
    rel = "W60_age70_sensitivity/w60_module_meta.tsv"
    m = read(rel).rename(columns={"pooled_d": "effect", "ci_low": "effect_ci_low", "ci_high": "effect_ci_high", "p": "other_p"})
    m["test_id"] = m["estimate"] + ":" + m["pool"]
    m["design"] = ("Ulm corrected with beta_age from " + m["estimate"] + "; REML random effects; k = " + m["k"].astype(str)
                   + "; I2 = " + m["I2"].round(1).astype(str) + "%; criterion met: " + m["criterion_met"].astype(str))
    m["gene_set"] = gene_set_label(m["module"])
    return rows(m, rel, analysis="W60", phase="Prespecified sensitivity (Amendment 6, run after W58)", cohort=m["pool"],
                contrast="ALS vs control", effect_type="pooled d", other_p_type="meta-analysis P")


def w61() -> pd.DataFrame:
    rel = "W61_small_k_meta_inference/w61_small_k_meta_inference.tsv"
    m = read(rel).rename(columns={"pooled_d": "effect", "hksj_ci_low": "effect_ci_low", "hksj_ci_high": "effect_ci_high",
                                  "hksj_p": "other_p"})
    m["test_id"] = m["pool"] + ":HKSJ"
    m["design"] = ("REML with Knapp-Hartung adjustment; k = " + m["k"].astype(str) + "; z-test P = " + m["z_p"].map("{:.3g}".format)
                   + "; criterion met (HKSJ): " + m["criterion_met_hksj"].astype(str))
    m["gene_set"] = gene_set_label(m["module"])
    return rows(m, rel, analysis="W61", phase="Post hoc, review-prompted (Amendment 7; observed before writing)", cohort=m["pool"],
                contrast="ALS vs control", effect_type="pooled d (HKSJ 95% CI)", other_p_type="HKSJ P")


def w62() -> pd.DataFrame:
    rel = "W62_discovery_confounder_audit/w62_module_set_tests.tsv"
    d = read(rel)
    d = d[d["model"] != "W34_primary_reproduction"].rename(
        columns={"model": "test_id", "n_als": "n_case", "median_logFC": "effect"})
    d["design"] = d["test_id"].map({"CGND_only": "CGND programme only; ~ disease + age + sex",
                                    "programme_adjusted": "all libraries; ~ disease + age + sex + PMI + programme"})
    d["gene_set"] = gene_set_label(d["module"])
    return rows(d, rel, analysis="W62", phase="Post hoc, review-prompted (Amendment 7; fixed before running)",
                cohort="NYGC GSE330130 (discovery)", region="MTC", contrast="ALS vs control", effect_type="median logFC")


def w64() -> pd.DataFrame:
    rel = "W64_bulk_published_composition/w64_set_tests.tsv"
    d = read(rel).rename(columns={"design": "test_id", "n_als_donors": "n_case", "n_control_donors": "n_control",
                                  "median_logFC": "effect"})
    d["design"] = d["test_id"].map({"a_st7_replaces_markers": "W55 design; marker scores replaced by NYGC ST7 proportions",
                                    "b_markers_plus_st7": "W55 design plus NYGC ST7 proportions"}) \
        + " (" + d["n_composition_covariates"].astype(str) + " covariates)"
    d["region"] = d["region"].str.replace("_", " ")
    d["gene_set"] = gene_set_label(d["gene_set"])
    return rows(d, rel, analysis="W64", phase="Exploratory (Amendment 8)", cohort="NYGC bulk hold-out",
                contrast="ALS vs control", effect_type="median logFC")


def w66() -> pd.DataFrame:
    rel = "W66_spatial_context/w66_set_tests.tsv"
    d = read(rel).rename(columns={"n_als": "n_case", "median_logFC": "effect", "fry_BH_primary_family": "fdr"})
    d["test_id"] = d["analysis"] + ":" + d["unit_type"] + ":" + d["stratum"]
    d["cohort"] = "GSE330224 Visium " + d["tissue"] + d["analysis"].map(lambda a: " (incl. GSE330130 donors)" if a == "sens_all_donors" else "")
    d["region"] = d["stratum"]
    d["phase"] = d["analysis"].map(lambda a: "Exploratory (Amendment 9)" if a == "primary" else "Exploratory sensitivity (Amendment 9)")
    d["gene_set"] = gene_set_label(d["gene_set"])
    d["fdr_family"] = d["fdr"].map(lambda x: "" if pd.isna(x) else "fry BH over 4 tests (frozen modules, whole unit, both tissues)")
    return rows(d, rel, analysis="W66", contrast="ALS vs non-neurological control", effect_type="median logFC")


def w67() -> pd.DataFrame:
    rel = "W67_motor_selectivity/w67_motor_selectivity_tests.tsv"
    d = read(rel).rename(columns={"n_als": "n_case", "median_logFC": "effect", "fry_BH_primary_family": "fdr"})
    d["test_id"] = d["comparison"] + ":" + d["contrast"]
    d["cohort"] = "GSE330224 Visium " + d["tissue"]
    d["region"] = d["comparison"]
    d["phase"] = d["family"].map({"primary": "Exploratory (Amendment 10 B, primary family)", "secondary": "Exploratory (Amendment 10 B, secondary)"})
    d["design"] = "~ territory * disease + sex + age, donor-blocked; contrast " + d["contrast"]
    d["gene_set"] = gene_set_label(d["gene_set"])
    d["fdr_family"] = d["fdr"].map(lambda x: "" if pd.isna(x) else "fry BH over 4 interaction tests")
    return rows(d, rel, analysis="W67", contrast="ALS vs control (territory interaction)", effect_type="median logFC")


def w68() -> pd.DataFrame:
    rel = "W68_GSE288365_replication/w68_replication_tests.tsv"
    d = read(rel).rename(columns={"n_als": "n_case", "median_logFC": "effect", "fry_BH_primary_family": "fdr"})
    d["test_id"] = d["analysis"] + ":" + d["comparison"] + ":" + d["contrast"]
    d["region"] = d["comparison"]
    d["phase"] = d["analysis"].map(lambda a: "Exploratory replication (Amendment 11, primary)" if a == "primary"
                                   else "Exploratory replication (Amendment 11, " + a + ")")
    d["design"] = "~ territory * disease, donor-blocked; contrast " + d["contrast"] + d["analysis"].map(lambda a: " + sex + age" if a == "sens_sex_age" else "")
    d["gene_set"] = gene_set_label(d["gene_set"])
    d["fdr_family"] = d["fdr"].map(lambda x: "" if pd.isna(x) else "fry BH over 4 interaction tests")
    return rows(d, rel, analysis="W68", cohort="GSE288365 Visium CytAssist cervical cord", contrast="ALS vs control (territory interaction)",
                effect_type="median logFC")


def w69() -> pd.DataFrame:
    rel_a = "W69_dose_response_and_combination/w69a_dose_response_by_cohort.tsv"
    a = read(rel_a).rename(columns={"n_als": "n_case", "rho": "effect", "p_one_sided": "other_p"})
    a["n_control"] = pd.NA
    a["test_id"] = "dose_response:" + a["cohort"]
    a["gene_set"] = gene_set_label(a["gene_set"])
    a = rows(a, rel_a, analysis="W69", phase="Exploratory (Amendment 12 A)", region="motor minus sensory grey",
             contrast="Within ALS: motor-neuron index vs module response", design="Spearman, within ALS donors",
             effect_type="Spearman rho (predicted < 0)", other_p_type="one-sided Spearman P")
    rel_b = "W69_dose_response_and_combination/w69b_interaction_combination.tsv"
    b = read(rel_b).rename(columns={"stouffer_z": "effect", "stouffer_p_one_sided": "other_p"})
    b["test_id"] = "stouffer:W67+W68"
    b["gene_set"] = gene_set_label(b["gene_set"])
    b = rows(b, rel_b, analysis="W69", phase="Post hoc summary of known results (Amendment 12 B)", cohort="GSE330224 + GSE288365",
             region="spinal grey matter", contrast="Territory x disease interaction, two cohorts combined",
             design="unweighted Stouffer of one-sided fry P", effect_type="Stouffer Z", other_p_type="combined one-sided P")
    return pd.concat([a, b])


def w70() -> pd.DataFrame:
    rel = "W70_microglia_adjusted_territory/w70_adjusted_interaction_tests.tsv"
    d = read(rel).rename(columns={"n_als": "n_case", "median_logFC_interaction": "effect"})
    d["test_id"] = d["cohort"] + ":" + d["score"] + "_adjusted:interaction"
    d["design"] = "W67/W68 interaction model + unit-level " + d["score"] + " marker score; attenuation vs unadjusted = " + d["attenuation"].round(2).astype(str)
    d["gene_set"] = gene_set_label(d["gene_set"])
    d["cohort"] = d["cohort"].map({"lumbar_GSE330224": "GSE330224 Visium LSC", "cervical_GSE288365": "GSE288365 Visium CytAssist cervical cord"})
    return rows(d, rel, analysis="W70", phase="Exploratory (Amendment 13)", region="motor vs sensory grey matter",
                contrast="Territory x disease interaction", effect_type="median logFC (interaction)")


def w71() -> pd.DataFrame:
    rel = "W71_maniatis_spatial_time_course/w71_mouse_tests.tsv"
    d = pd.concat([read(rel), read("W71_maniatis_spatial_time_course/w71_mouse_microglia_adjusted.tsv")], ignore_index=True).rename(columns={"n_G93A": "n_case", "n_WT": "n_control", "median_logFC": "effect", "fry_BH_family": "fdr"})
    d["test_id"] = d["age"] + ":" + d["comparison"] + ":" + d["contrast"] + d["adjusted"].map({True: ":microglia_adjusted", False: ""})
    d["phase"] = d.apply(lambda r: "Exploratory (Amendment 14, primary family)" if (r["gene_set"] == "up" and r["comparison"] == "grey"
                         and r["contrast"] == "interaction" and not r["adjusted"]) else "Exploratory (Amendment 14, secondary)", axis=1)
    d["region"] = d["comparison"].map({"grey": "ventral vs dorsal horn", "white": "ventral-lateral vs dorsal-medial white matter"})
    d["design"] = "~ territory * genotype + sex, animal-blocked; contrast " + d["contrast"] + d["adjusted"].map({True: "; + microglial score", False: ""})
    d["gene_set"] = d["gene_set"].map({"up": "ALS-up module, mouse orthologs (11 genes)", "up_consensus": "Consensus up, mouse orthologs (74 genes)"})
    d["fdr_family"] = d["fdr"].map(lambda x: "" if pd.isna(x) else "fry BH over 4 ages")
    out = rows(d, rel, analysis="W71", cohort="Maniatis 2019 ST, SOD1-G93A vs SOD1-WT mouse lumbar cord",
               contrast="G93A vs WT", effect_type="median logFC")
    out.loc[d["adjusted"].astype(bool).to_numpy(), "source_file"] = "results/W71_maniatis_spatial_time_course/w71_mouse_microglia_adjusted.tsv"
    return out


def w73() -> pd.DataFrame:
    rel = "W73_generic_signature_comparison/w73_module_adjusted_tests.tsv"
    d = read(rel).rename(columns={"n_als": "n_case", "median_logFC_interaction": "effect"})
    d["test_id"] = d["cohort"] + ":" + d["adjusted_for"] + "_adjusted:interaction"
    d["design"] = "W67/W68 interaction model + unit-level " + d["adjusted_for"] + " score; attenuation vs unadjusted = " + d["attenuation"].round(2).astype(str)
    d["gene_set"] = gene_set_label(d["gene_set"])
    d["cohort"] = d["cohort"].map({"lumbar_GSE330224": "GSE330224 Visium LSC", "cervical_GSE288365": "GSE288365 Visium CytAssist cervical cord"})
    return rows(d, rel, analysis="W73", phase="Exploratory (Amendment 15)", region="motor vs sensory grey matter",
                contrast="Territory x disease interaction", effect_type="median logFC (interaction)")


def build_s2() -> pd.DataFrame:
    parts = [w08(), w50(), w51(), w52(), w53(), w54(), w55(), w56(), w57(), w58(), w59(), w60(), w61(), w62(), w64(), w66(), w67(), w68(), w69(), w70(), w71(), w73()]
    return pd.concat([p.astype(object) for p in parts], ignore_index=True)


# ---------------------------------------------------------------- S3
S3_ROWS = [
    ("Plan 5.2 (effect size)", "Hedges' g", "d = beta_disease / sigma_residual from the donor-level model",
     "Changed by Amendment 4 before any pooled estimate", "W58", "docs/analysis_plans/" + PLAN_STEM + "_AMENDMENT_4.md"),
    ("Plan 5.1 (Zelic astrocytes)", "Author-annotated subset plus marker rule", "Marker rule applied to all nuclei, one pseudobulk per donor",
     "Changed by Amendment 4 before any pooled estimate", "W58", "docs/analysis_plans/" + PLAN_STEM + "_AMENDMENT_4.md"),
    ("Plan 5.3 (age-correction sensitivity)", "Re-estimate beta_age in donors aged 70 or older",
     "Omitted from Amendment 4 and W58; run later as W60 under Amendment 6, with only Pineda BA4 and BA9 eligible",
     "Unplanned omission found while compiling this table; restored after W58 results were known, with the reading fixed in Amendment 6 before running",
     "W58, W60", "docs/analysis_plans/" + PLAN_STEM + "_AMENDMENT_6.md"),
    ("Amendment 6 (gene-level age correction)", "Age >= 70 sensitivity of the gene-level (H3) Ulm correction",
     "Not run", "A4.3 pools gene-level age coefficients over at least three cohorts; only two cohorts were eligible", "W60",
     "docs/analysis_plans/" + PLAN_STEM + "_AMENDMENT_6.md"),
    ("Plan 5.5 (reproducibility null)", "Expression- and detection-matched random gene sets",
     "Matched on source-cohort AveExpr deciles only", "Simplified by Amendment 4", "W58",
     "docs/analysis_plans/" + PLAN_STEM + "_AMENDMENT_4.md"),
    ("Amendment 4 A4.4 (consensus run)", "Run (i), all six cohorts including discovery, is the primary consensus (316 genes)",
     "The manuscript reports run (ii), discovery excluded (159 genes), as the consensus signature; runs (i) and (iii) are in the W58 tables",
     "Run (ii) was chosen for the second confirmatory test in Amendment 5, before testing, because it is independent of the discovery cohort",
     "W58, W59", "results/W58_multicohort_meta/w58_consensus_genes.tsv"),
    ("Plan 5.1-5.2 (independence)", "Pineda BA4 and BA9 treated as independent cohorts",
     "They share 29 of 30 donors; the secondary pool and consensus run (ii) count these donors twice",
     "Found after W58; disclosed. The primary motor-cortex pool and run (iii) are unaffected", "W58, W59",
     "results/W58_multicohort_meta/README.md"),
    ("Plan Section 4 (ALS-FTD sensitivity)", "ALS-FTD vs control as a sensitivity analysis", "Not run",
     "No ALS-FTD sample survived QC and discovery-donor exclusion", "W55", "results/W55_NYGC_bulk_holdout/README.md"),
    ("Amendment 3 (log rho covariate)", "log(rho) as a covariate", "log(rho + half the minimum positive rho) in all cohorts",
     "One Ulm unit had rho = 0", "W57", "results/W57_ambient_contamination_test/README.md"),
    ("Amendment 3 (Ulm correction)", "Contamination correction in Ulm", "Run, but not interpretable",
     "Estimated rho exceeded 1 in many Ulm units (median 0.95, maximum 3.1), so the estimator's assumption fails", "W57",
     "results/W57_ambient_contamination_test/README.md"),
    ("Amendment 2 E3 (astrocyte-enriched subset)", "Test the astrocyte-enriched subset of the ALS-up module in bulk",
     "Only EGF was testable in bulk", "AL133482.1 was not measured after filtering; the subset test is a single-gene test", "W56",
     "results/W56_cell_source_check/w56_e3_nygc_composition_retest.tsv"),
    ("Plan Section 6 (E1 ATAC, E2 spatial, E3 down-module biology, E4 moderators)", "Exploratory analyses", "Not run",
     "Not run before submission; no results are reported or implied", "-", "docs/analysis_plans/" + PLAN_STEM + ".md"),
    ("Amendment 4 A4.2 (pooling inference)", "REML random effects with a z test",
     "Also reported with the Knapp-Hartung adjustment (W61); no pool meets the A4.2 criterion under HKSJ",
     "Post hoc, prompted by an internal simulated review; HKSJ results were observed before Amendment 7 was written. The prespecified z-test result is retained",
     "W58, W60, W61", "results/W61_small_k_meta_inference/README.md"),
    ("Discovery confounder audit (not in plan)", "None; the plan audited only external cohorts",
     "W62 audit and tissue-programme sensitivity (5 of 7 controls from a U54 programme with no ALS)",
     "Post hoc, review-prompted; design and reading fixed in Amendment 7 before running", "W62",
     "results/W62_discovery_confounder_audit/README.md"),
    ("Amendment 8 W65 (complex score)", "Ligand and receptor log2FC per sender-receiver pair",
     "For multi-subunit ligands or receptors, the log2FC is the mean over subunits (expression uses the minimum)",
     "Implementation detail not fixed in Amendment 8; set when writing the script, before any result", "W65",
     "results/W65_ligand_receptor_context/README.md"),
    ("Amendment 8 W64 (data handling)", "ST7 proportions as covariates",
     "Text-stored workbook values coerced to numeric; subregion levels left empty after the ST7 restriction dropped",
     "Needed for the model to run; no outcome involved", "W64", "results/W64_bulk_published_composition/README.md"),
    ("Amendment 9 W66 (lineage niches)", "Niche = lineage score >= section 75th percentile; astrocyte_only also requires vascular and microglia scores below the section median",
     "Niches use spots with >= 500 total counts; 'rich' also requires a score > 0; 'below the median' became 'at or below the median'",
     "Markers are zero in most low-count spots (median vascular and microglia scores 0 in 95 of 137 LSC sections), so the written rule left astrocyte_only nearly empty; changed after building the inputs and before any disease test",
     "W66", "results/W66_spatial_context/README.md"),
    ("W66 ligand-receptor co-detection null (not in Amendment 9)", "Observed / expected under independence of the two genes",
     "Also compared with 200 detection-matched random gene pairs per pair (empirical P); claims restricted to empirical P < 0.05",
     "Post hoc, prompted by an internal simulated review; computed and observed before Amendment 10 fixed its implementation", "W66",
     "results/W66_spatial_context/post_hoc/README.md"),
    ("W66 co-localisation depth sensitivity (not in Amendment 9)", "Spearman correlation across spots",
     "Also with ranks residualised on log spot depth", "Post hoc, review-prompted; observed before Amendment 10", "W66",
     "results/W66_spatial_context/post_hoc/README.md"),
    ("W69 B cross-cohort combination", "Not planned",
     "Unweighted Stouffer combination of the W67 and W68 one-sided territory-interaction P values",
     "Post hoc summary; both inputs were known when Amendment 12 was written; reported as a summary, not a new test", "W69",
     "results/W69_dose_response_and_combination/README.md"),
    ("W72 within-lineage model", "Listed outputs: one W72 script",
     "The within-lineage ALS vs control model is run by a second script, scripts/run_w72_lineage_de.R, because it uses limma",
     "Implementation detail; no analytic change", "W72", "results/W72_zelic_disease_cell_source/README.md"),
    ("Pre-plan W50 (M1 controls)", "BA4 plus M1 controls sensitivity", "Dropped",
     "The five M1 control samples are C57BL/6 mouse tissue", "W50", "results/W50_Pineda_final15_validation/README.md"),
]


def build_s3() -> pd.DataFrame:
    return pd.DataFrame(S3_ROWS, columns=["item", "planned", "done", "reason_and_timing", "analysis", "source_file"])


# ---------------------------------------------------------------- S4
def build_s4() -> dict[str, pd.DataFrame]:
    """Exploratory context tables, copied from the committed W63-W65 outputs."""
    gsea = pd.concat([read("W63_pathway_context/w63_themed_gsea.tsv"), read("W63_pathway_context/w63_hallmark_gsea.tsv")])
    return {"S4a_W63_pathway_gsea": gsea,
            "S4b_W63_consensus_ora": read("W63_pathway_context/w63_consensus_ora.tsv"),
            "S4c_W64_modules_in_NYGC_ST4": read("W64_bulk_published_composition/w64_module_genes_in_nygc_st4.tsv"),
            "S4d_W65_pairs_both_cohorts": read("W65_ligand_receptor_context/w65_pairs_in_both_cohorts.tsv"),
            "S4e_W66_region_share": read("W66_spatial_context/w66_region_share.tsv"),
            "S4f_W66_colocalisation": read("W66_spatial_context/w66_colocalisation_summary.tsv"),
            "S4g_W66_lr_codetection": read("W66_spatial_context/w66_lr_codetection_summary.tsv").merge(
                read("W66_spatial_context/post_hoc/w66_lr_codetection_empirical_null.tsv").drop(columns="observed"), on=["tissue", "pair"]),
            "S4h_W67_territory_genes": read("W67_motor_selectivity/w67_module_gene_results.tsv"),
            "S4i_W68_territory_genes": read("W68_GSE288365_replication/w68_module_gene_results.tsv"),
            "S4j_W68_protein_coloc": read("W68_GSE288365_replication/w68_protein_colocalisation_summary.tsv"),
            "S4k_W69_dose_response_donors": read("W69_dose_response_and_combination/w69a_donor_values.tsv"),
            "S4l_W70_microglia_adjusted": read("W70_microglia_adjusted_territory/w70_combined_reading.tsv"),
            "S4m_W70_microglial_score": read("W70_microglia_adjusted_territory/w70_microglial_score_by_territory.tsv"),
            "S4n_W71_mouse_gene_effects": read("W71_maniatis_spatial_time_course/w71_mouse_gene_effects.tsv"),
            "S4o_W71_human_coloc": read("W71_maniatis_spatial_time_course/w71_human_colocalisation_summary.tsv"),
            "S4p_W72_lineage_share": read("W72_zelic_disease_cell_source/w72_lineage_share.tsv"),
            "S4q_W72_within_lineage_DE": read("W72_zelic_disease_cell_source/w72_within_lineage_de.tsv"),
            "S4r_W73_module_adjusted": read("W73_generic_signature_comparison/w73_module_adjusted_combined.tsv"),
            "S4s_W73_generic_sets": read("W73_generic_signature_comparison/w73_generic_set_territory_tests.tsv").merge(
                read("W73_generic_signature_comparison/w73_generic_set_combined.tsv"), on="gene_set")}


README_ROWS = [
    ("Supplementary Table S1", "Frozen analysis plan and fifteen amendments: version-control commit, commit time, SHA-256 and the analysis each governs. Each amendment's commit precedes its governed result commit."),
    ("Supplementary Table S2", "Registry of every test of the frozen modules and of the consensus signature, including the superseded W08 module test, pre-plan checks (W50-W54) and governed analyses (W55-W62, W64, W66-W71, W73). phase gives the status of each test; null and post hoc results are included."),
    ("Supplementary Table S3", "Deviations from the frozen plan and amendments, with reasons and timing."),
    ("Supplementary Table S4", "Exploratory biological context (Amendments 8-15): S4a themed and Hallmark preranked GSEA on the discovery-agnostic meta statistic (W63); S4b over-representation in the consensus arms (W63); S4c frozen module genes in the NYGC consortium's bulk DE, Flagship ST4, which includes discovery donors (W64); S4d ligand-receptor pairs involving module genes, evaluable in both Pineda BA9 and Ulm (W65; no P values); S4e-S4g GSE330224 Visium region share, spot co-localisation and ligand-receptor co-detection with its matched-random-pair null (W66, Amendments 9 and 10); S4h gene-level ALS effects in motor and sensory spinal-cord territories (W67, Amendment 10); S4i the same in the independent GSE288365 cohort and S4j spot-level protein co-localisation there (W68, Amendment 11); S4k donor-level motor-neuron index and module response (W69, Amendment 12); S4l-S4m microglia-adjusted territory interaction and microglial score by territory (W70, Amendment 13); S4n gene-level SOD1-G93A effects by age and territory and S4o human co-localisation in Maniatis et al. 2019 (W71, Amendment 14); S4p disease-state lineage share and S4q within-lineage ALS vs control effects of the module genes in GSE287257 (W72, Amendment 14); S4r the ALS-up module and consensus up arm adjusted for generic inflammation, TYROBP-network and ageing scores, and S4s those generic sets in the territory test (W73, Amendment 15)."),
    ("Repository", f"Generated by scripts/build_v2_supplementary_tables.py at commit {git('rev-parse', '--short', 'HEAD')}."),
]


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    tables = {"S1_plan_and_amendments": build_s1(), "S2_test_registry": build_s2(), "S3_deviations": build_s3(), **build_s4()}
    for name, df in tables.items():
        df.to_csv(OUT / f"Supplementary_Table_{name}.tsv", sep="\t", index=False)
        logger.info("%s: %d rows", name, len(df))
    readme = pd.DataFrame(README_ROWS, columns=["table", "description"])
    with pd.ExcelWriter(OUT / "Supplementary_Tables_S1-S4.xlsx") as xl:
        readme.to_excel(xl, sheet_name="README", index=False)
        for name, df in tables.items():
            df.to_excel(xl, sheet_name=name.split("_")[0], index=False)
    logger.info("written to %s", OUT)


if __name__ == "__main__":
    main()
