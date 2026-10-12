# W53 Masrori 2025 astrocyte DE lookup

Date: 2026-10-10
Script: `scripts/run_w53_masrori_astrocyte_lookup.py`

## Source

Masrori et al. 2025, *Nat Neurosci* (PMC12586167; doi 10.1038/s41593-025-02075-1). The cohort is Leuven/NBB/IDIBAPS, with matched motor cortex (MC) and spinal cord (SC) from 5 sALS, 5 C9-ALS and 5 controls. The local copy of the Supplementary Tables 1-32 workbook (32.4 MB) is in `Data/Other_Datasets/Masrori2025_NatNeuro/files/`. It was downloaded from Springer static content because the PMC link returned an HTML interstitial.

The raw human data are restricted (EGA EGAS00001006711). `GSE302295_normalized_counts.txt.gz` (local `Data/GSE302295/`) contains 30 samples from the ALS *models* series (iPSC/xenograft microglia) and is not human tissue, so it is not used.

## Caveats

- The tables are Seurat-style cell-level DE (`p_val`, `avg_log2FC`, `pct.1`, `pct.2`). Their p values treat nuclei as replicates, so only effect directions are used.
- `pct.1` and `pct.2` equal 1 for every gene, which suggests the test was run on smoothed or imputed values.
- Positive `avg_log2FC` means higher in healthy, so the ALS effect is `-avg_log2FC` (sign check on SLC1A2, NTRK2, GJA1, SERPINA3, CHI3L1 and VIM).

## Results (`w53_module_direction_summary.tsv`)

| Comparison | Up genes concordant (sign P) | Down genes concordant (sign P) |
|---|---|---|
| MC sALS vs healthy (ST25) | 8/15 (0.50) | 13/14 (0.0009) |
| SC sALS vs healthy (ST28) | 6/15 (0.85) | 12/14 (0.0065) |
| MC C9-ALS vs healthy (ST24) | 9/15 (0.30) | 13/14 (0.0009) |
| SC C9-ALS vs healthy (ST27) | 5/15 (0.94) | 10/14 (0.090) |

Up-module effects are near zero. The largest is C3 +0.25 in SC sALS, and LTF, GPNMB, EGF and ANXA1 are all within ±0.1.

## Interpretation boundary

This is summary-level, 5-vs-5, cell-level evidence. It gives direction only and cannot count as replication. It agrees with Pineda in pointing to the ALS-down module as the more reproducible component, while the ALS-up module (C3, LTF, ...) is not supported.
