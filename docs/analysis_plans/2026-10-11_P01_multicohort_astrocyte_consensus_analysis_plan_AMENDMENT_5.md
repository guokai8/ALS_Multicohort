# Amendment 5: confirmatory test of the discovery-agnostic consensus signature in the NYGC bulk hold-out

Date: 2026-10-11
Amends: plan v1.0 (sha256 da63b81a...bfb3bd); adds a second confirmatory test
Status: committed **before** the consensus signature was tested in NYGC bulk expression data. The NYGC bulk hold-out (W55 sample set) has been used only to test the frozen W34 modules. It played no part in defining this signature.

## Rationale

W58 run (ii) defined 159 consensus genes from five independent single-nucleus astrocyte cohorts (Li, Ulm age-corrected, Pineda BA4 2019, Pineda BA9 batch 210526, Zelic) without the discovery cohort. NYGC bulk did not contribute, so testing this signature there is prospective.

## Frozen signature

- File: `results/W59_NYGC_consensus_signature_holdout/w59_frozen_consensus_signature.tsv`
- sha256: `9a4bd65b0d17995081f5ff7f4071522b447d4aea3f103258eced23ef0c1480ce`
- Source: `results/W58_multicohort_meta/w58_consensus_genes.tsv`, `analysis == ii_without_discovery`.
- Direction comes from the sign of the meta logFC.
- Genes that are also W55 or W56-E3a composition markers are removed before testing so that no covariate is tested as an outcome. This removed SPP1.
- Final: **85 up and 73 down genes (158)**.

## Tests (NYGC bulk, W55 inputs and helpers)

1. **Primary.** W55 primary complete-case design with blocked fry, testing `up_consensus` and `down_consensus` in motor cortex and spinal cord. BH is applied over the 4 tests.
   - **Confirmation criterion** (as in plan Section 4): a set is confirmed if fry BH < 0.05 in its hypothesised direction (up: Up; down: Down) in at least one region.
   - Whether both regions are confirmed is also reported.
2. **Secondary.** The same test with the W56-E3a design (adding activated-myeloid and vascular scores).
3. **Prespecified sensitivity for donor overlap with signature-defining cohorts.**
   - Li, Pineda and Zelic donors could also be NYGC consortium donors, and IDs cannot be linked. Ulm (European banks) cannot overlap.
   - Using the Amendment 1 fingerprint rule (diagnosis class, sex, age at death ±1; "90 or older" matches 89 or above), exclude every NYGC donor that matches any Li (W54), Pineda (W51 BA4 and W56 BA9 samples used in W58) or Zelic (W58) donor.
   - Then repeat test 1.
4. **Gene level.** Within each set and region, report the fraction of genes whose NYGC logFC has the hypothesised sign, with a one-sided binomial test.

## Interpretation

If the primary criterion is met for at least one set, the manuscript may state that a single-nucleus-derived cross-cohort consensus signature was validated in an independent bulk hold-out. Tests 2 and 3 qualify that statement (composition, donor overlap); they do not replace it.
