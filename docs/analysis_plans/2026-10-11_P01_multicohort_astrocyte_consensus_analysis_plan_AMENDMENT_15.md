# Amendment 15: is the ventral-horn response generic neuroinflammation? (W73)

Date: 2026-10-14
Amends: plan v1.0 (sha256 da63b81a...bfb3bd)
Status: committed **before** W73 was run. Exploratory. Prompted by external pre-submission review, which asked whether the module represents generic neuroinflammation or microglial activation.

## Inspected before writing (no outcome)

- The gene-set sizes in msigdbr 26.1.0.
- The overlap of each set with the frozen ALS-up module (15 genes) and the consensus up arm (85 genes):

| Set | Genes | ALS-up overlap | Consensus-up overlap |
|---|---|---|---|
| HALLMARK_INFLAMMATORY_RESPONSE | 200 | EMP3 | 0 |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | 200 | none | 2 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 200 | none | 0 |
| HALLMARK_COMPLEMENT | 200 | LTF, ANG, C3 | 3 |
| WP_TYROBP_CAUSAL_NETWORK_IN_MICROGLIA | 60 | C3 | 0 |
| LU_AGING_BRAIN_UP | 263 | none | 5 |

- NYGC disease-group counts. The metadata list only 11 "Other MND" and 5 "Other" motor-cortex donors, too few for a disease-specificity test.

No expression summary has been computed for these sets in any cohort.

## Gene sets

The six sets above (msigdbr 26.1.0, human).
- The TYROBP network is the available proxy for disease-associated microglial activation. MSigDB does not contain the Keren-Shaul DAM list.
- LU_AGING_BRAIN_UP is the ageing set.
- Before testing, module and consensus-up genes are removed from every generic set, so that the generic sets and the module do not share genes.

## Design

Units, donors and models are identical to the primary territory tests:
- **Lumbar cord** (GSE330224): ventral vs dorsal horn, `~ territory * disease + sex + age`.
- **Cervical cord** (GSE288365): anterior horn vs other grey matter, `~ territory * disease`.

Both are donor-blocked (duplicateCorrelation, two passes), with fry as the primary test and camera as the secondary.

**A. Generic sets in the territory test.** The territory × disease interaction is tested for each generic set in each cohort.
- Within each cohort, BH correction is applied across the six sets.
- The two cohorts are combined by unweighted Stouffer on one-sided P values in the Up direction.

**B. Module adjusted for each generic score.**
- The frozen ALS-up interaction is refitted with the centred unit-level mean log2-CPM (TMM, prior count 1) of each generic set added as a covariate, as in W70.
- Attenuation is defined as in W70, against the W67/W68 primary interaction.
- The two cohorts are combined by Stouffer.

## Readings

**Per generic set (A).** "Shares the ventral-horn pattern" if the combined one-sided P < 0.05 and the interaction is Up in both cohorts. Otherwise "does not share it".

**Per generic set (B), for the module:**

| Combined adjusted P | Interaction | Attenuation | Reading |
|---|---|---|---|
| < 0.05 | Up in both cohorts | < 50% in both | Beyond this programme |
| < 0.05 | Up in both cohorts | ≥ 50% in either | Partly shared |
| ≥ 0.05, or not Up in both | | | Largely shared |

**Overall.** The module is "not reducible to generic inflammation or ageing" if reading B is "beyond" for all four Hallmark sets and the ageing set. The TYROBP network is reported separately, because W70 already showed that much of the rise is shared with activated-myeloid markers.

**Secondary.** The same B analysis for the consensus up arm. Every outcome is reported.

## Limits stated in advance

These analyses test generic activation and ageing programmes. They do not test disease specificity, that is, ALS against other neurodegenerative diseases.

## Outputs

- `results/W73_generic_signature_comparison/`
- `scripts/run_w73_generic_signature_comparison.R`
