# W73: generic inflammation, microglial activation and ageing signatures in the territory test (Amendment 15)

Exploratory. Amendment 15 was committed before this analysis (1fca433).

Script: `scripts/run_w73_generic_signature_comparison.R`.

**Inputs.** Units, donors and models are identical to the primary territory tests (W67 lumbar, W68 cervical) and to W70.

**Gene sets.** Six msigdbr 26.1.0 sets, with module and consensus-up genes removed. The size after removal is given in brackets:
- HALLMARK_INFLAMMATORY_RESPONSE (199)
- HALLMARK_TNFA_SIGNALING_VIA_NFKB (198)
- HALLMARK_INTERFERON_GAMMA_RESPONSE (200)
- HALLMARK_COMPLEMENT (194)
- WP_TYROBP_CAUSAL_NETWORK_IN_MICROGLIA (59)
- LU_AGING_BRAIN_UP (258)

## A. Generic sets in the territory × disease test

| Set | Lumbar P (FDR) | Cervical P (FDR) | Combined one-sided P | Reading |
|---|---|---|---|---|
| TYROBP microglial network | 0.0045 (0.027) | 0.0056 (0.026) | 3.7 × 10⁻⁵ | Shares the ventral-horn pattern |
| Complement | 0.026 (0.052) | 0.013 (0.026) | 4.4 × 10⁻⁴ | Shares |
| Inflammatory response | 0.063 (0.094) | 0.017 (0.026) | 0.0014 | Shares |
| Ageing brain (up) | 0.010 (0.031) | 0.094 (0.094) | 0.0014 | Shares |
| TNFα–NF-κB | 0.17 (0.17) | 0.016 (0.026) | 0.0037 | Shares |
| IFN-γ response | 0.10 (0.13) | 0.032 (0.038) | 0.0039 | Shares |

All six programmes rise more in motor than in sensory grey matter in ALS. The largest interaction is for the TYROBP network (median log₂FC 0.38 lumbar, 0.62 cervical).

## B. ALS-up module adjusted for each generic score

| Adjusted for | Attenuation, lumbar / cervical | Combined one-sided P | Reading |
|---|---|---|---|
| Inflammatory response | 31% / −1% | 0.0055 | Beyond this programme |
| IFN-γ response | 38% / 3% | 0.0032 | Beyond |
| TNFα–NF-κB | 28% / 24% | 0.0069 | Beyond |
| Complement | 64% / 35% | 0.018 | Partly shared |
| Ageing brain (up) | 81% / 29% | 0.038 | Partly shared |
| TYROBP microglial network | 83% / 68% | 0.075 | Largely shared |

**Overall prespecified reading: no.** The module is not reducible to inflammatory, IFN-γ or TNFα signalling. It is partly shared with complement and ageing programmes, and largely shared with the TYROBP microglial network (as in W70).

**Consensus up arm.** Largely shared with five of the six programmes. It is partly shared with the TYROBP network (combined P = 0.033, maximum attenuation 55%) (Supplementary Table S4r).

## Interpretation

- **The region.** In ALS spinal cord, a broad glial and immune activation (complement, inflammatory, interferon and TNFα, TYROBP-network and ageing-associated genes) is concentrated in the ventral horn. This is consistent with known ventral-horn microgliosis.
- **The module.** The module is one readout of that state. Its regional rise goes beyond generic inflammatory signalling but is mostly shared with TYROBP-network microglial activation.
- **Not tested.** Disease specificity (ALS against other neurodegenerative diseases) was not tested. NYGC has only 11 "Other MND" and 5 "Other" motor-cortex donors.

## Files

- `w73_generic_set_territory_tests.tsv`
- `w73_generic_set_combined.tsv`
- `w73_module_adjusted_tests.tsv`
- `w73_module_adjusted_combined.tsv`
- `w73_reading.tsv`
- `w73_gene_sets.tsv`
- `w73_run.log`
