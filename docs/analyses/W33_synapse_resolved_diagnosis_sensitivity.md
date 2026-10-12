# W33 Synapse/NYGC-resolved diagnosis sensitivity

W33 does not replace frozen W01/W03. It reuses their lead-context counts, feature filtering, TMM, limma-voom and covariate-selection rules.
For exact NYGC ST3 RNA-ID matches, resolved diagnosis is authoritative: ALS and Control remain primary; ALS-FTD and Other are excluded from an ALS-versus-Control sensitivity model.
For libraries without an exact match, original GEO labels are retained. All changed/excluded libraries, model membership and frozen-versus-W33 feature comparisons are written here.
This is required because the frozen analyses used source-object GEO labels while later clinical reconciliation identified exact-match diagnosis conflicts.
