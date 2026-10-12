#!/usr/bin/env python3
"""Check the W54 marker-dominance astrocyte rule against Ulm WNN labels.

Applies `w54_marker_annotation.call_astrocytes` to every Ulm/DZNE nucleus
(Zenodo 18370646) and compares it with the authors' WNN_L2 == "Astrocytes"
label. This validates the annotation rule on an independently annotated
motor-cortex cohort before it is used on unannotated Li 2023 data.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from p01_paths import project_root, source_data_root
from w54_marker_annotation import MARKERS, call_astrocytes, lineage_scores

SOURCE = source_data_root() / "Other_Datasets/Multiome_Dataset/files"
OUT = project_root() / "results/W54_Li2023_motor_cortex_validation"
CHUNK_ROWS = 25_000_000


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(SOURCE / "Multiome_Dataset_RNA_features.tsv", sep="\t", header=None)[0].astype(str).tolist()
    meta = pd.read_csv(SOURCE / "Multiome_Dataset_Metadata.txt", sep="\t", usecols=["CellId", "nCount_RNA", "WNN_L2"])
    marker_genes = sorted({g for gs in MARKERS.values() for g in gs})
    row_to_marker = {genes.index(g) + 1: i for i, g in enumerate(marker_genes)}
    lookup = np.full(len(genes) + 1, -1, dtype=np.int64)
    for row, i in row_to_marker.items():
        lookup[row] = i

    rows, cols, vals = [], [], []
    reader = pd.read_csv(
        SOURCE / "Multiome_Dataset_RNA_counts_raw.mtx", sep=" ", header=None, skiprows=2,
        names=["row", "col", "value"], dtype={"row": np.int64, "col": np.int64, "value": np.float64},
        chunksize=CHUNK_ROWS, engine="c",
    )
    for chunk in reader:
        m = lookup[chunk["row"].to_numpy()]
        keep = m >= 0
        rows.append(chunk["col"].to_numpy()[keep] - 1)
        cols.append(m[keep])
        vals.append(chunk["value"].to_numpy()[keep])
    marker_counts = sp.csr_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(len(meta), len(marker_genes))
    )
    # Only marker rows were read, so full-transcriptome totals come from nCount_RNA.
    totals = meta["nCount_RNA"].to_numpy(dtype=float)
    scores, lineages = lineage_scores(marker_counts, marker_genes, totals)
    called = call_astrocytes(scores, lineages, totals)
    truth = meta["WNN_L2"].eq("Astrocytes").to_numpy()
    tp = int((called & truth).sum())
    summary = pd.DataFrame([{
        "n_nuclei": len(meta), "n_wnn_astrocytes": int(truth.sum()), "n_rule_astrocytes": int(called.sum()),
        "true_positive": tp, "precision": tp / max(called.sum(), 1), "recall": tp / max(truth.sum(), 1),
    }])
    summary.to_csv(OUT / "w54_annotation_rule_check_ulm.tsv", sep="\t", index=False)
    fp = pd.Series(meta.loc[called & ~truth, "WNN_L2"]).value_counts().rename("rule_astro_but_wnn_label")
    fp.to_csv(OUT / "w54_annotation_rule_check_ulm_false_positive_labels.tsv", sep="\t")
    print(summary.round(4).to_string(index=False))
    print(fp.head(8).to_string())


if __name__ == "__main__":
    main()
