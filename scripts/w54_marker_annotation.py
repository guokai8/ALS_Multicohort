"""Marker-dominance broad cell-type rule shared by W54 builders.

The rule was fixed before any W54 disease testing. Lineage score = mean
log1p(counts per 10k) over that lineage's markers. A nucleus is called an
astrocyte when the astrocyte score is the highest score and exceeds the
second-highest by at least MARGIN. No final W34 module gene is used as a marker.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

MARKERS: dict[str, list[str]] = {
    "Astrocyte": ["AQP4", "GJA1", "SLC1A3", "SLC1A2", "ALDH1L1", "SLC4A4"],
    "Neuron": ["SNAP25", "RBFOX3", "SYT1", "GAD1", "SLC17A7", "NRGN"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOBP", "MOG", "ST18"],
    "OPC": ["PDGFRA", "VCAN", "CSPG4"],
    "Microglia": ["CSF1R", "P2RY12", "CX3CR1", "PTPRC"],
    "Vascular": ["CLDN5", "PECAM1", "VWF", "PDGFRB"],
}
MARGIN = 0.5
MIN_UMI = 500


def lineage_scores(
    counts_cells_by_genes: sp.csr_matrix, gene_names: list[str], totals: np.ndarray | None = None
) -> tuple[np.ndarray, list[str]]:
    """Return a cells x lineages score matrix and the lineage order.

    `totals` are per-cell full-transcriptome UMIs; pass them when the matrix
    holds marker genes only.
    """
    index = {g: i for i, g in enumerate(gene_names)}
    if totals is None:
        totals = np.asarray(counts_cells_by_genes.sum(axis=1)).ravel()
    scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
    lineages = list(MARKERS)
    scores = np.zeros((counts_cells_by_genes.shape[0], len(lineages)))
    for j, lineage in enumerate(lineages):
        cols = [index[g] for g in MARKERS[lineage] if g in index]
        if not cols:
            raise ValueError(f"No markers for {lineage} found in feature list")
        sub = counts_cells_by_genes[:, cols].toarray() * scale[:, None]
        scores[:, j] = np.log1p(sub).mean(axis=1)
    return scores, lineages


def call_astrocytes(scores: np.ndarray, lineages: list[str], totals: np.ndarray) -> np.ndarray:
    astro = lineages.index("Astrocyte")
    ordered = np.sort(scores, axis=1)
    best_is_astro = scores.argmax(axis=1) == astro
    margin_ok = (ordered[:, -1] - ordered[:, -2]) >= MARGIN
    return best_is_astro & margin_ok & (totals >= MIN_UMI)
