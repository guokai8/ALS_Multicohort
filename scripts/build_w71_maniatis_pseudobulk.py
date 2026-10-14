#!/usr/bin/env python3
"""Build W71 inputs from the Maniatis et al. 2019 ALS-ST data (Amendment 14).

Mouse: raw counts summed per animal and annotated territory (spots with exactly one region),
for SOD1-G93A and SOD1-WT animals. Human: a spot table (spot totals, lineage and module
marker sums) for the descriptive co-localisation in W71 B.

Spot identifiers in the count files are `x_y` array coordinates; they are matched to the
annotation `xPos`/`yPos` columns. Arrays with two sections are split into `_1`/`_2` files.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from p01_paths import project_root, source_data_root

logger = logging.getLogger(__name__)

SRC = source_data_root() / "Maniatis2019"
COUNTS = SRC / "counts" / "count_matrices"
ANNOT = SRC / "annotations" / "annotations"
META = SRC / "metadata" / "metadata"
OUT = project_root() / "results" / "W71_maniatis_spatial_time_course" / "inputs"
REGIONS = ["Vent_Horn", "Dors_Horn", "Vent_Lat_White", "Dors_Med_White"]
GENOTYPES = {"B6SJLSOD1-G93A": "G93A", "B6SJLSOD1-WT": "WT"}
# W66 marker sets (no module gene) and the frozen ALS-up module (human symbols)
LINEAGES = {
    "astrocyte": ["AQP4", "GJA1", "SLC1A3", "SLC1A2", "ALDH1L1", "SLC4A4"],
    "microglia": ["P2RY12", "CX3CR1", "CSF1R", "PTPRC", "CD68", "TREM2", "TYROBP", "CD163", "ITGAX"],
    "endothelial": ["CLDN5", "PECAM1", "VWF", "FLT1"],
    "mural_fibroblast": ["PDGFRB", "DCN", "COL1A2", "LUM"],
    "oligodendrocyte": ["MBP", "PLP1", "MOBP", "MOG", "ST18"],
    "neuron": ["SNAP25", "RBFOX3", "SYT1", "NRGN", "SLC17A7", "GAD1"],
}
SEED = project_root() / "results/W34_synapse_resolved_regulatory_reanalysis/inputs/rna_seed_lists"


def read_counts(sample: str) -> pd.DataFrame:
    """Genes x spots raw counts; duplicate gene rows are summed."""
    f = COUNTS / f"{sample}_stdata_aligned_counts_IDs.txt"
    x = pd.read_csv(f, sep="\t", index_col=0)
    return x.groupby(level=0).sum()


def read_annotation(sample: str) -> pd.Series:
    """Region label per spot id for spots with exactly one region (others dropped)."""
    array = sample[:-2] if sample.endswith(("_1", "_2")) else sample
    a = pd.read_csv(ANNOT / f"{array}.tsv.annotations.tsv", sep="\t", dtype={"xPos": str, "yPos": str})
    labels = [c for c in a.columns if c not in ("Slide", "Array", "xPos", "yPos")]
    one = a[labels].sum(axis=1) == 1
    a = a[one]
    return pd.Series(a[labels].idxmax(axis=1).to_numpy(), index=a["xPos"] + "_" + a["yPos"])


def mouse() -> None:
    meta = pd.read_csv(META / "mouse_sample_names_sra.tsv", sep="\t", dtype=str)
    meta = meta[meta["breed"].isin(GENOTYPES)].copy()
    meta["genotype"] = meta["breed"].map(GENOTYPES)
    sums: dict[str, pd.Series] = {}
    rows, spot_log = [], []
    for _, m in meta.iterrows():
        if not (COUNTS / f"{m.sample_name}_stdata_aligned_counts_IDs.txt").exists():
            spot_log.append({"sample": m.sample_name, "status": "no count file"})
            continue
        x, lab = read_counts(m.sample_name), read_annotation(m.sample_name)
        shared = x.columns.intersection(lab.index)
        spot_log.append({"sample": m.sample_name, "status": "ok", "spots": x.shape[1], "annotated": len(shared)})
        for region in REGIONS:
            cols = shared[lab[shared] == region]
            if len(cols) == 0:
                continue
            unit = f"{m.isolate}__{region}"
            s = x[cols].sum(axis=1)
            sums[unit] = sums[unit].add(s, fill_value=0) if unit in sums else s
            rows.append({"unit": unit, "animal": m.isolate, "genotype": m.genotype, "age": m.age, "sex": m.sex,
                         "region": region, "n_spots": len(cols)})
    counts = pd.DataFrame(sums).fillna(0).astype(np.int64)
    units = (pd.DataFrame(rows).groupby(["unit", "animal", "genotype", "age", "sex", "region"], as_index=False)
             .agg(n_spots=("n_spots", "sum"), n_sections=("n_spots", "size")))
    units["lib_size"] = units["unit"].map(counts.sum(axis=0))
    units = units.set_index("unit").loc[counts.columns].rename_axis("unit").reset_index()
    counts.to_csv(OUT / "w71_mouse_unit_counts.tsv.gz", sep="\t", compression="gzip")
    units.to_csv(OUT / "w71_mouse_units.tsv", sep="\t", index=False)
    pd.DataFrame(spot_log).to_csv(OUT / "w71_mouse_section_log.tsv", sep="\t", index=False)
    logger.info("mouse: %d units, %d animals, %d genes", counts.shape[1], units["animal"].nunique(), counts.shape[0])


def human() -> None:
    meta = pd.read_csv(META / "human_sample_names_sra.tsv", sep="\t", dtype=str)
    up = pd.read_csv(SEED / "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv", sep="\t")["common_name"].tolist()
    sets = {**LINEAGES, "module_up": up}
    rows = []
    for _, m in meta.iterrows():
        if not (COUNTS / f"{m.sample_name}_stdata_aligned_counts_IDs.txt").exists():
            continue
        x, lab = read_counts(m.sample_name), read_annotation(m.sample_name)
        d = pd.DataFrame({"spot": x.columns, "total": x.sum(axis=0).to_numpy()})
        d["region"] = lab.reindex(x.columns).to_numpy()
        for name, genes in sets.items():
            g = x.index.intersection(genes)
            d[f"sum_{name}"] = x.loc[g].sum(axis=0).to_numpy()
        d.insert(0, "section", m.sample_name)
        d.insert(1, "patient", m["patient #ID"].split("_")[0])
        d.insert(2, "site", m["sample site"])
        rows.append(d)
    spots = pd.concat(rows, ignore_index=True)
    spots.to_csv(OUT / "w71_human_spots.tsv.gz", sep="\t", index=False, compression="gzip")
    logger.info("human: %d spots, %d sections, %d patients", len(spots), spots["section"].nunique(), spots["patient"].nunique())


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    mouse()
    human()


if __name__ == "__main__":
    main()
