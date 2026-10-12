#!/usr/bin/env python3
"""Build W67 inputs (Amendment 10 B): donor-level pseudobulk per motor and sensory territory.

LSC: ventral horn, dorsal horn, lateral funiculus (Med_Lat_White + Vent_Lat_White) and dorsal
columns (Dors_Med_White). MTC: Layer_5 and Layer_2 + Layer_3. Raw counts summed per donor over
sections. No disease comparison is made here.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_w66_visium_pseudobulk import H5AD, discovery_donors, indicator  # noqa: E402
from p01_paths import project_root  # noqa: E402

logger = logging.getLogger(__name__)
OUT = project_root() / "results" / "W67_motor_selectivity" / "inputs"
TERRITORY = {
    "LSC": {"Vent_Horn": "ventral_horn", "Dors_Horn": "dorsal_horn", "Med_Lat_White": "lateral_funiculus",
            "Vent_Lat_White": "lateral_funiculus", "Dors_Med_White": "dorsal_columns"},
    "MTC": {"Layer_5": "layer5", "Layer_2": "layer2_3", "Layer_3": "layer2_3"},
}


def build(tissue: str) -> None:
    a = ad.read_h5ad(H5AD[tissue])
    obs = a.obs.reset_index(drop=True)
    terr = obs["annotation"].astype(str).map(TERRITORY[tissue])
    labels = (obs["Subject ID"].astype(str) + "|" + terr).where(terr.notna())
    ind, names = indicator(labels.astype(object).where(labels.notna(), None))
    pb = (ind @ a.X.tocsr()).T.tocsr()
    meta = pd.DataFrame({"unit": names, "donor": [n.split("|")[0] for n in names],
                         "territory": [n.split("|")[1] for n in names], "n_spots": np.asarray(ind.sum(axis=1)).ravel()})
    dm = obs[["Subject ID", "Level 1", "Diagnosis", "Sex", "Age", "PMI"]].drop_duplicates("Subject ID")
    dm.columns = ["donor", "group", "diagnosis", "sex", "age", "pmi"]
    meta = meta.merge(dm, on="donor", how="left")
    meta["gse330130_donor"] = meta["donor"].isin(discovery_donors())
    meta["lib_size"] = np.asarray(pb.sum(axis=0)).ravel()
    counts = pd.DataFrame(pb.toarray(), index=a.var_names, columns=meta["unit"])
    counts.to_csv(OUT / f"{tissue}_territory_counts.tsv.gz", sep="\t")
    meta.to_csv(OUT / f"{tissue}_territory_units.tsv", sep="\t", index=False)
    logger.info("%s: %d units", tissue, len(meta))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    for tissue in ("LSC", "MTC"):
        build(tissue)


if __name__ == "__main__":
    main()
