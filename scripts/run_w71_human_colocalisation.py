#!/usr/bin/env python3
"""W71 B (Amendment 14): descriptive spot co-localisation in Maniatis et al. 2019 human ALS spinal cord.

Seven sporadic ALS patients, no controls, so no disease test. Spot scores use the W66 transform,
log1p(1e4 x summed marker counts / spot total). Within each section (>= 30 spots with >= 500 counts),
Spearman correlations between the frozen ALS-up module score and each lineage score are computed raw
and with ranks residualised on the rank of log spot depth (W66 post hoc method). Also reported: the
per-patient ventral-to-dorsal-horn ratio of module CPM (descriptive).
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from p01_paths import project_root

logger = logging.getLogger(__name__)

OUT = project_root() / "results" / "W71_maniatis_spatial_time_course"
LINEAGES = ["astrocyte", "microglia", "endothelial", "mural_fibroblast", "oligodendrocyte", "neuron"]
MIN_COUNTS, MIN_SPOTS = 500, 30


def score(d: pd.DataFrame, name: str) -> pd.Series:
    return np.log1p(1e4 * d[f"sum_{name}"] / d["total"])


def colocalisation(s: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (patient, sec), d in s.groupby(["patient", "section"]):
        if len(d) < MIN_SPOTS:
            continue
        rd = np.log(d["total"]).rank().to_numpy()

        def resid(v: np.ndarray) -> np.ndarray:
            return v - np.polyval(np.polyfit(rd, v, 1), rd)

        y = score(d, "module_up").rank().to_numpy()
        for lin in LINEAGES:
            x = score(d, lin).rank().to_numpy()
            rows.append({"patient": patient, "section": sec, "n_spots": len(d), "lineage": lin,
                         "rho": np.corrcoef(y, x)[0, 1], "rho_depth_partial": np.corrcoef(resid(y), resid(x))[0, 1]})
    return pd.DataFrame(rows)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    s = pd.read_csv(OUT / "inputs" / "w71_human_spots.tsv.gz", sep="\t")
    q = s[s["total"] >= MIN_COUNTS]
    per = colocalisation(q)
    summary = per.groupby("lineage", as_index=False).agg(
        n_sections=("section", "nunique"), median_rho=("rho", "median"), median_rho_depth_partial=("rho_depth_partial", "median"),
        frac_sections_positive=("rho_depth_partial", lambda v: float((v > 0).mean())))
    by_patient = per.groupby(["patient", "lineage"], as_index=False)["rho_depth_partial"].median()
    mg = summary.set_index("lineage")["median_rho_depth_partial"]
    others = mg[["astrocyte", "endothelial", "mural_fibroblast"]]
    reading = "replicates W66" if mg["microglia"] > others.max() else "does not replicate W66"

    # descriptive: ventral-to-dorsal-horn ratio of module CPM per patient (all spots, summed counts)
    horn = s[s["region"].isin(["Vent_Horn", "Dors_Horn"])].groupby(["patient", "region"])[["sum_module_up", "total"]].sum()
    horn["cpm"] = 1e6 * horn["sum_module_up"] / horn["total"]
    ratio = horn["cpm"].unstack("region")
    ratio["log2_ventral_over_dorsal"] = np.log2(ratio["Vent_Horn"] / ratio["Dors_Horn"])

    per.to_csv(OUT / "w71_human_colocalisation_sections.tsv", sep="\t", index=False)
    summary.assign(reading=reading).to_csv(OUT / "w71_human_colocalisation_summary.tsv", sep="\t", index=False)
    by_patient.to_csv(OUT / "w71_human_colocalisation_by_patient.tsv", sep="\t", index=False)
    ratio.reset_index().to_csv(OUT / "w71_human_ventral_dorsal_ratio.tsv", sep="\t", index=False)
    logger.info("%d qualifying spots in %d sections, %d patients", len(q), per["section"].nunique(), per["patient"].nunique())
    logger.info("\n%s\nreading: %s", summary.round(3).to_string(index=False), reading)
    logger.info("\n%s", by_patient.pivot(index="patient", columns="lineage", values="rho_depth_partial").round(2).to_string())
    logger.info("\n%s", ratio.round(2).to_string())


if __name__ == "__main__":
    main()
