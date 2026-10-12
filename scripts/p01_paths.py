#!/usr/bin/env python3
"""Path helpers for the P01 reproducibility package.

Environment variables:
- P01_PROJECT_ROOT: repository/project root containing metadata/, results/ and scripts/.
- P01_SOURCE_DATA_ROOT: external source-data root containing GSE330130/, GSE290359/, etc.
- P01_PLAN_ROOT: manuscript/package planning directory.
"""

from __future__ import annotations

import os
from pathlib import Path


def _path_from_env(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    return Path(value).expanduser().resolve() if value else default.expanduser().resolve()


def project_root() -> Path:
    default = Path(__file__).resolve().parents[1]
    return _path_from_env("P01_PROJECT_ROOT", default)


def researchd_root() -> Path:
    default = project_root().parents[1]
    return _path_from_env("RESEARCHD_ROOT", default)


def source_data_root() -> Path:
    default = researchd_root() / "Data"
    return _path_from_env("P01_SOURCE_DATA_ROOT", default)


def plan_root() -> Path:
    default = (
        researchd_root()
        / "plan"
        / "ALS_project_readiness_2026-08-11"
        / "P01_ALS_CNS_cell_type_regulatory_map"
    )
    return _path_from_env("P01_PLAN_ROOT", default)
