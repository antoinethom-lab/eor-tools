"""
polymer_injectivity
===================
Analytical polymer injectivity toolkit — eppok EURL.

Public API
----------
from polymer_injectivity import make_model, normalize_config, export_csv

Models
------
- vertical_delamaide      : 3-zone radial shear-thinning (SPE-195513)
- vertical_uvim           : Full UVIM viscoelastic, no fracture (Abdullah 2023)
- vertical_uvim_fracture  : UVIM + PKN/KGD fracture propagation (SPE-215083)
- horizontal_aitkulov     : Radial-to-linear horizontal well (Aitkulov 2021)
"""

from .engine import (
    make_model,
    normalize_config,
    export_csv,
    make_plot,
    load_rate_schedule_from_excel,
    rate_schedule_summary,
    BASE_DEFAULTS,
    MODEL_DEFAULTS,
    DEMO_CASES,
    VerticalDelamaideModel,
    UVIMVerticalModel,
    VerticalUVIMFractureModel,
    HorizontalAitkulovModel,
)

__all__ = [
    "make_model",
    "normalize_config",
    "export_csv",
    "make_plot",
    "load_rate_schedule_from_excel",
    "rate_schedule_summary",
    "BASE_DEFAULTS",
    "MODEL_DEFAULTS",
    "DEMO_CASES",
    "VerticalDelamaideModel",
    "UVIMVerticalModel",
    "VerticalUVIMFractureModel",
    "HorizontalAitkulovModel",
]

__version__ = "1.0.0"
__author__ = "Antoine — eppok EURL"
