#!/usr/bin/env python3
"""
Polymer injectivity analytical toolkit — final audited version.

Models
------
1) vertical_delamaide
   Vertical well, 3-zone radial Lake/Delamaide shear-thinning model.
   Ref: SPE-195513 (Delamaide 2019), Lake (2010), Green & Willhite (1998)

2) vertical_uvim
   Vertical well, UVIM/UVM radial polymer pressure drop with numerical
   integration of the full apparent viscosity (shear-thinning + thickening).
   Ref: Abdullah et al. (2023), Geoenergy Sci Eng 221:111259 — Eqs. 1–12

3) vertical_uvim_fracture
   Vertical well, UVIM polymer pressure coupled to PKN or KGD fracture
   propagation through an equivalent wellbore radius (Prats 1961).
   Ref: SPE-215083 (Abdullah et al. 2023) — PKN/KGD coupling
   Ref: SPE-224998 (AlAbdullah et al. 2026) — Wara field validation

4) horizontal_aitkulov
   Horizontal well, radial-to-linear flow model with optional outer zone.
   Ref: Aitkulov et al. (2021), J Pet Sci Eng 108748

Audit fixes applied (vs v2 and reviewed versions)
--------------------------------------------------
- Fracture geometry parameter Af: hf/4 when 2Lf>hf, else Lf (SPE-215083 §)
- Average fracture width: π/5 for PKN, π/4 for KGD (Rahman & Rahman 2010)
- Fracture persistence: recomputes full state when persisting prev_Lf
- Sw in shear rate: uses Sw_flood=1-Sorw behind polymer front
- Prats equivalent wellbore radius: proper bounded correlation
- UVIM apparent viscosity: verified against Eq.4 of Abdullah (2023)

Usage
-----
    python polymer_injectivity_final.py --model vertical_uvim_fracture --demo --csv demo.csv --plot demo.png
    python polymer_injectivity_final.py --model vertical_delamaide --print-config
    python polymer_injectivity_final.py --model horizontal_aitkulov --config case.json --csv out.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# matplotlib is optional and headless-safe
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except Exception:
    HAS_MATPLOTLIB = False


# -----------------------------------------------------------------------------
# Unit conversions
# -----------------------------------------------------------------------------
MD2M2 = 9.869233e-16
BBL2M3 = 0.158987294928
PSI2PA = 6894.757293168361
FT2M = 0.3048
IN2M = 0.0254
CP2PAS = 1e-3
BPD2M3S = BBL2M3 / 86400.0
DAY2S = 86400.0


# -----------------------------------------------------------------------------
# Default configurations and compatibility aliases
# -----------------------------------------------------------------------------
ALIASES = {
    "Pres_psi": "boundary_pressure_psi",
    "K_pl": "powerlaw_K_cp",
    "n_pl": "powerlaw_n",
    "mu_p0_cp": "mu0_cp",
    "mu_max_cp": "uvim_mu_max_cp",
    "lam1_s": "uvim_lambda_s",
    "lam2_s": "uvim_lambda2",
    "Spl": "skin",
    "KIC": "KIC_psi_sqrt_in",
    "C_turb": "Cturb",
}

BASE_DEFAULTS: Dict[str, float | int | str | bool | None | list] = dict(
    # Common simulation
    q_bpd=1000.0,
    # Rate schedule: list of (day, q_bpd) tuples.  If provided, overrides q_bpd.
    # Example: [(0, 800), (30, 1000), (90, 1200)]
    rate_schedule=None,
    duration_d=365.0,
    dt_d=1.0,
    boundary_pressure_psi=1500.0,
    skin=0.0,
    # Reservoir and fluids
    k_md=1150.0,
    h_m=5.0,
    rw_m=0.10,
    re_ft=330.0,
    phi=0.19,
    Sw=0.15,
    Sorw=0.30,
    krwmax=0.86,
    kromax=1.0,
    outer_zone_mu_cp=1.0,
    outer_zone_kr=0.86,
    mu_w_cp=1.0,
    rho_w=1020.0,
    rho_r=2650.0,
    Rk=4.4,
    IPV=0.0,
    # Polymer mass balance and power-law rheology
    Cp_ppm=1500.0,
    D_ug_g=62.9,
    powerlaw_K_cp=5.5,
    powerlaw_n=0.835,
    mu0_cp=5.0,
    C_corr=4.0,
    # UVIM/UVM
    uvim_muinf_cp=None,      # defaults to mu_w_cp
    uvim_mu0_cp=None,        # defaults to mu0_cp
    uvim_mu_max_cp=100.0,
    uvim_n1=None,            # defaults to powerlaw_n unless explicitly set
    uvim_lambda_s=0.5,
    uvim_lambda2=0.01,
    uvim_n2=1.6,
    tau_s=0.5,
    uvim_npts=500,
    # Horizontal well
    L_horiz_m=244.0,
    well_dist_m=175.0,
    outer_zone_length_m=None,  # if None, use well_dist_m/2
    # Fracture model
    fracture_on=True,
    fracture_model="PKN",
    sigma_hmin_psi=2650.0,
    fracture_initiation_pressure_psi=None,
    E_psi=1.0e6,
    nu_pr=0.25,
    hf_m=None,  # if None, use h_m
    KIC_psi_sqrt_in=1000.0,
    Cturb=16.0 / (3.0 * math.pi),
    fracture_max_length_ft=500.0,
    fracture_persistence=True,
    fracture_tol_psi=5.0,
    fracture_growth_mode="criterion",  # "criterion" or "bhp_target"
    fracture_bhp_target_psi=None,
    # Misc
    use_krw_waterflood=True,
    # tertiary_flood=True  → use Sw_flood = 1-Sorw in shear-rate / UVIM B-factor
    #                        (polymer injected after waterflood — standard assumption).
    # tertiary_flood=False → use Sw_connate = cfg["Sw"] instead
    #                        (secondary flood — polymer injected from day one).
    tertiary_flood=True,
)

MODEL_DEFAULTS = {
    "vertical_delamaide": dict(),
    "vertical_uvim": dict(
        fracture_on=False,
    ),
    "vertical_uvim_fracture": dict(),
    "horizontal_aitkulov": dict(
        fracture_on=False,
        powerlaw_n=0.80,
        powerlaw_K_cp=7.0,
        mu0_cp=20.0,
        outer_zone_mu_cp=10.0,
        outer_zone_kr=0.30,
    ),
}


# -----------------------------------------------------------------------------
# Config normalization and validation
# -----------------------------------------------------------------------------

def normalize_config(raw: Dict) -> Dict:
    cfg = dict(raw)
    for old, new in ALIASES.items():
        if old in cfg and new not in cfg:
            cfg[new] = cfg[old]

    # Compatibility with legacy well_type labels
    if cfg.get("well_type"):
        wt = str(cfg["well_type"]).strip().lower()
        if wt.startswith("horiz") and "model" not in cfg:
            cfg["model"] = "horizontal_aitkulov"
        elif wt.startswith("vert") and "model" not in cfg:
            cfg["model"] = "vertical_uvim_fracture" if cfg.get("fracture_on", True) else "vertical_delamaide"

    model = cfg.get("model", "vertical_uvim_fracture")
    if model not in MODEL_DEFAULTS:
        raise ValueError(f"Unknown model '{model}'. Valid choices: {', '.join(MODEL_DEFAULTS)}")

    merged = dict(BASE_DEFAULTS)
    merged.update(MODEL_DEFAULTS[model])
    merged.update(cfg)
    merged["model"] = model

    # Cross-defaults
    if merged["uvim_muinf_cp"] is None:
        merged["uvim_muinf_cp"] = merged["mu_w_cp"]
    if merged["uvim_mu0_cp"] is None:
        merged["uvim_mu0_cp"] = merged["mu0_cp"]
    if merged["uvim_n1"] is None:
        merged["uvim_n1"] = merged["powerlaw_n"]
    if merged["hf_m"] is None:
        merged["hf_m"] = merged["h_m"]
    if merged["fracture_initiation_pressure_psi"] is None:
        merged["fracture_initiation_pressure_psi"] = merged["sigma_hmin_psi"]
    if merged["outer_zone_length_m"] is None:
        merged["outer_zone_length_m"] = 0.5 * merged["well_dist_m"]

    # Allow old booleans/strings
    merged["fracture_model"] = str(merged["fracture_model"]).upper()
    if merged["fracture_model"] not in {"PKN", "KGD"}:
        raise ValueError("fracture_model must be 'PKN' or 'KGD'")

    merged["fracture_growth_mode"] = str(merged.get("fracture_growth_mode", "criterion")).strip().lower()
    if merged["fracture_growth_mode"] not in {"criterion", "bhp_target"}:
        raise ValueError("fracture_growth_mode must be 'criterion' or 'bhp_target'")
    if merged["fracture_growth_mode"] == "bhp_target" and merged["fracture_bhp_target_psi"] is None:
        merged["fracture_bhp_target_psi"] = merged["fracture_initiation_pressure_psi"] + 100.0

    # Basic validation
    positive_fields = [
        "q_bpd", "duration_d", "dt_d", "k_md", "h_m", "rw_m", "re_ft", "phi",
        "krwmax", "mu_w_cp", "rho_w", "rho_r", "Cp_ppm", "powerlaw_K_cp",
        "powerlaw_n", "mu0_cp", "C_corr", "uvim_muinf_cp", "uvim_mu0_cp",
        "uvim_lambda_s", "uvim_lambda2", "tau_s", "outer_zone_mu_cp",
        "outer_zone_kr", "L_horiz_m", "outer_zone_length_m", "hf_m", "E_psi",
        "Cturb", "uvim_npts"
    ]
    for name in positive_fields:
        val = merged[name]
        if val is None:
            raise ValueError(f"{name} cannot be None")
        if float(val) <= 0:
            raise ValueError(f"{name} must be > 0")

    if merged["fracture_growth_mode"] == "bhp_target":
        if merged["fracture_bhp_target_psi"] is None:
            raise ValueError("fracture_bhp_target_psi cannot be None when fracture_growth_mode='bhp_target'")
        if float(merged["fracture_bhp_target_psi"]) <= 0:
            raise ValueError("fracture_bhp_target_psi must be > 0")

    bounded = {
        "phi": (0.0, 1.0),
        "Sw": (0.0, 1.0),
        "Sorw": (0.0, 1.0),
        "IPV": (0.0, 0.99),
        "powerlaw_n": (0.01, 0.9999),
        "uvim_n1": (0.01, 0.9999),
        "nu_pr": (0.0, 0.4999),
    }
    for name, (lo, hi) in bounded.items():
        val = float(merged[name])
        if not (lo <= val <= hi):
            raise ValueError(f"{name} must be between {lo} and {hi}")

    return merged


@dataclass
class Derived:
    q_m3s: float
    k_m2: float
    re_m: float
    mu_w_pas: float
    mu0_pas: float
    K_pas_n: float
    kp_m2: float
    boundary_pressure_pa: float
    C_mass_kg_m3: float
    D_kg_kg: float
    hf_m: float
    Eprime_pa: float
    sigma_hmin_pa: float
    pfi_pa: float
    KIC_si: float


# -----------------------------------------------------------------------------
# Core utilities
# -----------------------------------------------------------------------------

def build_derived(cfg: Dict) -> Derived:
    q_m3s = cfg["q_bpd"] * BPD2M3S
    k_m2 = cfg["k_md"] * MD2M2
    kp_m2 = k_m2 * cfg["krwmax"] / cfg["Rk"]
    return Derived(
        q_m3s=q_m3s,
        k_m2=k_m2,
        re_m=cfg["re_ft"] * FT2M,
        mu_w_pas=cfg["mu_w_cp"] * CP2PAS,
        mu0_pas=cfg["mu0_cp"] * CP2PAS,
        K_pas_n=cfg["powerlaw_K_cp"] * CP2PAS,
        kp_m2=kp_m2,
        boundary_pressure_pa=cfg["boundary_pressure_psi"] * PSI2PA,
        C_mass_kg_m3=cfg["Cp_ppm"] / 1000.0,
        D_kg_kg=cfg["D_ug_g"] / 1e6,
        hf_m=cfg["hf_m"],
        Eprime_pa=cfg["E_psi"] * PSI2PA / (1.0 - cfg["nu_pr"] ** 2),
        sigma_hmin_pa=cfg["sigma_hmin_psi"] * PSI2PA,
        pfi_pa=cfg["fracture_initiation_pressure_psi"] * PSI2PA,
        KIC_si=cfg["KIC_psi_sqrt_in"] * PSI2PA * math.sqrt(IN2M),
    )


class RateSchedule:
    """Piecewise-constant injection rate schedule.

    Constructed from a list of (day, q_bpd) pairs sorted by day.
    Between entries the rate is held constant at the last specified value.
    """

    def __init__(self, entries: Optional[list] = None, default_q_bpd: float = 1000.0):
        if entries and len(entries) > 0:
            # Sort by day, ensure floats
            self.entries = sorted([(float(d), float(q)) for d, q in entries], key=lambda x: x[0])
        else:
            self.entries = [(0.0, default_q_bpd)]

    def q_bpd_at(self, day: float) -> float:
        """Return injection rate (bpd) at a given day (piecewise-constant, left-continuous)."""
        q = self.entries[0][1]
        for d, qval in self.entries:
            if day >= d:
                q = qval
            else:
                break
        return q

    def q_m3s_at(self, day: float) -> float:
        return self.q_bpd_at(day) * BPD2M3S

    def cumulative_volume_m3(self, day: float, dt_d: float) -> float:
        """Compute cumulative injected volume [m3] from day 0 to 'day' using the schedule."""
        if day <= 0:
            return 0.0
        vol = 0.0
        t = 0.0
        for i, (d_start, q_bpd) in enumerate(self.entries):
            d_end = self.entries[i + 1][0] if i + 1 < len(self.entries) else day
            d_end = min(d_end, day)
            if t >= day:
                break
            seg_start = max(t, d_start)
            if seg_start < d_end:
                vol += q_bpd * BPD2M3S * (d_end - seg_start) * DAY2S
            t = d_end
        # If day extends beyond last entry
        if t < day:
            vol += self.entries[-1][1] * BPD2M3S * (day - t) * DAY2S
        return vol


def safe_log(a: float, b: float) -> float:
    a = max(a, 1e-12)
    b = max(b, a * 1.000000001)
    return math.log(b / a)


def trapezoid_logspace_integral(func, r1: float, r2: float, npts: int = 500) -> float:
    if r2 <= r1:
        return 0.0
    radii = np.geomspace(max(r1, 1e-12), max(r2, r1 * 1.000001), int(max(npts, 50)))
    vals = func(radii)
    return float(np.trapezoid(vals, radii))


class BaseModel:
    def __init__(self, cfg: Dict):
        self.cfg = normalize_config(cfg)
        self.d = build_derived(self.cfg)
        # Build rate schedule
        raw_sched = self.cfg.get("rate_schedule")
        self.rate_schedule = RateSchedule(raw_sched, default_q_bpd=self.cfg["q_bpd"])
        # Current-step rate (updated each timestep by run())
        self._current_q_m3s = self.d.q_m3s

    def _set_timestep_rate(self, day: float) -> None:
        """Update the working injection rate for the current timestep."""
        self._current_q_m3s = self.rate_schedule.q_m3s_at(day)

    @property
    def q_m3s(self) -> float:
        """Current timestep injection rate [m3/s]."""
        return self._current_q_m3s

    # ----- Common mass-balance terms -----
    @property
    def storage_coeff(self) -> float:
        # Rigorous Green & Willhite (1998) formulation — ρ_w multiplies BOTH terms
        # in numerator and denominator, so it cancels globally.  The simplified form
        # (after cancellation) is in consistent units of [kg polymer / m³ bulk rock]:
        #
        #   storage = (φ - φ_IPV) · Sw · C  +  (1-φ) · ρ_r · D
        #
        # where:
        #   Sw = 1 - Sorw  (tertiary) or cfg["Sw"] (secondary) — honours tertiary_flood flag
        #   C  [kg/m³]   = C_mass_kg_m3  (polymer concentration in solution)
        #   D  [kg/kg]   = D_kg_kg       (adsorption, kg polymer / kg rock)
        #   ρ_r [kg/m³]  = rock grain density
        #
        # This is dimensionally consistent: both terms have units [kg polymer / m³ bulk]
        # and ρ_w does not appear (it cancels), removing the ~1000× imbalance between
        # the fluid and adsorption terms present in Delamaide's Eq. 8.
        phi_eff = self.cfg["phi"] * (1.0 - self.cfg["IPV"])
        if self.cfg.get("tertiary_flood", True):
            Sw = 1.0 - self.cfg["Sorw"]
        else:
            Sw = self.cfg["Sw"]
        return (
            phi_eff * Sw * self.d.C_mass_kg_m3
            + (1.0 - self.cfg["phi"]) * self.cfg["rho_r"] * self.d.D_kg_kg
        )

    def cumulative_injected_volume(self, day: float) -> float:
        return self.rate_schedule.cumulative_volume_m3(day, self.cfg["dt_d"])

    def polymer_bank_radius_vertical(self, Wi_m3: float) -> float:
        if Wi_m3 <= 0:
            return self.cfg["rw_m"]
        denom = math.pi * self.cfg["h_m"] * max(self.storage_coeff, 1e-30)
        # Numerator: Wi [m³] × C [kg/m³] = mass of polymer injected [kg].
        # No ρ_w here — consistent with the simplified G&W form where ρ_w cancels.
        num = Wi_m3 * self.d.C_mass_kg_m3
        rp = math.sqrt(max(num / denom, 0.0))
        return min(max(rp, self.cfg["rw_m"]), self.d.re_m)

    def water_dp_radial(self, r1_m: float, r2_m: float, mu_cp: Optional[float] = None, kr: Optional[float] = None, skin: float = 0.0) -> float:
        if r2_m <= r1_m:
            return 0.0
        mu = (self.cfg["mu_w_cp"] if mu_cp is None else mu_cp) * CP2PAS
        kr_eff = self.cfg["krwmax"] if kr is None else kr
        lam = self.d.k_m2 * kr_eff / mu
        return self.q_m3s / (2.0 * math.pi * self.cfg["h_m"] * lam) * (safe_log(r1_m, r2_m) + skin)

    def injectivity_index(self, bhp_pa: float) -> float:
        dp = bhp_pa - self.d.boundary_pressure_pa
        if dp <= 0:
            return 0.0
        q_bpd = self.q_m3s / BPD2M3S
        return q_bpd / (dp / PSI2PA)

    def _outer_zone_kr(self) -> float:
        if self.cfg.get("use_krw_waterflood", True):
            return self.cfg["krwmax"]
        return self.cfg["kromax"]


class PowerLawCommon(BaseModel):
    """Lake/Delamaide 3-zone power-law polymer injectivity model (SPE-195513)."""

    def __init__(self, cfg: Dict):
        super().__init__(cfg)
        n = self.cfg["powerlaw_n"]
        self.alpha_powerlaw = self.cfg["C_corr"] * ((1.0 + 3.0 * n) / (4.0 * n)) ** (n / (n - 1.0))
        # Sw behind the polymer front for the shear-rate denominator.
        # tertiary_flood=True  (default): polymer after waterflood → Sw = 1-Sorw
        # tertiary_flood=False           : secondary flood → Sw = connate Sw
        if self.cfg.get("tertiary_flood", True):
            Sw_flood = 1.0 - self.cfg["Sorw"]
        else:
            Sw_flood = self.cfg["Sw"]
        self.denom_powerlaw = math.sqrt(8.0 * self.d.k_m2 * self.cfg["krwmax"] * self.cfg["phi"] * Sw_flood)
        self.gamma_uc = max((self.d.mu_w_pas / self.d.K_pas_n) ** (1.0 / (n - 1.0)), 1e-9)
        self.gamma_lc = max((self.d.mu0_pas / self.d.K_pas_n) ** (1.0 / (n - 1.0)), 1e-12)
        self.u_uc = self.gamma_to_u_powerlaw(self.gamma_uc)
        self.u_lc = self.gamma_to_u_powerlaw(self.gamma_lc)
        self.lambda_uN = self.d.k_m2 * self.cfg["krwmax"] / (self.cfg["Rk"] * self.d.mu_w_pas)
        self.lambda_lN = self.d.k_m2 * self.cfg["krwmax"] / (self.cfg["Rk"] * self.d.mu0_pas)
        alpha4 = self.alpha_powerlaw * 4.0
        self.Hpl = self.d.K_pas_n * (alpha4 ** (n - 1.0)) * (self.denom_powerlaw ** (1.0 - n))

    def u_to_gamma_powerlaw(self, u_m_s: float) -> float:
        if self.denom_powerlaw <= 0:
            return 0.0
        return self.alpha_powerlaw * 4.0 * u_m_s / self.denom_powerlaw

    def gamma_to_u_powerlaw(self, gamma_s_inv: float) -> float:
        if self.alpha_powerlaw <= 0:
            return 0.0
        return gamma_s_inv * self.denom_powerlaw / (self.alpha_powerlaw * 4.0)

    def radial_critical_radii(self, q_m3s: Optional[float] = None, char_length_m: Optional[float] = None) -> Tuple[float, float]:
        q = self.q_m3s if q_m3s is None else q_m3s
        char = self.cfg["h_m"] if char_length_m is None else char_length_m
        r_uc = q / (2.0 * math.pi * max(char, 1e-12) * max(self.u_uc, 1e-30))
        r_lc = q / (2.0 * math.pi * max(char, 1e-12) * max(self.u_lc, 1e-30))
        r_uc = max(r_uc, self.cfg["rw_m"])
        r_lc = max(r_lc, r_uc)
        return r_uc, r_lc

    def dp_upper_newtonian_radial(self, r1: float, r2: float, char_length_m: float) -> float:
        if r2 <= r1:
            return 0.0
        return self.q_m3s / (2.0 * math.pi * char_length_m * self.lambda_uN) * safe_log(r1, r2)

    def dp_lower_newtonian_radial(self, r1: float, r2: float, char_length_m: float) -> float:
        if r2 <= r1:
            return 0.0
        return self.q_m3s / (2.0 * math.pi * char_length_m * self.lambda_lN) * safe_log(r1, r2)

    def dp_shear_thinning_radial(self, r1: float, r2: float, char_length_m: float) -> float:
        if r2 <= r1:
            return 0.0
        n = self.cfg["powerlaw_n"]
        coeff = self.cfg["Rk"] * self.Hpl / (self.d.k_m2 * self.cfg["krwmax"] * (1.0 - n))
        return coeff * (self.q_m3s / (2.0 * math.pi * char_length_m)) ** n * (r2 ** (1.0 - n) - r1 ** (1.0 - n))


class VerticalDelamaideModel(PowerLawCommon):
    """
    3-zone radial injectivity model per SPE-195513 (Delamaide 2019).

    Zone 1 (Upper Newtonian): rw → r_uc,  μ = μ_w            — Eq. 4
    Zone 2 (Shear-thinning):  r_uc → r_lc, μ = K·γ^(n-1)    — Eq. 5
    Zone 3 (Lower Newtonian): r_lc → rp,   μ = μ_p0          — Eq. 6
    Total: ΔP = ΔP1 + ΔP2 + ΔP3                             — Eq. 7
    Polymer bank radius from material balance                 — Eq. 8
    """

    def polymer_dp(self, rp_m: float, rw_eff_m: Optional[float] = None) -> Tuple[float, float, float, float]:
        rw_eff = self.cfg["rw_m"] if rw_eff_m is None else max(rw_eff_m, self.cfg["rw_m"] * 1e-6)
        rp = max(min(rp_m, self.d.re_m), rw_eff)
        r_uc, r_lc = self.radial_critical_radii(char_length_m=self.cfg["h_m"])

        dp1 = dp2 = dp3 = 0.0

        z1_end = min(rp, r_uc)
        if z1_end > rw_eff:
            dp1 = self.dp_upper_newtonian_radial(rw_eff, z1_end, self.cfg["h_m"])

        z2_start = max(rw_eff, r_uc)
        z2_end = min(rp, r_lc)
        if z2_end > z2_start:
            dp2 = self.dp_shear_thinning_radial(z2_start, z2_end, self.cfg["h_m"])

        z3_start = max(rw_eff, r_lc)
        z3_end = rp
        if z3_end > z3_start:
            dp3 = self.dp_lower_newtonian_radial(z3_start, z3_end, self.cfg["h_m"])

        return dp1, dp2, dp3, dp1 + dp2 + dp3

    def run(self) -> List[Dict]:
        rows: List[Dict] = []
        nsteps = int(round(self.cfg["duration_d"] / self.cfg["dt_d"])) + 1
        for i in range(nsteps):
            day = i * self.cfg["dt_d"]
            self._set_timestep_rate(day)
            Wi = self.cumulative_injected_volume(day)
            rp = self.polymer_bank_radius_vertical(Wi)
            dp1, dp2, dp3, dpp = self.polymer_dp(rp)
            dpw = self.water_dp_radial(rp, self.d.re_m, kr=self._outer_zone_kr())
            bhp = self.d.boundary_pressure_pa + dpp + dpw
            r_uc, r_lc = self.radial_critical_radii(char_length_m=self.cfg["h_m"])
            rows.append(dict(
                day=day,
                q_bpd=self.q_m3s / BPD2M3S,
                model=self.cfg["model"],
                rp_ft=rp / FT2M,
                r_uc_ft=r_uc / FT2M,
                r_lc_ft=r_lc / FT2M,
                dp_poly_psi=dpp / PSI2PA,
                dp1_psi=dp1 / PSI2PA,
                dp2_psi=dp2 / PSI2PA,
                dp3_psi=dp3 / PSI2PA,
                dp_outer_psi=dpw / PSI2PA,
                bhp_psi=bhp / PSI2PA,
                injectivity_bpd_per_psi=self.injectivity_index(bhp),
                frac="NO",
                Lf_ft=0.0,
                wfmax_in=0.0,
                wfavg_in=0.0,
                rwe_ft=self.cfg["rw_m"] / FT2M,
                NDe=0.0,
                Sorp_ratio=1.0,
            ))
        return rows


class UVIMVerticalModel(BaseModel):
    """
    UVIM (Unified Viscoelastic Injectivity Model) per Abdullah et al. (2023).

    Numerically integrates the full UVM apparent viscosity (Eq. 4 of the UVIM paper)
    including both shear-thinning (Carreau) and shear-thickening (elongational)
    contributions through the Darcy radial flow equation (Eq. 12).

    The pressure drop ΔPp is computed via trapezoidal quadrature over log-spaced
    radii from the effective wellbore radius to the polymer front.
    """

    def __init__(self, cfg: Dict):
        super().__init__(cfg)
        n1 = self.cfg["uvim_n1"]
        self.alpha_uvim = ((3.0 * n1 + 1.0) / (4.0 * n1)) ** (n1 / (n1 - 1.0)) * (4.0 / math.sqrt(8.0))
        # Per UVIM paper Eq.3 (Cannella 1988): γ_eff uses Sw in the polymer zone.
        # tertiary_flood=True  (default): Sw = max(1-Sorw, Sw) — post-waterflood
        # tertiary_flood=False           : Sw = connate Sw — secondary flood
        if self.cfg.get("tertiary_flood", True):
            Sw_eff = max(1.0 - self.cfg["Sorw"], self.cfg["Sw"])
        else:
            Sw_eff = self.cfg["Sw"]
        self.B_uvim = self.cfg["C_corr"] / math.sqrt(
            max(self.d.k_m2 * self.cfg["krwmax"] * Sw_eff * self.cfg["phi"], 1e-30)
        )

    def gamma_eff(self, u_m_s: np.ndarray | float) -> np.ndarray | float:
        return self.alpha_uvim * self.B_uvim * u_m_s

    def apparent_viscosity_cp(self, u_m_s: np.ndarray | float) -> np.ndarray | float:
        gamma = self.gamma_eff(u_m_s)
        mu_inf = self.cfg["uvim_muinf_cp"]
        mu0 = self.cfg["uvim_mu0_cp"]
        mu_w = self.cfg["mu_w_cp"]
        lam = self.cfg["uvim_lambda_s"]
        n1 = self.cfg["uvim_n1"]
        mu_max = self.cfg["uvim_mu_max_cp"]
        lam2 = self.cfg["uvim_lambda2"]
        tau = self.cfg["tau_s"]
        n2 = self.cfg["uvim_n2"]

        gamma_arr = np.asarray(gamma)
        shear_thinning = (mu0 - mu_w) * np.power(1.0 + np.power(lam * gamma_arr, 2.0), (n1 - 1.0) / 2.0)
        thick_arg = np.power(np.maximum(lam2 * tau * gamma_arr, 0.0), max(n2 - 1.0, 1e-8))
        shear_thickening = mu_max * (1.0 - np.exp(-thick_arg))
        return mu_inf + shear_thinning + shear_thickening

    def polymer_dp(self, rp_m: float, rw_eff_m: Optional[float] = None, skin: float = 0.0) -> float:
        rw_eff = self.cfg["rw_m"] if rw_eff_m is None else max(rw_eff_m, self.cfg["rw_m"] * 1e-6)
        rp = max(min(rp_m, self.d.re_m), rw_eff)
        if rp <= rw_eff:
            return 0.0

        kp = self.d.kp_m2
        q = self.q_m3s
        h = self.cfg["h_m"]
        npts = int(self.cfg["uvim_npts"])

        def integrand(r: np.ndarray) -> np.ndarray:
            ur = q / (2.0 * math.pi * h * r)
            mu_pas = np.asarray(self.apparent_viscosity_cp(ur)) * CP2PAS
            return q * mu_pas / (2.0 * math.pi * h * kp * r)

        dpp = trapezoid_logspace_integral(integrand, rw_eff, rp, npts=npts)
        if skin != 0.0:
            # Approximate additional skin using water viscosity.
            dpp += self.water_dp_radial(rw_eff, rw_eff * math.e, mu_cp=self.cfg["mu_w_cp"], kr=self.cfg["krwmax"], skin=skin)
        return dpp

    def fracture_flux_gamma(self, Lf_m: float) -> float:
        if Lf_m <= 0:
            return 0.0
        hf = max(self.d.hf_m, 1e-12)
        # same practical flux estimate used in many quick-look implementations
        flux = self.q_m3s / max(4.0 * hf * Lf_m, 1e-12)
        return float(self.gamma_eff(flux))

    def nde_and_sorp(self, Lf_m: float) -> Tuple[float, float]:
        gamma = self.fracture_flux_gamma(Lf_m)
        nde = gamma * self.cfg["tau_s"]
        if nde <= 1.0:
            return nde, 1.0
        return nde, max(1.0 - 0.133 * math.log10(max(nde, 1e-12)), 0.3)

    def run(self) -> List[Dict]:
        rows: List[Dict] = []
        nsteps = int(round(self.cfg["duration_d"] / self.cfg["dt_d"])) + 1
        for i in range(nsteps):
            day = i * self.cfg["dt_d"]
            self._set_timestep_rate(day)
            Wi = self.cumulative_injected_volume(day)
            rp = self.polymer_bank_radius_vertical(Wi)
            dpp = self.polymer_dp(rp, skin=self.cfg["skin"])
            dpw = self.water_dp_radial(rp, self.d.re_m, mu_cp=self.cfg["outer_zone_mu_cp"], kr=self._outer_zone_kr())
            bhp = self.d.boundary_pressure_pa + dpp + dpw
            rows.append(dict(
                day=day,
                q_bpd=self.q_m3s / BPD2M3S,
                model=self.cfg["model"],
                rp_ft=rp / FT2M,
                dp_poly_psi=dpp / PSI2PA,
                dp_outer_psi=dpw / PSI2PA,
                bhp_psi=bhp / PSI2PA,
                injectivity_bpd_per_psi=self.injectivity_index(bhp),
                frac="NO",
                Lf_ft=0.0,
                wfmax_in=0.0,
                wfavg_in=0.0,
                rwe_ft=self.cfg["rw_m"] / FT2M,
                NDe=0.0,
                Sorp_ratio=1.0,
            ))
        return rows


class FractureCommon(PowerLawCommon):
    """
    PKN and KGD fracture propagation models for non-Newtonian polymer.

    Ref: SPE-215083 (Abdullah et al. 2023), Appendices B–C
    PKN: Perkins-Kern-Nordgren — Eq. 19 (net pressure), Eq. 20 (width)
    KGD: Kristianovich-Geertsma-de Klerk — Eq. 30 (net pressure), Eq. 31 (width)
    Equivalent wellbore radius: Prats (1961) — Eq. 6–9
    Fracture propagation criterion — Eq. 18
    """

    def pkn_solution(self, Lf_m: float, q_m3s: Optional[float] = None) -> Tuple[float, float]:
        """PKN net pressure and max width at wellbore (SPE-215083 Eqs. 19, 20)."""
        q = self.q_m3s if q_m3s is None else q_m3s
        if Lf_m <= 0:
            return 0.0, 0.0
        hf = self.d.hf_m
        Ep = self.d.Eprime_pa
        n = self.cfg["powerlaw_n"]
        K = self.d.K_pas_n
        Ct = self.cfg["Cturb"]
        inner = (Ep / (2.0 * hf)) ** (2.0 * n + 1.0) * (((4.0 * q * (2.0 * n + 1.0) / n) / (math.pi * hf)) ** n) * (4.0 * Ct * Lf_m * K)
        pnet = inner ** (1.0 / (2.0 * n + 2.0))
        wfmax = 2.0 * hf / Ep * pnet
        return pnet, wfmax

    def kgd_solution(self, Lf_m: float, q_m3s: Optional[float] = None) -> Tuple[float, float]:
        """KGD net pressure and max width at wellbore (SPE-215083 Eqs. 30, 31)."""
        q = self.q_m3s if q_m3s is None else q_m3s
        if Lf_m <= 0:
            return 0.0, 0.0
        hf = self.d.hf_m
        Ep = self.d.Eprime_pa
        n = self.cfg["powerlaw_n"]
        K = self.d.K_pas_n
        Ct = self.cfg["Cturb"]
        inner = ((q / hf * (2.0 * n + 1.0) / n) ** n) * (4.0 * K * Ct * Lf_m) * (Ep / (4.0 * Lf_m)) ** (2.0 * n + 1.0)
        pnet = inner ** (1.0 / (2.0 * n + 2.0))
        wfmax = 4.0 * Lf_m / Ep * pnet
        return pnet, wfmax

    def fracture_geometry_parameter(self, Lf_m: float) -> float:
        """Af per SPE-215083: if 2Lf > hf then Af = hf/4 (PKN regime), else Af = Lf (KGD regime)."""
        hf = self.d.hf_m
        if 2.0 * Lf_m > hf:
            return hf / 4.0
        return Lf_m

    def fracture_width_average(self, wfmax_m: float) -> float:
        """Average fracture width: π/5 for PKN (Rahman & Rahman 2010), π/4 for KGD."""
        if self.cfg["fracture_model"] == "PKN":
            return math.pi / 5.0 * wfmax_m
        return math.pi / 4.0 * wfmax_m

    def equivalent_wellbore_radius(self, Lf_m: float, wfavg_m: float) -> float:
        """
        Prats (1961) equivalent wellbore radius (SPE-215083 Eqs. 6–9).

        kf = wf²/12 (parallel-plate permeability)
        a  = π·k·Lf / (2·kf·wf) (dimensionless relative capacity)
        rwD mapped via bounded Prats correlation:
            a < 0.16  → rwD = 0.5  (infinite conductivity)
            a < 1.0   → rwD = 0.5·exp(-π·a) (transitional)
            a ≥ 1.0   → rwD = 0.5/(π·a)  (finite conductivity)
        rwe = Lf · rwD
        """
        if Lf_m <= 0 or wfavg_m <= 0:
            return self.cfg["rw_m"]
        kf = wfavg_m ** 2 / 12.0
        a = math.pi * self.d.k_m2 * Lf_m / max(2.0 * kf * wfavg_m, 1e-30)
        a = max(a, 1e-12)
        # Practical Prats-style bounded approximation with correct asymptotes.
        if a < 0.16:
            rwD = 0.5
        elif a < 1.0:
            rwD = 0.5 * math.exp(-math.pi * a)
        else:
            rwD = 0.5 / (math.pi * a)
        return max(Lf_m * max(rwD, 1e-3), self.cfg["rw_m"])

    def fracture_state_from_length(self, Lf_m: float) -> Dict:
        if self.cfg["fracture_model"] == "PKN":
            pnet, wfmax = self.pkn_solution(Lf_m)
        else:
            pnet, wfmax = self.kgd_solution(Lf_m)
        wfavg = self.fracture_width_average(wfmax)
        rwe = self.equivalent_wellbore_radius(Lf_m, wfavg)
        Af = self.fracture_geometry_parameter(Lf_m)
        pf = self.d.sigma_hmin_pa + pnet
        return dict(Lf_m=Lf_m, pnet_pa=pnet, wfmax_m=wfmax, wfavg_m=wfavg, rwe_m=rwe, Af_m=Af, pf_pa=pf)



class VerticalUVIMFractureModel(UVIMVerticalModel, FractureCommon):
    def __init__(self, cfg: Dict):
        UVIMVerticalModel.__init__(self, cfg)
        FractureCommon.__init__(self, cfg)

    def bhp_with_rwe(self, rp_m: float, rwe_m: float) -> Tuple[float, float, float]:
        dpp = self.polymer_dp(rp_m, rw_eff_m=rwe_m, skin=self.cfg["skin"])
        dpw = self.water_dp_radial(
            rp_m,
            self.d.re_m,
            mu_cp=self.cfg["outer_zone_mu_cp"],
            kr=self._outer_zone_kr(),
        )
        bhp = self.d.boundary_pressure_pa + dpp + dpw
        return bhp, dpp, dpw

    def _criterion_target_pa(self, fr: Dict) -> float:
        toughness = self.d.KIC_si / math.sqrt(math.pi * max(fr["Af_m"], 1e-12))
        return fr["pf_pa"] + toughness

    def _active_target_pa(self, fr: Dict) -> float:
        if self.cfg["fracture_growth_mode"] == "bhp_target":
            return float(self.cfg["fracture_bhp_target_psi"]) * PSI2PA
        return self._criterion_target_pa(fr)

    def residual_for_length(self, rp_m: float, Lf_m: float) -> Tuple[float, Dict]:
        fr = self.fracture_state_from_length(Lf_m)
        bhp, dpp, dpw = self.bhp_with_rwe(rp_m, fr["rwe_m"])
        toughness = self.d.KIC_si / math.sqrt(math.pi * max(fr["Af_m"], 1e-12))
        criterion_target = fr["pf_pa"] + toughness
        active_target = self._active_target_pa(fr)
        residual = bhp - active_target
        fr.update(
            bhp_pa=bhp,
            dp_poly_pa=dpp,
            dp_outer_pa=dpw,
            toughness_pa=toughness,
            criterion_target_pa=criterion_target,
            target_pa=active_target,
            residual_pa=residual,
        )
        return residual, fr

    def _state_without_fracture(self, bhp_pa: float, dpp_pa: float, dpw_pa: float, prev_Lf_m: float = 0.0) -> Dict:
        if self.cfg["fracture_persistence"] and prev_Lf_m > 0.0:
            _, state = self.residual_for_length(self._pending_rp_m, prev_Lf_m)
            nde, sorp = self.nde_and_sorp(state["Lf_m"])
            state.update(frac="YES", NDe=nde, Sorp_ratio=sorp)
            return state

        nde, sorp = self.nde_and_sorp(0.0)
        return dict(
            frac="NO",
            Lf_m=0.0,
            wfmax_m=0.0,
            wfavg_m=0.0,
            rwe_m=self.cfg["rw_m"],
            bhp_pa=bhp_pa,
            dp_poly_pa=dpp_pa,
            dp_outer_pa=dpw_pa,
            criterion_target_pa=self.d.pfi_pa,
            target_pa=(float(self.cfg["fracture_bhp_target_psi"]) * PSI2PA
                       if self.cfg["fracture_growth_mode"] == "bhp_target"
                       else self.d.pfi_pa),
            residual_pa=bhp_pa - ((float(self.cfg["fracture_bhp_target_psi"]) * PSI2PA)
                                  if self.cfg["fracture_growth_mode"] == "bhp_target"
                                  else self.d.pfi_pa),
            NDe=nde,
            Sorp_ratio=sorp,
        )

    def solve_fracture_step(self, rp_m: float, prev_Lf_m: float = 0.0) -> Dict:
        self._pending_rp_m = rp_m
        bhp_no_frac, dpp0, dpw0 = self.bhp_with_rwe(rp_m, self.cfg["rw_m"])
        mode = self.cfg["fracture_growth_mode"]

        if mode == "bhp_target":
            target_pa = float(self.cfg["fracture_bhp_target_psi"]) * PSI2PA
            if bhp_no_frac <= target_pa:
                return self._state_without_fracture(bhp_no_frac, dpp0, dpw0, prev_Lf_m=prev_Lf_m)
        else:
            if bhp_no_frac < self.d.pfi_pa:
                return self._state_without_fracture(bhp_no_frac, dpp0, dpw0, prev_Lf_m=prev_Lf_m)

        tol = self.cfg["fracture_tol_psi"] * PSI2PA
        lo = max(1e-6, self.cfg["rw_m"])
        if self.cfg["fracture_persistence"] and prev_Lf_m > 0:
            lo = max(lo, prev_Lf_m)
        hi = max(lo, self.cfg["fracture_max_length_ft"] * FT2M)

        # Expand upper bound if needed (handles unusually high rate or soft rock).
        # Try 2×, 4×, 8× the configured max before giving up and capping at hi.
        for _exp in range(3):
            res_hi_check, _ = self.residual_for_length(rp_m, hi)
            if res_hi_check <= 0:
                break
            hi *= 2.0

        res_lo, state_lo = self.residual_for_length(rp_m, lo)
        if abs(res_lo) <= tol or res_lo <= 0:
            state = state_lo
        else:
            res_hi, state_hi = self.residual_for_length(rp_m, hi)

            if res_hi > 0:
                # BHP still exceeds target at expanded upper bound; cap at hi.
                state = state_hi
            else:
                state = state_hi
                for _ in range(80):
                    mid = 0.5 * (lo + hi)
                    res_mid, state_mid = self.residual_for_length(rp_m, mid)
                    state = state_mid
                    if abs(res_mid) <= tol:
                        break
                    if res_mid > 0:
                        lo = mid
                    else:
                        hi = mid

        if self.cfg["fracture_persistence"] and prev_Lf_m > 0 and state["Lf_m"] < prev_Lf_m:
            _, state = self.residual_for_length(rp_m, prev_Lf_m)

        nde, sorp = self.nde_and_sorp(state["Lf_m"])
        state.update(frac="YES", NDe=nde, Sorp_ratio=sorp)
        return state

    def run(self) -> List[Dict]:
        rows: List[Dict] = []
        nsteps = int(round(self.cfg["duration_d"] / self.cfg["dt_d"])) + 1
        prev_Lf = 0.0
        for i in range(nsteps):
            day = i * self.cfg["dt_d"]
            self._set_timestep_rate(day)
            Wi = self.cumulative_injected_volume(day)
            rp = self.polymer_bank_radius_vertical(Wi)
            state = self.solve_fracture_step(rp, prev_Lf_m=prev_Lf)
            prev_Lf = state["Lf_m"]
            rows.append(dict(
                day=day,
                q_bpd=self.q_m3s / BPD2M3S,
                model=self.cfg["model"],
                fracture_growth_mode=self.cfg["fracture_growth_mode"],
                rp_ft=rp / FT2M,
                dp_poly_psi=state["dp_poly_pa"] / PSI2PA,
                dp_outer_psi=state["dp_outer_pa"] / PSI2PA,
                bhp_psi=state["bhp_pa"] / PSI2PA,
                injectivity_bpd_per_psi=self.injectivity_index(state["bhp_pa"]),
                frac=state["frac"],
                Lf_ft=state["Lf_m"] / FT2M,
                wfmax_in=state["wfmax_m"] / IN2M,
                wfavg_in=state["wfavg_m"] / IN2M,
                rwe_ft=state["rwe_m"] / FT2M,
                fracture_target_psi=state.get("target_pa", 0.0) / PSI2PA,
                fracture_criterion_target_psi=state.get("criterion_target_pa", 0.0) / PSI2PA,
                fracture_residual_psi=state.get("residual_pa", 0.0) / PSI2PA,
                NDe=state["NDe"],
                Sorp_ratio=state["Sorp_ratio"],
            ))
        return rows


class HorizontalAitkulovModel(PowerLawCommon):
    def __init__(self, cfg: Dict):
        super().__init__(cfg)
        self.radial_linear_boundary = self.cfg["h_m"] / math.pi

    def polymer_front(self, Wi_m3: float) -> Dict:
        # Storage coefficient same mass-balance logic as the vertical case.
        storage = max(self.storage_coeff, 1e-30)
        if Wi_m3 <= 0:
            return dict(mode="radial", rp_m=self.cfg["rw_m"], x_m=self.cfg["rw_m"])

        radial_mass_capacity = math.pi * (self.radial_linear_boundary ** 2 - self.cfg["rw_m"] ** 2) * self.cfg["L_horiz_m"] * storage
        # Numerator mass [kg] = Wi [m³] × C [kg/m³]  (ρ_w cancels — see storage_coeff)
        injected_polymer_mass = Wi_m3 * self.d.C_mass_kg_m3

        if injected_polymer_mass <= radial_mass_capacity:
            rp = math.sqrt(injected_polymer_mass / (math.pi * self.cfg["L_horiz_m"] * storage) + self.cfg["rw_m"] ** 2)
            return dict(mode="radial", rp_m=max(rp, self.cfg["rw_m"]), x_m=max(rp, self.cfg["rw_m"]))

        extra_mass = injected_polymer_mass - radial_mass_capacity
        x = self.radial_linear_boundary + extra_mass / (2.0 * self.cfg["h_m"] * self.cfg["L_horiz_m"] * storage)
        return dict(mode="linear", rp_m=self.radial_linear_boundary, x_m=max(x, self.radial_linear_boundary))

    def _dp_linear_segment(self, dx_m: float, regime: str) -> float:
        if dx_m <= 0:
            return 0.0
        L = self.cfg["L_horiz_m"]
        h = self.cfg["h_m"]
        n = self.cfg["powerlaw_n"]
        if regime == "upper":
            return self.q_m3s / (2.0 * h * L * self.lambda_uN) * dx_m
        if regime == "lower":
            return self.q_m3s / (2.0 * h * L * self.lambda_lN) * dx_m
        coeff = self.cfg["Rk"] * self.Hpl / (self.d.k_m2 * self.cfg["krwmax"])
        return coeff * (self.q_m3s / (2.0 * L * h)) ** n * dx_m

    def _linear_regime(self) -> str:
        ux = self.q_m3s / (2.0 * self.cfg["h_m"] * self.cfg["L_horiz_m"])
        if ux >= self.u_uc:
            return "upper"
        if ux > self.u_lc:
            return "shear"
        return "lower"

    def polymer_dp(self, front: Dict) -> Tuple[float, float, float, float, float]:
        L = self.cfg["L_horiz_m"]
        rw = self.cfg["rw_m"]
        r_boundary = self.radial_linear_boundary
        rp_radial = min(front["rp_m"], r_boundary)
        r_uc, r_lc = self.radial_critical_radii(char_length_m=L)

        dp_rad_u = dp_rad_st = dp_rad_l = dp_lin = 0.0

        # Radial near wellbore region
        z1_end = min(rp_radial, r_uc)
        if z1_end > rw:
            dp_rad_u = self.dp_upper_newtonian_radial(rw, z1_end, L)

        z2_start = max(rw, r_uc)
        z2_end = min(rp_radial, r_lc)
        if z2_end > z2_start:
            dp_rad_st = self.dp_shear_thinning_radial(z2_start, z2_end, L)

        z3_start = max(rw, r_lc)
        z3_end = rp_radial
        if z3_end > z3_start:
            dp_rad_l = self.dp_lower_newtonian_radial(z3_start, z3_end, L)

        if front["mode"] == "linear":
            regime = self._linear_regime()
            dp_lin = self._dp_linear_segment(front["x_m"] - r_boundary, regime)

        total = dp_rad_u + dp_rad_st + dp_rad_l + dp_lin
        return dp_rad_u, dp_rad_st, dp_rad_l, dp_lin, total

    def outer_zone_dp(self, front: Dict) -> float:
        L = self.cfg["L_horiz_m"]
        h = self.cfg["h_m"]
        mu = self.cfg["outer_zone_mu_cp"] * CP2PAS
        kr = self.cfg["outer_zone_kr"]
        lam = self.d.k_m2 * kr / mu
        x_boundary = self.cfg["outer_zone_length_m"]
        r_boundary = self.radial_linear_boundary

        if front["mode"] == "radial":
            dp_rad = 0.0
            if r_boundary > front["rp_m"]:
                dp_rad = self.q_m3s / (2.0 * math.pi * L * lam) * safe_log(front["rp_m"], r_boundary)
            dp_lin = 0.0
            if x_boundary > r_boundary:
                dp_lin = self.q_m3s / (2.0 * h * L * lam) * (x_boundary - r_boundary)
            return dp_rad + dp_lin

        if x_boundary <= front["x_m"]:
            return 0.0
        return self.q_m3s / (2.0 * h * L * lam) * (x_boundary - front["x_m"])

    def run(self) -> List[Dict]:
        rows: List[Dict] = []
        nsteps = int(round(self.cfg["duration_d"] / self.cfg["dt_d"])) + 1
        for i in range(nsteps):
            day = i * self.cfg["dt_d"]
            self._set_timestep_rate(day)
            Wi = self.cumulative_injected_volume(day)
            front = self.polymer_front(Wi)
            dp1, dp2, dp3, dplin, dpp = self.polymer_dp(front)
            dp_outer = self.outer_zone_dp(front)
            bhp = self.d.boundary_pressure_pa + dpp + dp_outer
            rows.append(dict(
                day=day,
                q_bpd=self.q_m3s / BPD2M3S,
                model=self.cfg["model"],
                front_mode=front["mode"],
                rp_ft=front["rp_m"] / FT2M,
                x_ft=front["x_m"] / FT2M,
                radial_linear_boundary_ft=self.radial_linear_boundary / FT2M,
                dp1_psi=dp1 / PSI2PA,
                dp2_psi=dp2 / PSI2PA,
                dp3_psi=dp3 / PSI2PA,
                dp_linear_psi=dplin / PSI2PA,
                dp_poly_psi=dpp / PSI2PA,
                dp_outer_psi=dp_outer / PSI2PA,
                bhp_psi=bhp / PSI2PA,
                injectivity_bpd_per_psi=self.injectivity_index(bhp),
                frac="NO",
                Lf_ft=0.0,
                wfmax_in=0.0,
                wfavg_in=0.0,
                rwe_ft=self.cfg["rw_m"] / FT2M,
                NDe=0.0,
                Sorp_ratio=1.0,
            ))
        return rows


# -----------------------------------------------------------------------------
# Factories, export, plotting, CLI
# -----------------------------------------------------------------------------

def make_model(cfg: Dict):
    model = normalize_config(cfg)["model"]
    if model == "vertical_delamaide":
        return VerticalDelamaideModel(cfg)
    if model == "vertical_uvim":
        return UVIMVerticalModel(cfg)
    if model == "vertical_uvim_fracture":
        return VerticalUVIMFractureModel(cfg)
    if model == "horizontal_aitkulov":
        return HorizontalAitkulovModel(cfg)
    raise ValueError(f"Unhandled model '{model}'")


def export_csv(rows: List[Dict], path: Path) -> None:
    if not rows:
        raise ValueError("No rows to export")
    keys = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


# -----------------------------------------------------------------------------
# Excel rate-schedule loader
# -----------------------------------------------------------------------------

def load_rate_schedule_from_excel(
    path,
    date_col: str = "date",
    rate_col: str = "q_bpd",
    reference_date=None,
    sheet_name=0,
) -> List[Tuple[float, float]]:
    """Load a piecewise-constant injection rate schedule from an Excel file.

    Parameters
    ----------
    path : str or Path
        Path to the .xlsx / .xls file on disk.
    date_col : str
        Name of the column containing dates (default "date").
        Accepts any format pandas can parse: "2023-01-15", "15/01/2023",
        "Jan-23", Excel serial numbers, etc.
    rate_col : str
        Name of the column containing injection rates in bbl/day (default "q_bpd").
    reference_date : str, datetime, or None
        The date that maps to day 0.  If None the first date in the file is used.
    sheet_name : int or str
        Sheet to read (default: first sheet).

    Returns
    -------
    list of (day: float, q_bpd: float) tuples
        Sorted by day, ready to pass directly as ``rate_schedule=`` in a model config.

    Raises
    ------
    ImportError  – if pandas or openpyxl are not installed.
    ValueError   – if required columns are missing or the schedule is empty.

    Example
    -------
    >>> schedule = load_rate_schedule_from_excel("inj_rates.xlsx")
    >>> cfg = {"model": "vertical_uvim_fracture", "rate_schedule": schedule, ...}
    >>> rows = make_model(cfg).run()

    Expected Excel layout (minimum two columns, extra columns ignored)::

        date          q_bpd
        2023-01-01    800
        2023-02-15    1 000
        2023-04-01    1 200

    Dates can be any format Excel or pandas understands.
    The rate on each row is held constant until the next date (piecewise-constant,
    left-continuous) — exactly matching the engine's RateSchedule behaviour.
    """
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "pandas is required to load Excel rate schedules.  "
            "Install it with:  pip install pandas openpyxl"
        ) from e

    # Accept either a filesystem path (str/Path) or a file-like object
    # (e.g. a Streamlit UploadedFile / BytesIO).  Only do disk checks for paths.
    is_path_like = isinstance(path, (str, Path))
    if is_path_like:
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Excel file not found: {path}")
        display_name = path.name
    else:
        display_name = getattr(path, "name", "uploaded file")

    try:
        df = pd.read_excel(path, sheet_name=sheet_name)
    except Exception as e:
        raise ValueError(f"Could not read Excel file '{display_name}': {e}") from e

    # Normalise column names: strip whitespace, lower-case for comparison
    df.columns = [str(c).strip() for c in df.columns]
    col_map = {c.lower(): c for c in df.columns}

    def _find_col(wanted: str) -> str:
        """Return the actual column name matching 'wanted' (case-insensitive)."""
        if wanted in df.columns:
            return wanted
        if wanted.lower() in col_map:
            return col_map[wanted.lower()]
        raise ValueError(
            f"Column '{wanted}' not found in '{display_name}'.  "
            f"Available columns: {list(df.columns)}"
        )

    actual_date_col = _find_col(date_col)
    actual_rate_col = _find_col(rate_col)

    # Parse dates
    df[actual_date_col] = pd.to_datetime(df[actual_date_col], dayfirst=False, errors="coerce")
    bad = df[actual_date_col].isna().sum()
    if bad > 0:
        raise ValueError(
            f"{bad} date(s) in column '{actual_date_col}' could not be parsed.  "
            "Check date format in the Excel file."
        )

    # Parse rates — strip thousand-separators robustly.
    # Use Python's re (not pandas .str.replace) so it works regardless of the
    # pandas string backend (object vs PyArrow), which differ in escape handling.
    import re as _re
    _sep_chars = " \t,\u202f\u00a0"  # space, tab, comma, thin-space, non-breaking space
    _strip_re = _re.compile("[" + _re.escape(_sep_chars) + "]")

    def _to_float(x):
        if isinstance(x, (int, float)):
            return float(x)
        s = _strip_re.sub("", str(x))
        return float(s)

    df[actual_rate_col] = df[actual_rate_col].map(_to_float)

    df = df[[actual_date_col, actual_rate_col]].dropna().sort_values(actual_date_col).reset_index(drop=True)

    if df.empty:
        raise ValueError(f"No valid rows found in '{display_name}' after parsing.")

    # Reference date → day 0
    if reference_date is not None:
        t0 = pd.to_datetime(reference_date, dayfirst=False)
    else:
        t0 = df[actual_date_col].iloc[0]

    df["day"] = (df[actual_date_col] - t0).dt.total_seconds() / DAY2S

    # Validate: no negative days if reference_date was provided explicitly
    if reference_date is not None and (df["day"] < 0).any():
        n_neg = (df["day"] < 0).sum()
        raise ValueError(
            f"{n_neg} date(s) are before the reference date ({t0.date()}).  "
            "Check reference_date or the file contents."
        )

    schedule = [(float(row["day"]), float(row[actual_rate_col])) for _, row in df.iterrows()]
    return schedule


def rate_schedule_summary(schedule: List[Tuple[float, float]], reference_date=None) -> str:
    """Return a human-readable one-line summary of a rate schedule.

    Parameters
    ----------
    schedule : list of (day, q_bpd) tuples
    reference_date : str or datetime, optional
        If provided, convert day numbers back to calendar dates in the summary.

    Returns
    -------
    str, e.g. "18 steps · 0–365 days · 800–1 500 bpd"
    or  "18 steps · 2023-01-01 → 2023-12-31 · 800–1 500 bpd"
    """
    if not schedule:
        return "Empty schedule"
    days = [s[0] for s in schedule]
    rates = [s[1] for s in schedule]
    n = len(schedule)
    rate_str = f"{min(rates):,.0f}–{max(rates):,.0f} bpd"
    if reference_date is not None:
        try:
            import pandas as pd
            t0 = pd.to_datetime(reference_date, dayfirst=False)
            d_start = (t0 + pd.Timedelta(days=days[0])).strftime("%Y-%m-%d")
            d_end = (t0 + pd.Timedelta(days=days[-1])).strftime("%Y-%m-%d")
            return f"{n} steps · {d_start} → {d_end} · {rate_str}"
        except Exception:
            pass
    return f"{n} steps · day {days[0]:.0f}–{days[-1]:.0f} · {rate_str}"


def make_plot(rows: List[Dict], path: Path) -> None:
    if not HAS_MATPLOTLIB:
        raise RuntimeError("matplotlib is not available in this environment")
    days = [row["day"] for row in rows]
    bhp = [row["bhp_psi"] for row in rows]
    rp = [row.get("rp_ft", 0.0) for row in rows]
    ii = [row.get("injectivity_bpd_per_psi", 0.0) for row in rows]
    Lf = [row.get("Lf_ft", 0.0) for row in rows]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    ax = axes.ravel()
    ax[0].plot(days, bhp)
    ax[0].set_title("BHP")
    ax[0].set_xlabel("day")
    ax[0].set_ylabel("psi")

    ax[1].plot(days, rp)
    ax[1].set_title("Polymer front extent")
    ax[1].set_xlabel("day")
    ax[1].set_ylabel("ft")

    ax[2].plot(days, ii)
    ax[2].set_title("Injectivity index")
    ax[2].set_xlabel("day")
    ax[2].set_ylabel("bpd/psi")

    ax[3].plot(days, Lf)
    ax[3].set_title("Fracture half-length")
    ax[3].set_xlabel("day")
    ax[3].set_ylabel("ft")

    for a in ax:
        a.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def print_summary(rows: List[Dict]) -> None:
    first = rows[0]
    last = rows[-1]
    print(f"Model: {first['model']}")
    print(f"Days: {first['day']} -> {last['day']}")
    print(f"BHP start/end [psi]: {first['bhp_psi']:.1f} / {last['bhp_psi']:.1f}")
    print(f"Injectivity start/end [bpd/psi]: {first['injectivity_bpd_per_psi']:.3f} / {last['injectivity_bpd_per_psi']:.3f}")
    if "rp_ft" in last:
        print(f"Polymer front end [ft]: {last['rp_ft']:.2f}")
    if "x_ft" in last:
        print(f"Horizontal x-front end [ft]: {last['x_ft']:.2f}")
    if max(row.get("Lf_ft", 0.0) for row in rows) > 0:
        print(f"Max fracture half-length [ft]: {max(row.get('Lf_ft', 0.0) for row in rows):.2f}")


# Demo cases tuned to be numerically stable rather than history-matched.
DEMO_CASES = {
    "vertical_delamaide": dict(
        model="vertical_delamaide",
        q_bpd=1000.0,
        duration_d=180.0,
        dt_d=2.0,
    ),
    "vertical_uvim": dict(
        model="vertical_uvim",
        q_bpd=2400.0,
        duration_d=30.0,
        dt_d=1.0,
        k_md=550.0,
        h_m=16.5 * FT2M,
        re_ft=984.0,
        boundary_pressure_psi=1552.0,
        phi=0.26,
        Sw=0.87,
        krwmax=0.30,
        mu_w_cp=0.6,
        Rk=1.0,
        Cp_ppm=1000.0,
        C_corr=1.0,
        uvim_muinf_cp=0.6,
        uvim_mu0_cp=7.8,
        uvim_n1=0.4,
        uvim_lambda_s=0.65,
        uvim_mu_max_cp=33.0,
        uvim_lambda2=0.01,
        uvim_n2=2.0,
        tau_s=0.6,
        powerlaw_K_cp=7.0,
        powerlaw_n=0.85,
        mu0_cp=7.8,
        fracture_on=False,
    ),
    "vertical_uvim_fracture": dict(
        model="vertical_uvim_fracture",
        q_bpd=2400.0,
        duration_d=30.0,
        dt_d=1.0,
        k_md=550.0,
        h_m=16.5 * FT2M,
        hf_m=16.5 * FT2M,
        re_ft=984.0,
        boundary_pressure_psi=1552.0,
        phi=0.26,
        Sw=0.87,
        krwmax=0.30,
        mu_w_cp=0.6,
        Rk=1.0,
        Cp_ppm=1000.0,
        C_corr=1.0,
        uvim_muinf_cp=0.6,
        uvim_mu0_cp=7.8,
        uvim_n1=0.4,
        uvim_lambda_s=0.65,
        uvim_mu_max_cp=33.0,
        uvim_lambda2=0.01,
        uvim_n2=2.0,
        tau_s=0.6,
        powerlaw_K_cp=7.0,
        powerlaw_n=0.85,
        mu0_cp=7.8,
        sigma_hmin_psi=2650.0,
        fracture_initiation_pressure_psi=2650.0,
        E_psi=1_060_000.0,
        nu_pr=0.365,
        KIC_psi_sqrt_in=0.0,
        Cturb=1.7,
        fracture_model="PKN",
    ),
    "horizontal_aitkulov": dict(
        model="horizontal_aitkulov",
        q_bpd=1200.0,
        duration_d=180.0,
        dt_d=2.0,
        k_md=2000.0,
        h_m=10.0,
        L_horiz_m=400.0,
        well_dist_m=300.0,
        outer_zone_length_m=150.0,
        mu_w_cp=1.0,
        outer_zone_mu_cp=12.0,
        outer_zone_kr=0.3,
        krwmax=0.8,
        Rk=1.5,
        powerlaw_K_cp=15.0,
        powerlaw_n=0.8,
        mu0_cp=25.0,
        C_corr=1.5,
        phi=0.24,
        Sw=0.25,
        Sorw=0.30,
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Polymer injectivity analytical toolkit")
    parser.add_argument("--model", choices=sorted(MODEL_DEFAULTS), default="vertical_uvim_fracture")
    parser.add_argument("--config", type=str, help="Path to JSON config")
    parser.add_argument("--csv", type=str, help="Optional CSV export path")
    parser.add_argument("--plot", type=str, help="Optional PNG export path")
    parser.add_argument("--print-config", action="store_true", help="Print normalized default config and exit")
    parser.add_argument("--demo", action="store_true", help="Run a built-in demo case for the selected model")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.demo:
        cfg = dict(DEMO_CASES[args.model])
    elif args.config:
        with open(args.config, "r", encoding="utf-8") as f:
            cfg = json.load(f)
        cfg["model"] = args.model
    else:
        cfg = {"model": args.model}

    cfg = normalize_config(cfg)

    if args.print_config:
        print(json.dumps(cfg, indent=2))
        return 0

    model = make_model(cfg)
    rows = model.run()
    print_summary(rows)

    if args.csv:
        export_csv(rows, Path(args.csv))
        print(f"CSV written to {args.csv}")
    if args.plot:
        make_plot(rows, Path(args.plot))
        print(f"Plot written to {args.plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
