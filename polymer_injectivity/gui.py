#!/usr/bin/env python3
"""
Polymer Injectivity Calculator — GUI (Final Audited Version)
═════════════════════════════════════════════════════════════

Wraps the audited polymer_injectivity_final_bhp_target.py engine in a tkinter GUI.

Adds support for pressure-targeted fracture growth:
  - criterion: original fracture-mechanics target
  - bhp_target: grow fracture to try to keep BHP near a user target

MODELS:
  1. Delamaide/Lake (SPE-195513): 3-zone radial shear-thinning, no fracture
  2. UVIM (Abdullah 2023): Full viscoelastic (shear-thin + thickening), no fracture
  3. UVIM + Fracture (SPE-215083): UVIM coupled to PKN/KGD fracture
  4. Horizontal (Aitkulov 2021): Radial-to-linear flow for horizontal wells

PRESETS:
  Belayim (Egypt), Matzen (Austria), Pelican Lake (Canada), Wara (Kuwait)

Usage:
  python polymer_injectivity_gui_final.py

Requirements:
  pip install matplotlib numpy
  polymer_injectivity_final_bhp_target.py should be in the same directory
  (falls back to polymer_injectivity_final.py if the new engine is absent)
"""

import math
import sys
import os
import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext, filedialog
import json

# Ensure engine module is importable from same directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
except ImportError:
    print("ERROR: matplotlib required. pip install matplotlib")
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("ERROR: numpy required. pip install numpy")
    sys.exit(1)

try:
    from polymer_injectivity.engine import (
        make_model, normalize_config, export_csv, BASE_DEFAULTS, MODEL_DEFAULTS,
        DEMO_CASES, ALIASES, FT2M, PSI2PA, IN2M, CP2PAS, MD2M2,
    )
    ENGINE_NAME = "polymer_injectivity.engine"
except ImportError:
    # Fallback: run gui.py directly from the polymer_injectivity/ directory
    from engine import (  # type: ignore
        make_model, normalize_config, export_csv, BASE_DEFAULTS, MODEL_DEFAULTS,
        DEMO_CASES, ALIASES, FT2M, PSI2PA, IN2M, CP2PAS, MD2M2,
    )
    ENGINE_NAME = "engine"

try:
    from polymer_injectivity.presets import (
        PARAM_GROUPS, INT_PARAMS, PRESETS, MODEL_LABELS,
    )
except ImportError:
    from presets import (  # type: ignore
        PARAM_GROUPS, INT_PARAMS, PRESETS, MODEL_LABELS,
    )


# PARAM_GROUPS, INT_PARAMS, PRESETS, MODEL_LABELS now live in presets.py
# (imported above) so the desktop GUI and the web app share one source of truth.


EQUATIONS_TEXT = """
══════════════════════════════════════════════════════════════════════════════════
             POLYMER INJECTIVITY CALCULATOR — THEORY, EQUATIONS & WORKFLOW
══════════════════════════════════════════════════════════════════════════════════

OVERVIEW
────────
This tool implements four analytical models for predicting polymer injectivity in
vertical and horizontal injection wells, based on peer-reviewed publications.
Each model computes bottom-hole pressure (BHP), injectivity index (II), polymer
front radius, and optionally fracture half-length as a function of time.

All models assume:
  • Steady-state (or quasi-steady-state) radial or linear flow
  • Incompressible fluids; single active layer
  • Constant injection rate (or a user-defined step schedule)
  • Reservoir pre-waterflooded to residual oil saturation (tertiary mode, default)
    OR polymer injected from day one into connate-water reservoir (secondary mode)

HOW TO CHOOSE A MODEL
──────────────────────
The four models differ in three key dimensions:
  (1) Rheology: power-law only vs. full viscoelastic (shear-thinning + thickening)
  (2) Fracture: matrix-only vs. coupled fracture propagation
  (3) Well geometry: vertical radial vs. horizontal radial-to-linear

  ┌─────────────────────────────┬──────────────────┬───────────────┬────────────────┐
  │ Model                       │ Rheology          │ Fracture?     │ Well geometry  │
  ├─────────────────────────────┼──────────────────┼───────────────┼────────────────┤
  │ 1. Delamaide/Lake           │ Power-law (3-zone)│ No            │ Vertical       │
  │ 2. UVIM                     │ Full viscoelastic │ No            │ Vertical       │
  │ 3. UVIM + Fracture          │ Full viscoelastic │ Yes (PKN/KGD) │ Vertical       │
  │ 4. Horizontal (Aitkulov)    │ Power-law (3-zone)│ No            │ Horizontal     │
  └─────────────────────────────┴──────────────────┴───────────────┴────────────────┘

► Use Model 1 (Delamaide) when:
    - You only have power-law viscometer data (K, n) and field injectivity data
    - BHP is well below fracture pressure (matrix injection confirmed)
    - Fast screening or history matching of field cases
    - BHP will plateau once shear-thinning zone fills in; the model has no built-in
      fracture relief, so BHP keeps rising if you exceed parting pressure

► Use Model 2 (UVIM, no fracture) when:
    - Your polymer exhibits measurable shear-thickening (high-MW HPAM, γ > 10 s⁻¹)
    - You want to assess how close BHP gets to the fracture initiation pressure
    - You need to see the full near-wellbore viscosity profile
    ⚠ WARNING: If UVIM predicts BHP > σ_hmin, fractures are expected. Use Model 3.
       The UVIM-only result becomes non-physical above the parting pressure.

► Use Model 3 (UVIM + Fracture) when:
    - BHP exceeds or approaches the minimum horizontal stress (σ_hmin)
    - You want to predict fracture half-length evolution and width
    - You are designing injection pressure limits or assessing fracture containment
    - You need to compare PKN vs. KGD fracture geometry assumptions
    - The fracture provides pressure relief: BHP flattens or decreases once fracture
      grows, unlike Models 1 and 2 where BHP rises without bound above σ_hmin

► Use Model 4 (Horizontal / Aitkulov) when:
    - The injector is a horizontal well (e.g., heavy oil polymer flood)
    - Flow transitions from radial (near wellbore) to linear (far field)
    - Low near-wellbore velocities make shear-thickening negligible
    ⚠ NOTE: This model has NO fracture mechanism. BHP rises linearly once the
      polymer front enters the linear flow regime and will not plateau. This is
      physically correct for matrix-only injection; in practice, either reduce rate
      or add a BHP limit before field application.

WHAT CHANGES WHEN YOU SWITCH MODELS
──────────────────────────────────────
Switching from Delamaide → UVIM:
  • The rheology changes from a simple 3-zone power-law to a full Carreau +
    elongational viscosity model. You gain two more zones (shear-thickening near
    the wellbore) and four new parameters (μ_max, λ₂, n₂, τ).
  • BHP will typically be HIGHER than Delamaide for the same rate, because the
    UVIM captures the viscosity upturn at high shear near the wellbore.
  • No change to the fracture logic (both are fracture-off in this comparison).

Switching from UVIM → UVIM + Fracture:
  • A fracture initiation criterion is now active: once BHP ≥ p_fi (≈ σ_hmin),
    a fracture opens and propagates.
  • The fracture is represented by an equivalent wellbore radius r_we (Prats 1961),
    which replaces r_w in the UVIM integral. Larger r_we → lower near-wellbore
    pressure drop → BHP relief.
  • BHP behaviour changes dramatically: instead of rising without limit, BHP
    tends to plateau or even decrease slightly as the fracture grows.
  • New active parameters: σ_hmin, E, ν, K_Ic, C_turb, h_f, fracture model (PKN/KGD).
  • The fracture growth mode (see below) further controls how pressure evolves.

Switching from Vertical → Horizontal (Aitkulov):
  • Flow geometry changes: near-wellbore is radial (with L replacing h in the
    flow equations), then transitions to 1-D linear flow at r = h/π.
  • The polymer front position x advances linearly with time once in linear flow,
    so the pressure drop in that zone also grows linearly → straight-line BHP rise.
  • No shear-thickening term (power-law only). Outer zone uses a separate viscosity
    and relative permeability for the waterflooded region ahead of the front.
  • No fracture; BHP is unbounded and will keep rising at constant injection rate.

FRACTURE GROWTH MODES (Model 3 only)
──────────────────────────────────────
Two modes control how the fracture half-length L_f evolves at each timestep:

  criterion (default):
    • Fracture grows when BHP ≥ p_fi + K_Ic / √(π·A_f)   (fracture-mechanics criterion)
    • L_f increases just enough to bring this criterion into balance
    • Physically: fracture opens when the stress intensity factor at the tip
      exceeds toughness K_Ic. This is the geomechanically rigorous approach.
    • BHP is determined by the fracture geometry — not pre-set by the user.
    • If K_Ic = 0, fracture grows as soon as BHP ≥ σ_hmin (zero-toughness limit).
    • Use this mode when you want the model to determine BHP from physics.
    • Key inputs: σ_hmin, K_Ic, E, ν, h_f.

  bhp_target:
    • L_f is adjusted iteratively at each timestep to keep BHP ≈ a user-specified
      target pressure (fracture_bhp_target_psi).
    • Physically: this simulates a field strategy of growing the fracture just enough
      to maintain injection at a fixed wellhead/BHP limit (e.g., regulatory cap or
      tubing pressure limit).
    • BHP is controlled by the user target — not determined from crack-tip mechanics.
    • If the target is below the pressure needed for matrix injection, the fracture
      cannot relieve enough to reach target; L_f grows to its maximum (fracture_max_length_ft).
    • Use this mode when you have a fixed maximum operating BHP and want to know
      what fracture size is needed to maintain a given injection rate.
    • Key input: fracture_bhp_target_psi (in addition to fracture geometry params).

  fracture_persistence (checkbox):
    • When checked: once L_f reaches a value, it never shrinks in subsequent steps
      (irreversible fracture — physically realistic for rock).
    • When unchecked: L_f can shrink if pressure drops (elastic closure — conservative).
    • Recommended: leave checked for field applications.


──────────────────────────────────────────────────────────────────────────────────
MODEL 1 — DELAMAIDE/LAKE 3-ZONE RADIAL (SPE-195513, Lake 2010, Green & Willhite 1998)
──────────────────────────────────────────────────────────────────────────────────
PHYSICS:
  Polymer solution flows radially outward from the injector through three zones,
  each characterized by a different viscosity behavior:

    Zone 1 — Upper Newtonian (near wellbore, r_w to r_uc):
      Shear rate γ > γ_uc. Viscosity ≈ μ_w (water viscosity).
      Dominant only at very high injection rates or low-permeability rocks.

    Zone 2 — Shear-Thinning (r_uc to r_lc):
      Intermediate shear rate. Viscosity follows a power-law:
        μ_pl(γ) = K · γ^(n−1)          [cp, s⁻¹]
      This is typically the dominant zone for most field cases.

    Zone 3 — Lower Newtonian (r_lc to r_p):
      Shear rate γ < γ_lc. Viscosity ≈ μ₀ (low-shear plateau).
      Occurs far from wellbore; important in long horizontal wells.

CRITICAL SHEAR RATES (transition between zones):
  γ_uc = (μ_w  / K)^(1/(n−1))   [s⁻¹]  — upper Newtonian boundary
  γ_lc = (μ₀   / K)^(1/(n−1))   [s⁻¹]  — lower Newtonian boundary

SHEAR RATE IN POROUS MEDIA (Cannella et al. 1988):
  γ_eff = C_corr · [(3n+1)/(4n)]^(n/(n−1)) · 4u / √(8k·krw·φ·Sw)   [s⁻¹]

  where:
    u       = Darcy velocity [m/s]
    k       = permeability [m²]
    krw     = water endpoint relative permeability [—]
    φ       = porosity [—]
    Sw      = water saturation [—]
    C_corr  = correction factor (≈1.0–6.0; tune to coreflood data)

  Note: C_corr = 4 (Cannella) ≈ 3× higher than Lake simplified form.
        Use C_corr ≈ 1.0–1.3 when matching Delamaide field cases.

CRITICAL RADII (from injection rate q and critical velocities):
  r_uc = q / (2π · h · u_uc)   [m]
  r_lc = q / (2π · h · u_lc)   [m]

  where u_uc, u_lc are Darcy velocities at γ_uc, γ_lc respectively.

PRESSURE DROP EQUATIONS (SPE-195513, Eqs. 4–7):
  ΔP₁ = (q / 2πh·λ_uN) · ln(r_uc / r_w)                                   [Pa]
  ΔP₂ = (Rk·H_pl / k·krw·(1−n)) · (q/2πh)ⁿ · (r_lc^(1−n) − r_uc^(1−n))  [Pa]
  ΔP₃ = (q / 2πh·λ_lN) · ln(r_p / r_lc)                                    [Pa]
  ΔP_w= (q / 2πh·λ_w)  · ln(r_e / r_p)   — waterflooded outer zone         [Pa]
  BHP  = P_res + ΔP₁ + ΔP₂ + ΔP₃ + ΔP_w                                   [psi]

  Polymer mobilities:
    λ_uN = k·krw / (Rk·μ_w)    [m²/Pa·s]   — upper Newtonian
    λ_lN = k·krw / (Rk·μ₀)    [m²/Pa·s]   — lower Newtonian
    H_pl = K · (4·α)^(n−1) · (√(8k·krw·φ·Sw))^(1−n)   — power-law integral coeff.
    where α = C_corr · [(3n+1)/(4n)]^(n/(n−1))

  Injectivity Index:
    II = q / (BHP − P_res)    [bbl/d/psi]

POLYMER BANK RADIUS (material balance, Green & Willhite 1998 Eq. 8):
  r_p = √[ Wi·ρ_w·C / (π·h·(ρ_w·C·φ_eff·S_wf + (1−φ)·ρ_r·D)) ]   [m]

  where:
    Wi      = cumulative injected polymer solution volume [m³]
    C       = polymer concentration [kg/m³]  (= ppm/1000)
    φ_eff   = φ − φ_IPV (effective porosity, accounting for inaccessible pore volume)
    S_wf    = 1 − Sorw (water saturation behind front)
    ρ_w     = brine density [kg/m³]
    ρ_r     = rock grain density [kg/m³]
    D       = adsorption [kg_polymer/kg_rock]  (= μg/g × 10⁻⁶)
    IPV     = inaccessible pore volume fraction [—]

PARAMETERS NEEDED (Model 1):
  Reservoir: k [md], h [m], φ [—], Sw [—], Sorw [—], krwmax [—], re [ft], rw [m], P_res [psi]
  Polymer:   μ_w [cp], μ₀ [cp], K [cp·sⁿ⁻¹], n [—], C_p [ppm], D [μg/g], IPV [—], Rk [—], C_corr [—]
  Injection: q [bbl/d], duration [days]


──────────────────────────────────────────────────────────────────────────────────
MODEL 2 — UVIM: UNIFIED VISCOELASTIC MODEL (Abdullah et al. 2023, Geoenergy Sci Eng 221:111259)
──────────────────────────────────────────────────────────────────────────────────
PHYSICS:
  Extends Model 1 to account for shear-thickening (elongational viscoelasticity)
  near the wellbore, where polymer chains cannot relax fast enough and viscosity
  increases with shear rate. This is critical for HPAM polymers at near-wellbore
  velocities, and can dramatically increase BHP beyond the parting pressure.

APPARENT VISCOSITY (UVM, Delshad et al. 2008, Eq. 4):
  μ_app(γ) = μ∞
            + (μ₀ᵖ − μ_w) · [1 + (λ·γ)²]^((n₁−1)/2)        ← shear-thinning term (Carreau)
            + μ_max · [1 − exp(−(λ₂·τ·γ)^(n₂−1))]           ← shear-thickening term

  Parameters:
    μ∞        = viscosity at high shear rate ≈ μ_w   [cp]
    μ₀ᵖ       = zero-shear polymer viscosity          [cp]
    μ_max     = maximum elongational viscosity         [cp]
    λ         = shear-thinning time constant           [s]
    λ₂        = shear-thickening parameter (≈ 0.01)   [—]
    τ         = in-situ viscoelastic relaxation time   [s]
    n₁        = Carreau shear-thinning index           [—]
    n₂        = shear-thickening exponent              [—]
    γ = γ_eff = effective shear rate in porous media   [s⁻¹]

PRESSURE DROP INTEGRAL (Eq. 12, Abdullah et al. 2023):
  ΔP_p = ∫[r_w → r_p]  (q · μ_app(γ_eff(r))) / (2π · r · h · k_p)  dr   [Pa]

  This integral is evaluated numerically (Runge-Kutta) or via the closed-form
  asymptotic solutions using Gauss hypergeometric functions and incomplete gamma
  functions (see Appendices A–C in Abdullah et al. 2023).

  k_p = k · krw / Rk   (polymer permeability)   [m²]

TOTAL BHP:
  ΔP_T = ΔP_p + ΔP_w   (polymer zone + waterflooded outer zone)
  BHP  = P_e + ΔP_T

FRACTURE INITIATION INDICATOR:
  If BHP > parting pressure (≈ σ_hmin), fractures are expected to initiate.
  The UVIM-only model (no fracture) will predict unrealistically high BHP;
  use Model 3 to couple fracture propagation.

ADDITIONAL UVIM PARAMETERS (beyond Model 1):
  μ∞ [cp], μ₀ᵖ [cp], μ_max [cp], λ [s], λ₂ [—], n₁ [—], n₂ [—], τ [s]


──────────────────────────────────────────────────────────────────────────────────
MODEL 3 — UVIM + FRACTURE (SPE-215083, Abdullah et al. 2023 ATCE)
──────────────────────────────────────────────────────────────────────────────────
PHYSICS:
  Once BHP exceeds fracture initiation pressure (p_fi ≈ σ_hmin), a vertical
  fracture propagates from the wellbore. The fracture is represented using
  classical 2-D fracture geometry models (PKN or KGD). The fracture is
  incorporated into the UVIM via an equivalent wellbore radius (Prats 1961),
  which effectively reduces the near-wellbore pressure drop and models
  injectivity enhancement due to fracture opening.

FRACTURE INITIATION (Eq. 10, Haimson & Fairhurst 1967):
  p_fi = (3σ_hmin − σ_hmax − 2η·ΔP_res + T₀) / (2 − 2η)   [psi]

  In practice (when σ_hmax and stress tests are unavailable):
    p_fi ≈ σ_hmin   (Eqs. 11–12)

  Fracture propagates while:
    BHP ≥ p_f + K_Ic / √(π·A_f)   — Eq. 18

PKN FRACTURE MODEL (Nordgren 1972, Perkins & Kern 1961):
  Net pressure at wellbore (Eq. 19/29):
    p_net = { (E'/2h_f)^(2n+1) · [4q·(2n+1)/(n·π·h_f)]ⁿ · 4·C_turb·L_f·K }^(1/(2n+2))

  Maximum fracture width at wellbore (Eq. 20):
    w_f,max = (2h_f/E') · p_net(max)   [m]
    w̄_f = (π/5) · w_f,max             — PKN average width

  Valid when:  2·L_f > h_f   (fracture length > fracture height)

KGD FRACTURE MODEL (Geertsma & De Klerk 1969):
  p_net at wellbore (Eq. 30/Appendix B):
    p_net = { [q·(2n+1)/(n·h_f)]ⁿ · 4·K·C_turb·L_f · [E'/(4L_f)]^(2n+1) }^(1/(2n+2))

  Maximum fracture width (Eq. 31):
    w_f,max = (4·L_f / E') · p_net(max)
    w̄_f = (π/4) · w_f,max             — KGD average width

  Valid when:  2·L_f < h_f   (fracture height > fracture length)

  In both models:  E' = E / (1 − ν²)   (plane-strain modulus)   [psi]
  Polymer in fracture: power-law μ_p = K · γ^(n_p−1)   (shear-thinning only)
  C_turb = 16/(3π) for turbulent flow, 1.0 for laminar flow

EQUIVALENT WELLBORE RADIUS (Prats 1961, Eqs. 6–9):
  r_we = L_f · r_wD(a)

  Relative capacity parameter:
    a = π·k·L_f / (2·k_f·w̄_f)

  Fracture permeability (cubic law):
    k_f = β · w̄_f²  / 12    [mD]  where β = 9.413×10¹³ (ft² to mD)

  Dimensionless wellbore radius:
    r_wD = 0.511 + (a^0.95)^0.95

  For infinitely conductive fractures (F_CD > 10, a < 0.16):
    r_we ≈ 0.5 · L_f

FRACTURE GROWTH MODES (GUI options):
  criterion  — grow L_f while BHP ≥ fracture pressure criterion (Eq. 18)
  bhp_target — grow L_f iteratively to maintain BHP ≈ user-specified target

ELASTIC DESATURATION (Qi et al. 2018, Eqs. 21–23):
  Deborah number:  N_De = γ_eff · τ
  If N_De ≥ 1:  S_orp / S_orw = 1 − 0.133 · log₁₀(N_De)
  Flux through fracture faces:  flux = q_i / (4 · h_f · L_f)   [m/s]
  Note: fracture initiation typically reduces N_De sharply → minimal desaturation.

ADDITIONAL FRACTURE PARAMETERS:
  σ_hmin [psi], E [psi], ν [—], K_Ic [psi√in], C_turb [—], h_f [m],
  max L_f [ft], solver tolerance [psi], BHP target [psi] (bhp_target mode only)


──────────────────────────────────────────────────────────────────────────────────
MODEL 4 — HORIZONTAL WELL (Aitkulov et al. 2021, J. Pet. Sci. Eng. 204:108748)
──────────────────────────────────────────────────────────────────────────────────
PHYSICS:
  Flow is radial near the wellbore (r_w to h/π), then transitions to linear flow
  (from h/π outward to the polymer front distance x). The same three flow zones
  (upper Newtonian, shear-thinning, lower Newtonian) apply in both regimes.

RADIAL ZONE (r_w to h/π):
  Same pressure drop equations as Model 1 (Eqs. 4–6), but using horizontal
  wellbore length L instead of interval thickness h.

  ΔP_radial = (Rk·H_pl / k·krw·(1−n)) · (q/2πL)ⁿ · (r_lc^(1−n) − r_uc^(1−n))

LINEAR ZONE (h/π to polymer front x):
  ΔP_linear,ST = (Rk·H_pl / k·krw) · (q / 2hL)ⁿ · (x − h/π)   [Pa]
  ΔP_linear,lN = (q / 2hL·λ_lN) · (x − h/π)

POLYMER FRONT DISTANCE x (including adsorption):
  x = [ Wi·ρ_w·C − (h/2π)²·π·L·φ·(1−Sorw)·C·ρ_w ] / [2hL·(φ·ρ_w·C·(1−Sorw) + (1−φ)·ρ_r·D) ]
                                                                + h/π   [m]

  Nine theoretical cases exist depending on which flow regimes are active
  at the wellbore (r_w) and at the polymer front (r_p or x). The GUI
  automatically selects the correct case (flowchart in Aitkulov 2021, Fig. A-1).

NOTE: Shear-thickening is negligible for horizontal wells (low near-wellbore
velocities). This model focuses on shear-thinning and lower Newtonian regimes.

ADDITIONAL PARAMETERS:
  L [m] (horizontal well length), well spacing [m]


──────────────────────────────────────────────────────────────────────────────────
UNIT REFERENCE TABLE
──────────────────────────────────────────────────────────────────────────────────
  Parameter                   GUI Unit        SI Equivalent
  ─────────────────────────── ─────────────── ───────────────────────────────────
  Permeability (k)            md              × 9.869×10⁻¹⁶ → m²
  Pay thickness (h)           m               (SI input)
  Wellbore radius (rw)        m               (SI input)
  Drainage radius (re)        ft              × 0.3048 → m
  Injection rate (q)          bbl/d           × 1.840×10⁻⁶ → m³/s
  Pressure (BHP, σ, P_res)    psi             × 6894.76 → Pa
  Viscosity (μ)               cp              × 10⁻³ → Pa·s
  Power-law K                 cp·sⁿ⁻¹         × 10⁻³ → Pa·sⁿ
  Concentration (C_p)         ppm             ÷ 1000 → kg/m³
  Adsorption (D)              μg/g            × 10⁻⁶ → kg/kg
  Young's modulus (E)         psi             × 6894.76 → Pa
  Fracture toughness (K_Ic)   psi·√in         × 1099.7 → Pa·√m
  Fracture width (w_f)        inches          × 0.0254 → m
  Fracture half-length (L_f)  ft              × 0.3048 → m
  Relaxation time (τ, λ)      s               (SI input)
  Density (ρ)                 kg/m³           (SI input)
  Skin (s)                    —               applied as r_we = r_w · e^(−s)


──────────────────────────────────────────────────────────────────────────────────
WORKFLOW GUIDE
──────────────────────────────────────────────────────────────────────────────────
STEP 1 — Select a model (left panel):
  • Start with "Delamaide/Lake" for a quick check using power-law rheology
    and field data from viscometer/coreflood.
  • Use "UVIM" when the polymer shows shear-thickening (high MW HPAM,
    elevated shear rates near wellbore, e.g. > 10 s⁻¹).
  • Use "UVIM + Fracture" for vertical well designs above parting pressure,
    or to assess fracture containment risk.
  • Use "Horizontal Well" for horizontal injectors (heavy oil, e.g. Pelican Lake).

STEP 2 — Load a preset or enter parameters:
  • Use a preset to reproduce a published field case as a sanity check.
  • Then adjust to your reservoir/polymer.
  • Key parameters to tune first: k, h, krwmax, K, n, μ₀, C_corr, Rk.

STEP 3 — Match rheology:
  • Fit K and n to lab viscosity–shear rate data (shear-thinning branch).
  • Fit λ, n₁ to complete the Carreau shear-thinning curve.
  • Fit μ_max, n₂, τ to coreflood shear-thickening data (if UVIM).
  • C_corr converts bulk viscometer shear rate to in-situ porous media shear
    rate; tune to match injectivity at known rate (typical range: 1–6).

STEP 4 — Run and interpret plots:
  • BHP plot: compare to measured wellhead + hydrostatic head.
    Flat BHP → fracturing; rising BHP → matrix flow dominates.
  • Pressure components: identify which zone drives the pressure drop.
    Shear-thinning usually dominant; lower Newtonian grows with rp.
  • Fracture plot: monitor L_f vs. inter-well distance.
    Keep L_f < 1/3 of well spacing to avoid early breakthrough.
  • Polymer front (rp): compare to inter-well distance to estimate breakthrough.

STEP 5 — Sensitivity and design:
  • Vary q: lower rate → lower BHP but longer project.
  • Vary C_p: higher concentration → higher viscosity → higher BHP.
  • Vary σ_hmin: determines when fractures initiate.
  • BHP target mode: size fracture to keep BHP below a target (safe injection).

STEP 6 — Export results:
  • Export CSV for further analysis or history matching.
  • Export PNG for reports.
  • Save JSON config to reproduce the run later.


──────────────────────────────────────────────────────────────────────────────────
KEY REFERENCES
──────────────────────────────────────────────────────────────────────────────────
  [1] Delamaide E. (2019). SPE-195513-MS. How to Use Analytical Tools to Forecast
      Injectivity in Polymer Floods. SPE Europec / 81st EAGE.
  [2] Abdullah M.B. et al. (2023). Unified Viscoelastic Injectivity Model: Analytical
      Solutions Predicting Polymer Excess Pressure and Fracture Initiation.
      Geoenergy Science and Engineering 221: 111259.
  [3] Abdullah M.B. et al. (2023). SPE-215083-MS. An Analytical Tool to Predict
      Fracture Extension and Elastic Desaturation for Polymer Field Projects.
      SPE ATCE, San Antonio, TX.
  [4] Aitkulov A. et al. (2021). An Analytical Tool to Forecast Horizontal Well
      Injectivity in Viscous Oil Polymer Floods.
      J. Petroleum Science and Engineering 204: 108748.
  [5] Lake L.W. (2010). Enhanced Oil Recovery. Society of Petroleum Engineers.
  [6] Cannella W.J., Huh C., Seright R.S. (1988). Prediction of Xanthan Rheology
      in Porous Media. SPE-18089-MS. SPE ATCE, Houston.
  [7] Delshad M. et al. (2008). Mechanistic Interpretation and Utilization of
      Viscoelastic Behavior of Polymer Solutions. SPE-113620-MS.
  [8] Qi P. et al. (2018). Simulation of Viscoelastic Polymer Flooding — From Lab
      to Field. SPE-191498-MS. SPE ATCE, Dallas.
  [9] Prats M. (1961). Effect of Vertical Fractures on Reservoir Behavior.
      SPE Journal 1(02): 105–118.
═══════════════════════════════════════════════════════════════════════════
"""


# ══════════════════════════════════════════════════════════════════════════════
#  GUI APPLICATION
# ══════════════════════════════════════════════════════════════════════════════
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Polymer Injectivity Calculator — Final Audited + BHP Target")
        self.geometry("1650x970")
        self.minsize(1300, 750)

        self.entries = {}
        self.results = None
        self.model_var = tk.StringVar(value="vertical_uvim_fracture")
        self.frac_model_var = tk.StringVar(value="PKN")
        self.frac_on_var = tk.BooleanVar(value=True)
        self.krw_wf_var = tk.BooleanVar(value=True)
        self.tertiary_var = tk.BooleanVar(value=True)
        self.persist_var = tk.BooleanVar(value=True)
        self.frac_growth_mode_var = tk.StringVar(value="criterion")

        self._load_defaults()
        self._build_ui()
        self.after(200, self._run_calc)

    def _load_defaults(self):
        self.params = dict(BASE_DEFAULTS)

    # ──────────────────────────────────────────────────────────────────
    #  UI BUILD
    # ──────────────────────────────────────────────────────────────────
    def _build_ui(self):
        nb = ttk.Notebook(self)
        nb.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        tab1 = ttk.Frame(nb); nb.add(tab1, text=" Input & Plots ")
        tab2 = ttk.Frame(nb); nb.add(tab2, text=" Equations & Theory ")
        tab3 = ttk.Frame(nb); nb.add(tab3, text=" Data Table ")

        self._build_main_tab(tab1)
        self._build_equations_tab(tab2)
        self._build_data_tab(tab3)

    def _build_main_tab(self, parent):
        pw = ttk.PanedWindow(parent, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True, padx=2, pady=2)

        # ── LEFT PANEL ──────────────────────────────────────────────
        left = ttk.Frame(pw, width=480)
        pw.add(left, weight=0)

        # Title
        ttk.Label(left, text="Polymer Injectivity Calculator",
                  font=("Helvetica", 13, "bold")).pack(pady=(6,0))
        ttk.Label(left, text=f"Engine: {ENGINE_NAME} — 4 models, paper-verified",
                  font=("Helvetica", 8, "italic")).pack(pady=(0,4))

        # Model selection
        mf = ttk.LabelFrame(left, text="Model", padding=4)
        mf.pack(fill=tk.X, padx=4, pady=2)
        for key, label in MODEL_LABELS.items():
            ttk.Radiobutton(mf, text=label, variable=self.model_var,
                            value=key).pack(anchor="w", pady=1)

        # Preset selector
        pf = ttk.Frame(left)
        pf.pack(fill=tk.X, padx=4, pady=2)
        ttk.Label(pf, text="Preset:").pack(side=tk.LEFT, padx=(0,4))
        self.preset_combo = ttk.Combobox(pf, values=list(PRESETS.keys()),
                                         state="readonly", width=32)
        self.preset_combo.pack(side=tk.LEFT, padx=2)
        ttk.Button(pf, text="Load", command=self._load_preset).pack(side=tk.LEFT, padx=2)

        # Scrollable parameters
        canvas = tk.Canvas(left, width=455)
        sb = ttk.Scrollbar(left, orient=tk.VERTICAL, command=canvas.yview)
        sf = ttk.Frame(canvas)
        sf.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0,0), window=sf, anchor="nw")
        canvas.configure(yscrollcommand=sb.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        def _wheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind_all("<MouseWheel>", _wheel)

        for group_name, params in PARAM_GROUPS:
            grp = ttk.LabelFrame(sf, text=group_name, padding=4)
            grp.pack(fill=tk.X, padx=3, pady=2)
            for label, key, unit, tip in params:
                rf = ttk.Frame(grp)
                rf.pack(fill=tk.X, pady=1)
                lbl = ttk.Label(rf, text=label, width=22, anchor="w", font=("Helvetica",9))
                lbl.pack(side=tk.LEFT, padx=(2,3))
                default_val = self.params.get(key, "")
                if default_val is None:
                    default_val = ""
                var = tk.StringVar(value=str(default_val))
                ent = ttk.Entry(rf, textvariable=var, width=11, justify="center",
                                font=("Helvetica",10,"bold"))
                ent.pack(side=tk.LEFT, padx=2)
                ttk.Label(rf, text=unit, width=7, font=("Helvetica",8,"italic"),
                          foreground="#666").pack(side=tk.LEFT)
                self.entries[key] = var

        # Rate Schedule
        rsf = ttk.LabelFrame(sf, text="Rate Schedule (optional — overrides Rate q)", padding=4)
        rsf.pack(fill=tk.X, padx=3, pady=2)

        # --- Excel loader row ---
        xl_row = ttk.Frame(rsf)
        xl_row.pack(fill=tk.X, pady=(0, 4))
        ttk.Button(xl_row, text="Load from Excel...",
                   command=self._load_rate_schedule_excel).pack(side=tk.LEFT, padx=(0, 6))
        ttk.Button(xl_row, text="Download template",
                   command=self._save_rate_template).pack(side=tk.LEFT)

        # Status label shows file name + summary after a successful load
        self.excel_status_var = tk.StringVar(value="")
        self.excel_status_label = ttk.Label(rsf, textvariable=self.excel_status_var,
                                            font=("Helvetica", 8, "italic"),
                                            foreground="#1F4E79", wraplength=340,
                                            justify="left")
        self.excel_status_label.pack(anchor="w", padx=2, pady=(0, 2))

        # Column-name configuration
        col_cfg_frame = ttk.Frame(rsf)
        col_cfg_frame.pack(fill=tk.X)
        ttk.Label(col_cfg_frame, text="Date column:",
                  font=("Helvetica", 8), foreground="#555").pack(side=tk.LEFT)
        self.xl_date_col_var = tk.StringVar(value="date")
        ttk.Entry(col_cfg_frame, textvariable=self.xl_date_col_var,
                  width=10, font=("Courier", 9)).pack(side=tk.LEFT, padx=(2, 10))
        ttk.Label(col_cfg_frame, text="Rate column:",
                  font=("Helvetica", 8), foreground="#555").pack(side=tk.LEFT)
        self.xl_rate_col_var = tk.StringVar(value="q_bpd")
        ttk.Entry(col_cfg_frame, textvariable=self.xl_rate_col_var,
                  width=10, font=("Courier", 9)).pack(side=tk.LEFT, padx=2)

        # Divider
        ttk.Separator(rsf, orient="horizontal").pack(fill=tk.X, pady=4)

        # Manual entry (fallback)
        ttk.Label(rsf, text="Or type day, rate pairs (one per line):",
                  font=("Helvetica", 8, "italic"), foreground="#666").pack(anchor="w")
        ttk.Label(rsf, text="e.g.:  0, 800\n       30, 1000\n       90, 1200",
                  font=("Courier", 8), foreground="#888").pack(anchor="w", padx=10)
        self.rate_schedule_text = tk.Text(rsf, height=5, width=30,
                                          font=("Courier", 9), bg="#FFFFF0",
                                          relief="sunken", bd=1)
        self.rate_schedule_text.pack(fill=tk.X, padx=4, pady=2)

        btn_row = ttk.Frame(rsf)
        btn_row.pack(fill=tk.X, padx=4, pady=(0, 2))
        ttk.Button(btn_row, text="Clear schedule",
                   command=self._clear_rate_schedule).pack(side=tk.LEFT)

        # Internal storage for the Excel-loaded schedule (overrides text box when set)
        self._excel_schedule = []

        # Options
        of = ttk.LabelFrame(sf, text="Options", padding=4)
        of.pack(fill=tk.X, padx=3, pady=2)

        row = ttk.Frame(of); row.pack(fill=tk.X, pady=1)
        ttk.Label(row, text="Fracture model:").pack(side=tk.LEFT)
        ttk.Radiobutton(row, text="PKN", variable=self.frac_model_var, value="PKN").pack(side=tk.LEFT, padx=4)
        ttk.Radiobutton(row, text="KGD", variable=self.frac_model_var, value="KGD").pack(side=tk.LEFT)

        row2 = ttk.Frame(of); row2.pack(fill=tk.X, pady=1)
        ttk.Label(row2, text="Growth mode:").pack(side=tk.LEFT)
        ttk.Radiobutton(row2, text="Criterion", variable=self.frac_growth_mode_var, value="criterion").pack(side=tk.LEFT, padx=4)
        ttk.Radiobutton(row2, text="BHP target", variable=self.frac_growth_mode_var, value="bhp_target").pack(side=tk.LEFT)

        ttk.Checkbutton(of, text="Enable fracture", variable=self.frac_on_var).pack(anchor="w")
        ttk.Checkbutton(of, text="Fracture persistence (monotonic growth)",
                        variable=self.persist_var).pack(anchor="w")
        ttk.Checkbutton(of, text="Use krw in waterflood zone (else kro)",
                        variable=self.krw_wf_var).pack(anchor="w")
        ttk.Checkbutton(of, text="Tertiary flood (post-waterflood Sw in shear rate/UVIM B); uncheck for secondary flood",
                        variable=self.tertiary_var).pack(anchor="w")

        # Run button
        bf = ttk.Frame(sf); bf.pack(fill=tk.X, padx=4, pady=6)
        btn = tk.Button(bf, text="▶  CALCULATE", font=("Helvetica",12,"bold"),
                        bg="#1F4E79", fg="white", activebackground="#2E6DA4",
                        activeforeground="white", relief="raised", bd=2,
                        command=self._run_calc, height=2)
        btn.pack(fill=tk.X, padx=8)

        # Export
        ef = ttk.Frame(sf); ef.pack(fill=tk.X, padx=4, pady=2)
        ttk.Button(ef, text="Export CSV", command=self._export_csv).pack(side=tk.LEFT, padx=4)
        ttk.Button(ef, text="Export PNG", command=self._export_png).pack(side=tk.LEFT, padx=4)
        ttk.Button(ef, text="Save Config JSON", command=self._save_json).pack(side=tk.LEFT, padx=4)

        # Status
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(sf, textvariable=self.status_var,
                  font=("Helvetica",9,"italic"), foreground="#333").pack(pady=3)

        # Summary
        self.summary_frame = ttk.LabelFrame(sf, text="Key Results", padding=4)
        self.summary_frame.pack(fill=tk.X, padx=3, pady=2)
        self.summary_text = tk.Text(self.summary_frame, height=12, width=55,
                                    font=("Courier",8), state="disabled",
                                    bg="#F5F5F0", relief="sunken", bd=1)
        self.summary_text.pack(fill=tk.X)

        # ── RIGHT PANEL: PLOTS ──────────────────────────────────────
        right = ttk.Frame(pw)
        pw.add(right, weight=1)
        self.fig = Figure(figsize=(10,8), dpi=100, facecolor="#FAFAFA")
        self.canvas_plot = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas_plot.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(self.canvas_plot, right)
        toolbar.update()

    def _build_equations_tab(self, parent):
        txt = scrolledtext.ScrolledText(parent, font=("Courier",9),
                                        wrap=tk.WORD, state="normal")
        txt.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        txt.insert("1.0", EQUATIONS_TEXT)
        txt.config(state="disabled")

    def _build_data_tab(self, parent):
        self.data_text = scrolledtext.ScrolledText(parent, font=("Courier",8),
                                                    wrap=tk.NONE, state="disabled")
        self.data_text.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

    # ──────────────────────────────────────────────────────────────────
    #  PARAMETER READ / PRESET
    # ──────────────────────────────────────────────────────────────────
    def _read_params(self):
        p = {}
        for key, var in self.entries.items():
            raw = var.get().strip()
            if raw == "" or raw.lower() == "none":
                continue
            try:
                if key in INT_PARAMS:
                    p[key] = int(float(raw))
                else:
                    p[key] = float(raw)
            except ValueError:
                messagebox.showerror("Input Error", f"Invalid value for '{key}': {raw}")
                return None

        p["model"] = self.model_var.get()
        p["fracture_model"] = self.frac_model_var.get()
        p["fracture_on"] = self.frac_on_var.get()
        p["fracture_persistence"] = self.persist_var.get()
        p["fracture_growth_mode"] = self.frac_growth_mode_var.get()
        p["use_krw_waterflood"] = self.krw_wf_var.get()
        p["tertiary_flood"] = self.tertiary_var.get()

        # Parse rate schedule — Excel schedule takes priority over manual text
        if self._excel_schedule:
            p["rate_schedule"] = self._excel_schedule
        else:
            raw_sched = self.rate_schedule_text.get("1.0", tk.END).strip()
            if raw_sched:
                schedule = []
                for line in raw_sched.splitlines():
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.replace(";", ",").split(",")
                    if len(parts) >= 2:
                        try:
                            day_val = float(parts[0].strip())
                            q_val = float(parts[1].strip())
                            schedule.append((day_val, q_val))
                        except ValueError:
                            messagebox.showerror("Rate Schedule Error",
                                                 f"Invalid rate schedule line: '{line}'\n"
                                                 f"Expected format: day, rate_bpd")
                            return None
                if schedule:
                    p["rate_schedule"] = schedule

        return p

    # ──────────────────────────────────────────────────────────────────
    #  RATE SCHEDULE — EXCEL LOADING
    # ──────────────────────────────────────────────────────────────────

    def _clear_rate_schedule(self):
        """Clear both the Excel-loaded schedule and the manual text box."""
        self._excel_schedule = []
        self.rate_schedule_text.delete("1.0", tk.END)
        self.excel_status_var.set("")

    def _load_rate_schedule_excel(self):
        """Open a file picker, load the Excel rate schedule, validate, and store."""
        try:
            import pandas as pd
        except ImportError:
            messagebox.showerror(
                "Missing dependency",
                "pandas is required to load Excel files.\n\n"
                "Install it with:\n    pip install pandas openpyxl\n\n"
                "Then restart the application."
            )
            return

        path = filedialog.askopenfilename(
            title="Select injection rate schedule",
            filetypes=[
                ("Excel files", "*.xlsx *.xls"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return  # user cancelled

        date_col = self.xl_date_col_var.get().strip() or "date"
        rate_col = self.xl_rate_col_var.get().strip() or "q_bpd"

        try:
            from polymer_injectivity.engine import load_rate_schedule_from_excel, rate_schedule_summary
        except ImportError:
            from engine import load_rate_schedule_from_excel, rate_schedule_summary  # type: ignore

        try:
            schedule = load_rate_schedule_from_excel(
                path,
                date_col=date_col,
                rate_col=rate_col,
            )
        except Exception as e:
            messagebox.showerror("Excel load error", str(e))
            return

        if not schedule:
            messagebox.showwarning("Empty schedule", "No valid rows were found in the file.")
            return

        self._excel_schedule = schedule

        # Show a preview in the text box so the user can see what was loaded
        self.rate_schedule_text.delete("1.0", tk.END)
        preview_lines = []
        for day, q in schedule[:20]:
            preview_lines.append(f"{day:8.1f}, {q:8.1f}")
        if len(schedule) > 20:
            preview_lines.append(f"... ({len(schedule) - 20} more rows)")
        self.rate_schedule_text.insert("1.0", "\n".join(preview_lines))
        self.rate_schedule_text.config(state="disabled", bg="#E8F4E8")

        # Update duration to match the schedule span
        days = [s[0] for s in schedule]
        total_duration = days[-1] - days[0]
        if total_duration > 0 and "duration_d" in self.entries:
            self.entries["duration_d"].set(f"{total_duration:.0f}")

        import os
        fname = os.path.basename(path)
        summary = rate_schedule_summary(schedule)
        self.excel_status_var.set(f"Loaded: {fname}\n{summary}")

    def _save_rate_template(self):
        """Save a minimal Excel template that the user can fill in."""
        try:
            import openpyxl
        except ImportError:
            messagebox.showerror(
                "Missing dependency",
                "openpyxl is required to create Excel files.\n\n"
                "Install it with:\n    pip install openpyxl"
            )
            return

        path = filedialog.asksaveasfilename(
            title="Save rate schedule template",
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx")],
            initialfile="injection_rate_schedule_template.xlsx",
        )
        if not path:
            return

        import openpyxl
        from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
        import datetime

        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = "Rate Schedule"

        # ── Header styling ──
        header_fill = PatternFill("solid", fgColor="1F4E79")
        header_font = Font(bold=True, color="FFFFFF", size=11)
        thin = Side(style="thin", color="CCCCCC")
        border = Border(left=thin, right=thin, top=thin, bottom=thin)

        headers = ["date", "q_bpd", "notes"]
        col_widths = [16, 12, 30]
        for col_idx, (h, w) in enumerate(zip(headers, col_widths), start=1):
            cell = ws.cell(row=1, column=col_idx, value=h)
            cell.font = header_font
            cell.fill = header_fill
            cell.alignment = Alignment(horizontal="center")
            cell.border = border
            ws.column_dimensions[cell.column_letter].width = w

        # ── Example rows ──
        note_fill = PatternFill("solid", fgColor="FFF9E6")
        example_rows = [
            (datetime.date(2024, 1, 1),  800,  "Injection start — ramp-up"),
            (datetime.date(2024, 2, 1), 1000,  "Step increase"),
            (datetime.date(2024, 4, 1), 1200,  "Full rate"),
            (datetime.date(2024, 7, 1), 1500,  "Plateau"),
            (datetime.date(2025, 1, 1), 1200,  "Rate reduction"),
        ]
        for row_idx, (dt, q, note) in enumerate(example_rows, start=2):
            ws.cell(row=row_idx, column=1, value=dt).number_format = "YYYY-MM-DD"
            ws.cell(row=row_idx, column=2, value=q)
            note_cell = ws.cell(row=row_idx, column=3, value=note)
            note_cell.fill = note_fill
            for col_idx in range(1, 4):
                ws.cell(row=row_idx, column=col_idx).border = border

        # ── Instructions sheet ──
        ws2 = wb.create_sheet("Instructions")
        instructions = [
            ("Polymer Injectivity Calculator — Rate Schedule Template", None),
            ("", None),
            ("HOW TO USE THIS FILE", None),
            ("", None),
            ("1. Fill in the 'Rate Schedule' sheet with your injection history or plan.", None),
            ("2. Column 'date': any date format works (YYYY-MM-DD, DD/MM/YYYY, Jan-24, etc.).", None),
            ("3. Column 'q_bpd': injection rate in barrels per day at that date.", None),
            ("4. The rate on each row is held constant until the next date (step function).", None),
            ("5. Sort rows by date (oldest first). Extra columns (like 'notes') are ignored.", None),
            ("6. In the GUI, click 'Load from Excel' and select this file.", None),
            ("   If you renamed the columns, update 'Date column' and 'Rate column' in the GUI.", None),
            ("", None),
            ("RULES", None),
            ("- Minimum 2 rows required.", None),
            ("- Dates must be unique (no duplicate dates).", None),
            ("- Rates must be positive numbers (no text, no blanks in the rate column).", None),
            ("- You can add as many rows as needed (daily, monthly, or per-event).", None),
            ("", None),
            ("COLUMN NAMES", None),
            ("The GUI reads columns named 'date' and 'q_bpd' by default.", None),
            ("If you use different names (e.g. 'Date', 'Rate_bbl_d'), enter them in the", None),
            ("'Date column' and 'Rate column' fields in the GUI before loading.", None),
        ]
        ws2.column_dimensions["A"].width = 75
        title_font = Font(bold=True, size=12, color="1F4E79")
        section_font = Font(bold=True, size=10)
        for row_idx, (text, _) in enumerate(instructions, start=1):
            cell = ws2.cell(row=row_idx, column=1, value=text)
            if row_idx == 1:
                cell.font = title_font
            elif text.isupper() and text:
                cell.font = section_font

        wb.save(path)
        messagebox.showinfo(
            "Template saved",
            f"Template saved to:\n{path}\n\n"
            "Fill in the 'Rate Schedule' sheet and load it with 'Load from Excel'."
        )

    def _load_preset(self):
        name = self.preset_combo.get()
        if not name or name not in PRESETS:
            messagebox.showinfo("Preset", "Select a preset from the dropdown first.")
            return
        preset = PRESETS[name]
        self.model_var.set(preset.get("model", "vertical_uvim_fracture"))
        self.frac_model_var.set(preset.get("fracture_model", "PKN"))
        self.frac_on_var.set(preset.get("fracture_on", True))
        self.frac_growth_mode_var.set(preset.get("fracture_growth_mode", "criterion"))
        for key, var in self.entries.items():
            if key in preset:
                var.set(str(preset[key]))
        # Clear rate schedule on preset load
        self._clear_rate_schedule()

    # ──────────────────────────────────────────────────────────────────
    #  CALCULATION
    # ──────────────────────────────────────────────────────────────────
    def _run_calc(self):
        p = self._read_params()
        if p is None:
            return
        self.status_var.set("Calculating...")
        self.update_idletasks()

        try:
            model = make_model(p)
            rows = model.run()
            self.results = rows
            self._update_plots(rows, p)
            self._update_summary(rows, p)
            self._update_data_table(rows)

            last = rows[-1]
            self.status_var.set(
                f"Done — {len(rows)} steps | "
                f"BHP={last['bhp_psi']:.0f} psi | "
                f"II={last.get('injectivity_bpd_per_psi',0):.2f} bpd/psi"
            )
        except Exception as e:
            self.status_var.set(f"Error: {e}")
            import traceback
            messagebox.showerror("Calculation Error", traceback.format_exc())

    # ──────────────────────────────────────────────────────────────────
    #  PLOTS
    # ──────────────────────────────────────────────────────────────────
    def _update_plots(self, rows, p):
        self.fig.clear()
        days = np.array([r["day"] for r in rows])
        bhp  = np.array([r["bhp_psi"] for r in rows])
        II   = np.array([r.get("injectivity_bpd_per_psi", 0) for r in rows])
        rp   = np.array([r.get("rp_ft", 0) for r in rows])
        Lf   = np.array([r.get("Lf_ft", 0) for r in rows])
        wf   = np.array([r.get("wfmax_in", 0) for r in rows])
        q_arr = np.array([r.get("q_bpd", p.get("q_bpd", 0)) for r in rows])
        has_schedule = "rate_schedule" in p and p["rate_schedule"]
        NDe  = np.array([r.get("NDe", 0) for r in rows])

        sk = dict(linewidth=1.8)
        gk = dict(alpha=0.3, linestyle="--")

        # 1) BHP
        ax1 = self.fig.add_subplot(2, 3, 1)
        ax1.plot(days, bhp, color="#C0504D", **sk)
        shmin = p.get("sigma_hmin_psi", 0)
        target = p.get("fracture_bhp_target_psi", 0)
        if shmin > 0:
            ax1.axhline(shmin, color="grey", ls=":", lw=0.8, label=f"σhmin={shmin:.0f}")
        if p.get("fracture_growth_mode") == "bhp_target" and target and target > 0:
            ax1.axhline(target, color="#2E75B6", ls="--", lw=1.0, label=f"BHP target={target:.0f}")
        if has_schedule:
            ax1r = ax1.twinx()
            ax1r.step(days, q_arr, where="post", color="#2E75B6", alpha=0.4, lw=1.2, label="q (bpd)")
            ax1r.set_ylabel("Rate [bpd]", color="#2E75B6", fontsize=8)
            ax1r.tick_params(axis="y", labelcolor="#2E75B6", labelsize=7)
        if shmin > 0 or (p.get("fracture_growth_mode") == "bhp_target" and target and target > 0):
            ax1.legend(fontsize=7)
        ax1.set(xlabel="Time [days]", ylabel="BHP [psi]", title="Bottom-Hole Pressure")
        ax1.grid(**gk)

        # 2) Pressure components
        ax2 = self.fig.add_subplot(2, 3, 2)
        dp_poly = np.array([r.get("dp_poly_psi", 0) for r in rows])
        dp_outer = np.array([r.get("dp_outer_psi", 0) for r in rows])
        # Try to get sub-components if available
        dp1 = np.array([r.get("dp1_psi", 0) for r in rows])
        dp2 = np.array([r.get("dp2_psi", 0) for r in rows])
        dp3 = np.array([r.get("dp3_psi", 0) for r in rows])
        if np.max(dp1) > 0 or np.max(dp2) > 0 or np.max(dp3) > 0:
            ax2.stackplot(days, dp1, dp2, dp3, dp_outer,
                          labels=["ΔP Upper", "ΔP Shear", "ΔP Lower", "ΔP Water"],
                          colors=["#4472C4","#70AD47","#ED7D31","#FFC000"], alpha=0.8)
        else:
            ax2.stackplot(days, dp_poly, dp_outer,
                          labels=["ΔP polymer", "ΔP water"],
                          colors=["#4472C4","#FFC000"], alpha=0.8)
        ax2.set(xlabel="Time [days]", ylabel="ΔP [psi]", title="Pressure Components")
        ax2.legend(fontsize=6, loc="upper left")
        ax2.grid(**gk)

        # 3) Injectivity
        ax3 = self.fig.add_subplot(2, 3, 3)
        ax3.plot(days, II, color="#4472C4", **sk)
        ax3.set(xlabel="Time [days]", ylabel="II [bbl/d/psi]", title="Injectivity Index")
        ax3.grid(**gk)
        if len(II) > 1 and np.max(II) > 0 and np.max(II)/max(np.min(II[II>0]),1e-6) > 10:
            ax3.set_yscale("log")

        # 4) Polymer front / radii
        ax4 = self.fig.add_subplot(2, 3, 4)
        ax4.plot(days, rp, color="#70AD47", **sk, label="rp")
        # Horizontal well x_ft
        x_ft = np.array([r.get("x_ft", 0) for r in rows])
        if np.max(x_ft) > np.max(rp) * 1.1:
            ax4.plot(days, x_ft, "--", color="#ED7D31", lw=1.2, label="x front")
        re = p.get("re_ft", 0)
        if re > 0:
            ax4.axhline(re, color="grey", ls=":", lw=0.8)
        ax4.legend(fontsize=7)
        ax4.set(xlabel="Time [days]", ylabel="Radius / front [ft]", title="Polymer Front")
        ax4.grid(**gk)

        # 5) Fracture length
        ax5 = self.fig.add_subplot(2, 3, 5)
        if np.max(Lf) > 0.01:
            ax5.plot(days, Lf, color="#ED7D31", **sk)
            ax5.set(ylabel="Lf [ft]")
        else:
            ax5.text(0.5, 0.5, "No fracture\n(BHP < σhmin or fracture off)",
                     transform=ax5.transAxes, ha="center", va="center",
                     fontsize=10, color="#999", style="italic")
        ax5.set(xlabel="Time [days]", title="Fracture Half-Length")
        ax5.grid(**gk)

        # 6) Width or Deborah
        ax6 = self.fig.add_subplot(2, 3, 6)
        if np.max(wf) > 0:
            ax6.plot(days, wf, color="#7030A0", **sk, label="wf_max (in)")
            ax6r = ax6.twinx()
            ax6r.plot(days, NDe, "--", color="#4472C4", lw=1, label="NDe")
            ax6r.set_ylabel("Deborah number", color="#4472C4")
            ax6r.legend(fontsize=6, loc="center right")
            ax6.set(ylabel="wf_max [in]")
            ax6.legend(fontsize=6, loc="upper left")
        else:
            # Show rheology curve
            gammas = np.logspace(-2, 5, 200)
            K_cp = p.get("powerlaw_K_cp", 5.5)
            n_pl = p.get("powerlaw_n", 0.8)
            mu_pl = K_cp * gammas**(n_pl - 1)
            mu_w = p.get("mu_w_cp", 1.0)
            mu0 = p.get("mu0_cp", 5.0)
            mu_eff = np.clip(mu_pl, mu_w, mu0)
            ax6.loglog(gammas, mu_eff, color="#7030A0", **sk, label="Eff μ")
            ax6.loglog(gammas, mu_pl, "--", color="#999", lw=0.8, label="Power law")
            ax6.axhline(mu_w, color="#CCC", ls=":", lw=0.8)
            ax6.axhline(mu0, color="#CCC", ls=":", lw=0.8)
            ax6.set(xlabel="γ [s⁻¹]", ylabel="μ [cp]")
            ax6.legend(fontsize=6)
        ax6.set_title("Fracture Width / Rheology", fontsize=10, fontweight="bold")
        ax6.grid(**gk)

        self.fig.tight_layout(pad=2.0)
        self.canvas_plot.draw()

    # ──────────────────────────────────────────────────────────────────
    #  SUMMARY
    # ──────────────────────────────────────────────────────────────────
    def _update_summary(self, rows, p):
        r0 = rows[0]
        rm = rows[len(rows)//2]
        re = rows[-1]
        frac_days = sum(1 for r in rows if r.get("frac") == "YES")
        max_Lf = max(r.get("Lf_ft", 0) for r in rows)
        max_wf = max(r.get("wfmax_in", 0) for r in rows)

        hdr = f"  {'':30s} {'Day 0':>10s} {'Day '+str(int(rm['day'])):>10s} {'Day '+str(int(re['day'])):>10s}"
        sep = f"  {'─'*30} {'─'*10} {'─'*10} {'─'*10}"
        lines = [
            f"  Model: {p.get('model','')}",
            f"  Fracture growth mode: {p.get('fracture_growth_mode', 'criterion')}",
            "",
            hdr, sep,
            f"  {'BHP [psi]':30s} {r0['bhp_psi']:10.1f} {rm['bhp_psi']:10.1f} {re['bhp_psi']:10.1f}",
            f"  {'II [bbl/d/psi]':30s} {r0.get('injectivity_bpd_per_psi',0):10.2f} {rm.get('injectivity_bpd_per_psi',0):10.2f} {re.get('injectivity_bpd_per_psi',0):10.2f}",
            f"  {'rp [ft]':30s} {r0.get('rp_ft',0):10.1f} {rm.get('rp_ft',0):10.1f} {re.get('rp_ft',0):10.1f}",
            f"  {'Lf [ft]':30s} {r0.get('Lf_ft',0):10.1f} {rm.get('Lf_ft',0):10.1f} {re.get('Lf_ft',0):10.1f}",
            "",
            f"  Fracture: {'YES (' + str(frac_days) + ' steps)' if frac_days else 'NO'}",
            f"  BHP target: {p.get('fracture_bhp_target_psi', '—')} psi" if p.get('fracture_growth_mode') == 'bhp_target' else "",
            f"  Max Lf: {max_Lf:.1f} ft   Max wf: {max_wf:.4f} in" if max_Lf > 0 else "",
        ]

        self.summary_text.config(state="normal")
        self.summary_text.delete("1.0", tk.END)
        self.summary_text.insert("1.0", "\n".join(lines))
        self.summary_text.config(state="disabled")

    def _update_data_table(self, rows):
        if not rows:
            return
        keys = list(rows[0].keys())
        header = "\t".join(f"{k:>14s}" for k in keys)
        lines = [header]
        for r in rows:
            vals = []
            for k in keys:
                v = r[k]
                if isinstance(v, float):
                    vals.append(f"{v:14.4f}")
                else:
                    vals.append(f"{str(v):>14s}")
            lines.append("\t".join(vals))

        self.data_text.config(state="normal")
        self.data_text.delete("1.0", tk.END)
        self.data_text.insert("1.0", "\n".join(lines))
        self.data_text.config(state="disabled")

    # ──────────────────────────────────────────────────────────────────
    #  EXPORT
    # ──────────────────────────────────────────────────────────────────
    def _export_csv(self):
        if not self.results:
            messagebox.showinfo("Export", "Run a calculation first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv",
                                            filetypes=[("CSV", "*.csv")])
        if path:
            from pathlib import Path
            export_csv(self.results, Path(path))
            messagebox.showinfo("Export", f"CSV saved to {path}")

    def _export_png(self):
        path = filedialog.asksaveasfilename(defaultextension=".png",
                                            filetypes=[("PNG", "*.png")])
        if path:
            self.fig.savefig(path, dpi=150)
            messagebox.showinfo("Export", f"Plot saved to {path}")

    def _save_json(self):
        p = self._read_params()
        if not p:
            return
        path = filedialog.asksaveasfilename(defaultextension=".json",
                                            filetypes=[("JSON", "*.json")])
        if path:
            with open(path, "w") as f:
                json.dump(p, f, indent=2)
            messagebox.showinfo("Export", f"Config saved to {path}")


if __name__ == "__main__":
    app = App()
    app.mainloop()
