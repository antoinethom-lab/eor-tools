# Engine reference — polymer_injectivity/engine.py

**Who this is for:** anyone who wants to understand, verify, or modify the physics
in `engine.py`. Every formula is linked to its source equation, every class and
method is located by line number, and every modification point is explicitly called out.

**How to use this document:** use your editor's search (`Ctrl+F`) to jump to a
keyword — model name, parameter name, SPE paper number, or formula symbol.

---

## Table of contents

1. [Architecture overview](#1-architecture-overview)
2. [Unit conversions](#2-unit-conversions)
3. [Configuration system](#3-configuration-system)
4. [Shared base classes](#4-shared-base-classes)
5. [Model 1 — Delamaide/Lake (SPE-195513)](#5-model-1--delamaide--lake-spe-195513)
6. [Model 2 — UVIM no fracture (Abdullah 2023)](#6-model-2--uvim-no-fracture-abdullah-2023)
7. [Model 3 — UVIM + Fracture (SPE-215083)](#7-model-3--uvim--fracture-spe-215083)
8. [Model 4 — Horizontal well (Aitkulov 2021)](#8-model-4--horizontal-well-aitkulov-2021)
9. [Shared fracture geometry (FractureCommon)](#9-shared-fracture-geometry-fracturecommon)
10. [Rate schedule](#10-rate-schedule)
11. [Output columns reference](#11-output-columns-reference)
12. [How to make common modifications](#12-how-to-make-common-modifications)
13. [Known limitations and assumptions](#13-known-limitations-and-assumptions)

---

## 1. Architecture overview

```
engine.py
│
├── BaseModel                        line 380
│   ├── storage_coeff                line 401   ← Green & Willhite mass balance
│   ├── polymer_bank_radius_vertical line 430   ← rp from cumulative injection
│   └── water_dp_radial              line 440   ← outer zone ΔP
│
├── PowerLawCommon(BaseModel)        line 461
│   ├── 3-zone shear-rate logic      line 464   ← Lake/Delamaide α, γ_uc, γ_lc
│   ├── dp_upper_newtonian_radial    line 504   ← Zone 1 ΔP
│   ├── dp_shear_thinning_radial     line 514   ← Zone 2 ΔP
│   └── dp_lower_newtonian_radial    line 509   ← Zone 3 ΔP
│
├── VerticalDelamaideModel           line 522   ← Model 1
│   └── inherits PowerLawCommon
│
├── UVIMVerticalModel(BaseModel)     line 593   ← Model 2
│   ├── apparent_viscosity_cp        line 623   ← UVIM Eq. 4
│   └── polymer_dp (numerical ∫)     line 641   ← UVIM Eq. 12
│
├── FractureCommon(PowerLawCommon)   line 709
│   ├── pkn_solution                 line 720   ← SPE-215083 Eqs. 19–20
│   ├── kgd_solution                 line 735   ← SPE-215083 Eqs. 30–31
│   └── equivalent_wellbore_radius   line 763   ← Prats 1961
│
├── VerticalUVIMFractureModel        line 802   ← Model 3
│   └── inherits UVIMVerticalModel + FractureCommon
│
└── HorizontalAitkulovModel          line 964   ← Model 4
    └── inherits PowerLawCommon
```

**Inheritance is the key design choice.** Models 1 and 4 share the power-law
shear-rate logic (PowerLawCommon). Models 2 and 3 share the UVIM viscosity
(UVIMVerticalModel). Model 3 additionally inherits the fracture geometry
(FractureCommon). If you change a formula in a parent class, it affects all
children — check the table above before editing.

---

## 2. Unit conversions

**Location:** `engine.py` lines 63–70

All inputs arrive in field units (md, psi, ft, cp, bbl/d). The engine converts
everything to SI immediately in `build_derived()` (line 295) and works in SI
throughout. Outputs are converted back to field units only when writing result rows.

```python
MD2M2   = 9.869233e-16      # millidarcy → m²
BBL2M3  = 0.158987          # barrel → m³
PSI2PA  = 6894.757          # psi → Pa
FT2M    = 0.3048            # foot → metre
IN2M    = 0.0254            # inch → metre
CP2PAS  = 1e-3              # centipoise → Pa·s
BPD2M3S = BBL2M3 / 86400   # bbl/day → m³/s
DAY2S   = 86400             # day → second
```

**To add a new unit:** define a constant here, use it in `build_derived()`, and
convert outputs back at the point where you append to `rows`.

---

## 3. Configuration system

### Parameter defaults — line 92

`BASE_DEFAULTS` is the master list of all parameters and their default values.
`MODEL_DEFAULTS` (line 164) overrides specific defaults per model — e.g.
`vertical_uvim` has `fracture_on=False`.

**Special `None` defaults** — these parameters are intentionally `None` and the
engine computes them from other parameters in `normalize_config()`:

| Parameter | Inherits from | Line in normalize_config |
|---|---|---|
| `uvim_muinf_cp` | `mu_w_cp` | 239 |
| `uvim_mu0_cp` | `mu0_cp` | 241 |
| `uvim_n1` | `powerlaw_n` | 243 |
| `hf_m` | `h_m` | 245 |
| `fracture_initiation_pressure_psi` | `sigma_hmin_psi` | 247 |
| `outer_zone_length_m` | `well_dist_m / 2` | 249 |
| `fracture_bhp_target_psi` | `fracture_initiation_pressure_psi + 100` | 256 |

> ⚠️ **Important for web/GUI developers:** never pass `0.0` for these fields —
> always pass `None` or omit them. `0.0` is physically invalid and the engine
> will reject it. See `NULLABLE_PARAMS` in `streamlit_app.py`.

### normalize_config() — line 185

Called at the start of every model instantiation. It:
1. Resolves aliases (old parameter names → new ones)
2. Merges BASE_DEFAULTS + MODEL_DEFAULTS + user config
3. Fills in cross-defaults (the `None` fields above)
4. Validates all values (positive, bounded, model string valid)

**To add a new parameter:** add it to `BASE_DEFAULTS` with a sensible default,
add validation logic in the `positive_fields` or `bounded` dicts (lines 260–278),
and add it to `PARAM_GROUPS` in `presets.py` so it appears in the GUI.

### build_derived() — line 295

Converts field-unit config values to SI and bundles them in a `Derived` dataclass.
Used internally by all model classes as `self.d`.

---

## 4. Shared base classes

### BaseModel — line 380

Every model inherits from this. It holds:
- `self.cfg` — normalized config dict
- `self.d` — `Derived` SI values
- `self.rate_schedule` — `RateSchedule` object

#### storage_coeff — line 401

The polymer mass-balance storage coefficient from Green & Willhite (1998).

**Formula:**
```
A = (φ · (1 - IPV) · Sw_flood · C)  +  ((1 - φ) · ρ_r · D)
```

Where:
- `φ` = porosity
- `IPV` = inaccessible pore volume fraction
- `Sw_flood` = `1 - Sorw` (tertiary, default) or `Sw` (secondary)
- `C` = polymer concentration [kg/m³]
- `ρ_r` = rock grain density [kg/m³]
- `D` = adsorption [kg polymer / kg rock]

> **Audit note:** ρ_w cancels in the full G&W derivation and is NOT present here.
> Earlier versions incorrectly included it in only one term, creating a ~1000×
> imbalance. This is the corrected form.

**To modify:** edit lines 401–427. The `tertiary_flood` flag at line 407 controls
which `Sw` is used.

#### polymer_bank_radius_vertical — line 430

Computes the radial extent of the polymer bank from cumulative injection volume.

**Formula (from material balance):**
```
rp = sqrt( Wi · C / (π · h · A) )
```
Where `Wi` = cumulative injection [m³], `A` = `storage_coeff`.

#### water_dp_radial — line 440

Radial Darcy pressure drop in the outer (water/oil) zone:
```
ΔP = (q / (2π · h · λ)) · ln(r2 / r1) + skin term
```
Where `λ = k · kr / μ` is the fluid mobility.

---

## 5. Model 1 — Delamaide / Lake (SPE-195513)

**Class:** `VerticalDelamaideModel` — line 522  
**Inherits:** `PowerLawCommon` → `BaseModel`  
**Reference:** Delamaide (2019) SPE-195513; Lake (2010) *Enhanced Oil Recovery*

### Concept

The polymer zone is divided into three radial regions based on local shear rate:

```
Wellbore ──[Zone 1]──[Zone 2]──[Zone 3]── Polymer front (rp)── Water/oil ──
           r_w     r_uc      r_lc         rp
           Upper   Shear-    Lower
           Newton  thinning  Newton
           μ≈μw    K·γ^(n-1) μ≈μ₀
```

The boundaries `r_uc` and `r_lc` are where local shear rate equals the upper
and lower critical shear rates (transition between Newtonian plateaus and
power-law zone).

### PowerLawCommon.__init__() — line 464

Pre-computes the constants that define the three zones:

| Symbol | Variable | Formula | Description |
|---|---|---|---|
| α | `alpha_powerlaw` | `C_corr · ((1+3n)/(4n))^(n/(n-1))` | Cannella (1988) shear-rate correction |
| γ_uc | `gamma_uc` | `(μ_w / K)^(1/(n-1))` | Upper critical shear rate [s⁻¹] |
| γ_lc | `gamma_lc` | `(μ₀ / K)^(1/(n-1))` | Lower critical shear rate [s⁻¹] |
| u_uc | `u_uc` | `γ_uc · √(8kφSw) / (4α)` | Darcy velocity at upper transition |
| u_lc | `u_lc` | `γ_lc · √(8kφSw) / (4α)` | Darcy velocity at lower transition |
| H | `Hpl` | `K · (4α)^(n-1) / denom^(1-n)` | Power-law ΔP coefficient |

### Zone pressure drop formulas — lines 504–521

**Zone 1 (Upper Newtonian)** — SPE-195513 Eq. 4:
```
ΔP₁ = (q / (2π h λ_uN)) · ln(r_uc / r_w)
```
where `λ_uN = k · krw / (Rk · μ_w)`

**Zone 2 (Shear-thinning)** — SPE-195513 Eq. 5:
```
ΔP₂ = (Rk · H / (k · krw · (1-n))) · (q / (2π h))^n · (r_lc^(1-n) - r_uc^(1-n))
```

**Zone 3 (Lower Newtonian)** — SPE-195513 Eq. 6:
```
ΔP₃ = (q / (2π h λ_lN)) · ln(rp / r_lc)
```
where `λ_lN = k · krw / (Rk · μ₀)`

**Total:** ΔP_poly = ΔP₁ + ΔP₂ + ΔP₃ — SPE-195513 Eq. 7

**To modify the power-law formula:** edit `dp_shear_thinning_radial()` at line 514.  
**To change the zone boundary logic:** edit `radial_critical_radii()` at line 495.  
**To change the shear-rate correction:** edit `alpha_powerlaw` at line 469 and the `C_corr` parameter.

---

## 6. Model 2 — UVIM no fracture (Abdullah 2023)

**Class:** `UVIMVerticalModel` — line 593  
**Inherits:** `BaseModel`  
**Reference:** Abdullah et al. (2023), *Geoenergy Sci. Eng.* 221:111259, Eqs. 1–12

### Concept

Instead of the three-zone approximation, the UVIM numerically integrates the
**full apparent viscosity** (including both shear-thinning AND shear-thickening)
across the polymer bank. This is more accurate for high-MW HPAM that exhibits
significant viscoelastic thickening at high shear rates near the wellbore.

### UVIMVerticalModel.__init__() — line 605

Pre-computes two constants:

**α_UVIM** — shear-rate structure factor (Cannella, adapted for UVIM):
```
α = ((3n₁+1)/(4n₁))^(n₁/(n₁-1)) · (4/√8)
```

**B** — in-situ shear-rate scale factor (Eq. 3 of Abdullah 2023):
```
B = C_corr / √(k · krw · Sw_flood · φ)
```

Where `Sw_flood = 1 - Sorw` for tertiary floods.

### apparent_viscosity_cp() — line 623

The full UVIM apparent viscosity — **Abdullah 2023, Eq. 4**:
```
μ_app = μ∞  +  (μ₀ - μ_w) · (1 + (λ₁·γ_eff)²)^((n₁-1)/2)   ← shear-thinning (Carreau)
              +  μ_max · (1 - exp(-(λ₂·τ·γ_eff)^(n₂-1)))      ← shear-thickening
```

| Symbol | Parameter | Description |
|---|---|---|
| μ∞ | `uvim_muinf_cp` | Infinite-shear viscosity ≈ μ_w |
| μ₀ | `uvim_mu0_cp` | Zero-shear polymer viscosity |
| μ_max | `uvim_mu_max_cp` | Maximum elongational viscosity |
| λ₁ | `uvim_lambda_s` | Shear-thinning time constant [s] |
| λ₂ | `uvim_lambda2` | Shear-thickening parameter |
| n₁ | `uvim_n1` | Carreau shear-thinning exponent |
| n₂ | `uvim_n2` | Shear-thickening exponent |
| τ | `tau_s` | Viscoelastic relaxation time [s] |
| γ_eff | computed | Effective shear rate = α · B · u(r) |

**To modify the viscosity model:** edit `apparent_viscosity_cp()` at line 623.
Adding a new term (e.g. a different thickening model) means adding to the return
expression — both the shear-thinning and shear-thickening branches are additive.

### polymer_dp() — line 641

Numerically integrates the pressure drop across the polymer bank:

**Formula — Abdullah 2023, Eq. 12:**
```
ΔP_poly = ∫[r_w → rp]  q · μ_app(u(r)) / (2π · h · kp · r)  dr
```
where `u(r) = q / (2π · h · r)` and `kp = k · krw / Rk`.

Integration uses 500 log-spaced radial points (controlled by `uvim_npts`).
Uses `numpy.trapezoid` (trapezoidal rule in log-r space).

**To change integration accuracy:** modify `uvim_npts` in config (default 500).
More points = slower but more accurate; 200 is usually sufficient.  
**To change integration method:** edit the `trapezoid_logspace_integral()` helper
at line 372 — you could substitute Simpson's rule or Gaussian quadrature.

---

## 7. Model 3 — UVIM + Fracture (SPE-215083)

**Class:** `VerticalUVIMFractureModel` — line 802  
**Inherits:** `UVIMVerticalModel` + `FractureCommon`  
**References:**
- SPE-215083 (Abdullah et al. 2023) — fracture propagation
- SPE-224998 (AlAbdullah et al. 2026) — Wara field validation
- Prats (1961) — equivalent wellbore radius

### Concept

When BHP exceeds the fracture initiation pressure, the model switches to a
fracture-coupled mode. A hydraulic fracture is simulated as an equivalent
wellbore radius `rwe > rw`, which reduces the near-wellbore pressure drop and
allows higher injection rates at the same BHP.

At each timestep the solver finds the fracture half-length `Lf` such that:

```
BHP(rp, rwe(Lf)) = target pressure
```

where the target is either the fracture-mechanics criterion (toughness-based)
or a user-specified BHP target.

### Fracture growth modes — line 822

**`criterion` mode** (default): fracture grows when BHP equals the fracture
closure pressure plus a toughness resistance term:
```
target = σ_hmin + P_net + KIC / √(π · Af)
```

**`bhp_target` mode**: fracture grows to hold BHP at a user-specified pressure
(`fracture_bhp_target_psi`). Useful for modelling rate-controlled injection
where a surface BHP limit applies.

**To change the growth criterion:** edit `_criterion_target_pa()` at line 818.

### solve_fracture_step() — line 873

Bisection solver (80 iterations, tolerance = `fracture_tol_psi`).
At each timestep:
1. Compute BHP without fracture — if below threshold, no fracture this step.
2. Otherwise bisect on `Lf` until `BHP(rp, rwe(Lf)) = target`.
3. If `fracture_persistence=True`, `Lf` can only grow (monotonic).

**To change solver tolerance:** modify `fracture_tol_psi` in config (default 5 psi).  
**To change solver algorithm:** replace the bisection loop in `solve_fracture_step()`
at line 873 — e.g. Brent's method would converge faster.

---

## 8. Model 4 — Horizontal well (Aitkulov 2021)

**Class:** `HorizontalAitkulovModel` — line 964  
**Inherits:** `PowerLawCommon` → `BaseModel`  
**Reference:** Aitkulov et al. (2021), *J. Pet. Sci. Eng.* 108748

### Concept

Near the horizontal wellbore, flow is radial (cylindrical). Beyond a transition
boundary at `r_boundary = h/π`, flow becomes linear (piston-like). The polymer
front position depends on which regime dominates.

```
Wellbore ──[Radial zone]──[Transition]──[Linear zone]──[Outer zone]── boundary
            rw → h/π       r_boundary    r_boundary → x_front          → well spacing
            Cylindrical                  1D linear
```

### polymer_front() — line 969

Returns the front position as a dict with `mode` ("radial" or "linear") and
coordinates `rp_m` and `x_m`.

**Radial phase** (while front is inside `r_boundary`):
```
rp = √( Wi·C / (π · L · A) + rw² )
```

**Linear phase** (once front passes `r_boundary`):
```
x = r_boundary + (Wi·C - mass_radial) / (2 · h · L · A)
```

### polymer_dp() — line 1008

Computes pressure drop in the radial zone using the same three-zone power-law
logic as Model 1, then adds the linear zone contribution.

**Linear zone ΔP** — determined by the dominant flux regime at the transition:
```
ΔP_linear = (q / (2hL · λ)) · Δx    (Newtonian regime)
ΔP_linear = coeff · (q/(2hL))^n · Δx  (shear-thinning regime)
```

**To modify the linear flow formula:** edit `_dp_linear_segment()` at line 987.

---

## 9. Shared fracture geometry (FractureCommon)

**Class:** `FractureCommon` — line 709  
**Inherits:** `PowerLawCommon`

Used only by Model 3. Contains four independent sub-models:

### pkn_solution() — line 720

**PKN fracture** (Perkins-Kern-Nordgren). Use when `Lf >> hf`.

**Net pressure — SPE-215083 Eq. 19:**
```
P_net = [ (E'/(2hf))^(2n+1) · ((4q(2n+1))/(π·hf·n))^n · (4·Ct·Lf·K) ]^(1/(2n+2))
```

**Max width at wellbore — SPE-215083 Eq. 20:**
```
wf_max = 2hf / E' · P_net
```

### kgd_solution() — line 735

**KGD fracture** (Kristianovich-Geertsma-de Klerk). Use when `Lf ~ hf`.

**Net pressure — SPE-215083 Eq. 30:**
```
P_net = [ (q/hf · (2n+1)/n)^n · 4K·Ct·Lf · (E'/(4Lf))^(2n+1) ]^(1/(2n+2))
```

**Max width — SPE-215083 Eq. 31:**
```
wf_max = 4Lf / E' · P_net
```

### fracture_geometry_parameter() — line 750

**SPE-215083, geometry regime switch:**
```
if 2Lf > hf:  Af = hf/4    (PKN regime)
else:          Af = Lf      (KGD regime)
```
Af appears in the toughness term: `KIC / √(π·Af)`.

> **Audit fix:** earlier versions always used `Af = Lf`, which overestimates
> toughness resistance in the PKN regime and delays predicted fracture initiation.

### fracture_width_average() — line 757

**Rahman & Rahman (2010):**
```
PKN: w_avg = (π/5) · wf_max
KGD: w_avg = (π/4) · wf_max
```

**To change the width correction:** edit this method. Note that `w_avg` feeds
directly into the Prats equivalent wellbore radius — changing it affects BHP.

### equivalent_wellbore_radius() — line 763

**Prats (1961) bounded correlation — SPE-215083 Eqs. 6–9:**

Step 1 — fracture permeability (parallel plates):
```
kf = w_avg² / 12
```

Step 2 — dimensionless relative capacity:
```
a = π · k · Lf / (2 · kf · w_avg)
```

Step 3 — dimensionless wellbore radius:
```
a < 0.16 : rwD = 0.5              (infinite conductivity)
a < 1.0  : rwD = 0.5 · exp(-π·a) (transitional)
a ≥ 1.0  : rwD = 0.5 / (π·a)     (finite conductivity)
```

Step 4:
```
rwe = Lf · rwD
```

`rwe` replaces `rw` in the UVIM pressure drop integral, reducing near-well ΔP.

> **Audit fix:** earlier versions used a single linear approximation without the
> correct asymptotic limits. The three-branch correlation above is the standard
> Prats result.

---

## 10. Rate schedule

**Class:** `RateSchedule` — line 318

Accepts a list of `(day, q_bpd)` tuples. Injection rate is piecewise-constant
and left-continuous: the rate on each row holds until the next row's date.

**From Excel:** `load_rate_schedule_from_excel()` at line 1130 converts calendar
dates to day numbers (day 0 = first date in file) and returns the same tuple list.

**To add new schedule logic** (e.g. ramped rates, sinusoidal): subclass
`RateSchedule` and override `q_bpd_at()`. The rest of the engine calls only
`rate_schedule.q_m3s_at(day)` — you don't need to change anything else.

---

## 11. Output columns reference

Every `run()` method returns a list of dicts, one per timestep.
All models output at minimum:

| Column | Unit | Description |
|---|---|---|
| `day` | days | Simulation time |
| `q_bpd` | bbl/d | Injection rate at this timestep |
| `model` | string | Model name |
| `rp_ft` | ft | Polymer front radius |
| `dp_poly_psi` | psi | ΔP across the polymer bank |
| `dp_outer_psi` | psi | ΔP across the outer (water/oil) zone |
| `bhp_psi` | psi | Bottom-hole pressure = Pres + dp_poly + dp_outer |
| `injectivity_bpd_per_psi` | bpd/psi | II = q / (BHP − Pres) |
| `frac` | "YES"/"NO" | Whether fracture is active |
| `Lf_ft` | ft | Fracture half-length (0 if no fracture) |
| `wfmax_in` | in | Max fracture width at wellbore |
| `wfavg_in` | in | Average fracture width |
| `rwe_ft` | ft | Equivalent wellbore radius |
| `NDe` | — | Deborah number (viscoelastic measure) |
| `Sorp_ratio` | — | Screen-out resistance reduction factor |

Model 1 additionally outputs: `dp1_psi`, `dp2_psi`, `dp3_psi`, `r_uc_ft`, `r_lc_ft`  
Model 3 additionally outputs: `fracture_growth_mode`, `fracture_target_psi`, `fracture_residual_psi`  
Model 4 additionally outputs: `front_mode`, `x_ft`, `radial_linear_boundary_ft`, `dp_linear_psi`

---

## 12. How to make common modifications

### Change the polymer rheology model

The power-law model (Models 1 and 4) is in `dp_shear_thinning_radial()` at
line 514. To substitute a different viscosity function (e.g. Ellis model):
1. Modify `PowerLawCommon.__init__()` to pre-compute your model's constants.
2. Replace the formula in `dp_shear_thinning_radial()`.
3. Update `radial_critical_radii()` to use your model's transition velocities.

For the UVIM (Models 2 and 3), edit `apparent_viscosity_cp()` at line 623.
The integration in `polymer_dp()` calls it automatically — no other changes needed.

### Add a new fracture model (e.g. PKN-C with fluid loss)

1. Add a new method in `FractureCommon` following the pattern of `pkn_solution()`.
2. Add the new name to the `fracture_model` validation in `normalize_config()` (line 260).
3. Update `fracture_state_from_length()` at line 789 to call your method.
4. Add the new option to `frac_model` in `streamlit_app.py` and `gui.py`.

### Add a new output column

In the relevant `run()` method, add the key/value to the `rows.append(dict(...))` call.
That's all — it will automatically appear in CSV exports and the Streamlit data table.

### Change the fracture initiation criterion

Edit `_criterion_target_pa()` at line 818. Currently:
```python
target = pf + KIC / sqrt(π · Af)
```
Replace the right-hand side with your preferred criterion (e.g. net-pressure-based,
or a fixed override pressure).

### Add temperature dependence to viscosity

The current engine is isothermal. To add temperature:
1. Add `T_reservoir_C` to `BASE_DEFAULTS`.
2. In `PowerLawCommon.__init__()`, adjust `K` and `n` using a T-dependent
   correlation (e.g. Yasuda model).
3. Similarly adjust UVIM parameters in `UVIMVerticalModel.__init__()`.

### Change how the polymer bank radius is computed

Edit `polymer_bank_radius_vertical()` at line 430 in `BaseModel`.
For a non-piston-like displacement (e.g. fractional-flow-based front tracking),
this is the method to replace.

---

## 13. Known limitations and assumptions

| Assumption | Where it applies | Alternative if needed |
|---|---|---|
| Quasi-steady-state radial flow | All models | Transient solution (requires PDE solver) |
| Single layer, uniform properties | All models | Multilayer: sum ΔP across layers with separate rp per layer |
| Incompressible fluids | All models | Add compressibility term to storage_coeff |
| Piston-like displacement | All models (rp calculation) | Fractional flow + Buckley-Leverett |
| Constant injection rate per timestep | All models | Already handled by RateSchedule |
| No gravity / dip | All models | Add ρgh term to bhp calculation |
| No polymer degradation in situ | All models | Add time-decay to μ₀ or K in run() loop |
| Planar 2D fracture (PKN/KGD) | Model 3 | 3D fracture requires numerical simulator |
| No fluid loss (leakoff) in fracture | Model 3 | Add Carter leakoff to pkn_solution() / kgd_solution() |
| Constant Rk throughout bank | All models | Make Rk a function of local concentration |

---

*Last updated: June 2026 — Antoine, eppok EURL*  
*Engine version: 1.0.0 — all 6 audit fixes applied, 58 tests passing*
