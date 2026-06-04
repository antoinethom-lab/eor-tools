"""
polymer_injectivity.presets
============================
Single source of truth for the field-case presets and model labels,
shared by the desktop GUI (gui.py) and the web app (streamlit_app.py).

Each preset is a complete model config dict ready to pass to make_model().
"""

from .engine import FT2M


# ── Parameter layout: (label, config_key, unit, tooltip) ────────────────────
PARAM_GROUPS = [
    ("Reservoir", [
        ("Permeability k",           "k_md",               "md",    "Absolute permeability"),
        ("Pay thickness h",          "h_m",                "m",     "Net pay zone height"),
        ("Porosity φ",               "phi",                "—",     "Fractional porosity (0–1)"),
        ("Water sat. Sw",            "Sw",                 "—",     "Connate water saturation"),
        ("Residual oil Sorw",        "Sorw",               "—",     "Residual oil sat behind front"),
        ("Endpoint krw",             "krwmax",             "—",     "Max krw at residual oil"),
        ("Endpoint kro",             "kromax",             "—",     "Max kro at connate water"),
        ("Drainage radius re",       "re_ft",              "ft",    "Radial boundary distance"),
        ("Wellbore radius rw",       "rw_m",               "m",     "Physical wellbore radius"),
        ("Boundary pressure",        "boundary_pressure_psi", "psi", "Static reservoir pressure"),
    ]),
    ("Polymer (Power-Law)", [
        ("Water viscosity μw",       "mu_w_cp",            "cp",    "Aqueous phase viscosity"),
        ("Low-shear viscosity μ₀",   "mu0_cp",             "cp",    "Zero-shear polymer viscosity"),
        ("Consistency K",            "powerlaw_K_cp",      "cp·sⁿ⁻¹", "Power-law prefactor"),
        ("Flow index n",             "powerlaw_n",         "—",     "Power-law exponent (n<1 = shear-thin)"),
        ("Concentration Cp",         "Cp_ppm",             "ppm",   "Polymer concentration"),
        ("Adsorption D",             "D_ug_g",             "μg/g",  "Rock adsorption"),
        ("IPV",                      "IPV",                "—",     "Inaccessible pore volume fraction"),
        ("Cannella C",               "C_corr",             "—",     "Shear rate correction factor"),
        ("Resistance Rk",            "Rk",                 "—",     "Permeability reduction factor"),
    ]),
    ("UVIM (Viscoelastic)", [
        ("μ∞ (high shear)",          "uvim_muinf_cp",      "cp",    "Infinite-shear viscosity (≈ μw)"),
        ("μ₀ᵖ (zero shear)",        "uvim_mu0_cp",        "cp",    "Zero-shear polymer viscosity"),
        ("μmax (elongational)",      "uvim_mu_max_cp",     "cp",    "Max shear-thickening viscosity"),
        ("λ (shear-thin time)",      "uvim_lambda_s",      "s",     "Shear-thinning time constant"),
        ("λ₂ (thickening param)",    "uvim_lambda2",       "—",     "Shear-thickening parameter (≈0.01)"),
        ("n₁ (thin index)",          "uvim_n1",            "—",     "Shear-thinning index (Carreau n)"),
        ("n₂ (thick exponent)",      "uvim_n2",            "—",     "Shear-thickening exponent"),
        ("τ (relaxation time)",      "tau_s",              "s",     "In-situ viscoelastic relaxation time"),
    ]),
    ("Injection & Time", [
        ("Rate q",                   "q_bpd",              "bbl/d", "Constant injection rate"),
        ("Duration",                 "duration_d",         "days",  "Simulation duration"),
        ("Timestep",                 "dt_d",               "days",  "Output timestep"),
    ]),
    ("Fracture (SPE-215083)", [
        ("σhmin",                    "sigma_hmin_psi",     "psi",   "Min horizontal stress (fracture init)"),
        ("Frac init pressure",       "fracture_initiation_pressure_psi", "psi", "Pressure used as fracture-initiation baseline"),
        ("Young's modulus E",        "E_psi",              "psi",   "Rock elastic modulus"),
        ("Poisson's ratio ν",        "nu_pr",              "—",     "Poisson ratio (0.2–0.35)"),
        ("Toughness KIc",            "KIC_psi_sqrt_in",    "psi√in", "Fracture toughness"),
        ("Turbulence Cturb",         "Cturb",              "—",     "Turbulence correction"),
        ("Fracture height hf",       "hf_m",               "m",     "Fracture height (default = h)"),
        ("Max fracture Lf",          "fracture_max_length_ft", "ft", "Maximum fracture half-length allowed by solver"),
        ("Fracture tolerance",       "fracture_tol_psi",   "psi",   "Fracture-length solver tolerance"),
        ("BHP target",               "fracture_bhp_target_psi", "psi", "Used only in bhp_target growth mode"),
    ]),
    ("Horizontal Well", [
        ("Horiz. length L",          "L_horiz_m",          "m",     "Horizontal wellbore length"),
        ("Well spacing",             "well_dist_m",        "m",     "Distance to offset well"),
    ]),
    ("Densities & Other", [
        ("Water density ρw",         "rho_w",              "kg/m³", "Brine density"),
        ("Rock density ρr",          "rho_r",              "kg/m³", "Rock grain density"),
        ("Skin factor",              "skin",               "—",     "Near-wellbore skin"),
    ]),
]

INT_PARAMS = {"duration_d", "uvim_npts"}

# ── Field presets from the papers ───────────────────────────────────────────
FIELD_PRESETS = {
    "Belayim (SPE-195513)": dict(
        model="vertical_delamaide",
        k_md=1150, h_m=5.0, phi=0.19, Sw=0.15, Sorw=0.30,
        krwmax=0.86, mu_w_cp=1.0, rho_r=2650, rho_w=1020,
        boundary_pressure_psi=1500, re_ft=330, rw_m=0.10,
        Cp_ppm=1500, powerlaw_K_cp=5.5, powerlaw_n=0.835, mu0_cp=5.0,
        D_ug_g=62.9, IPV=0.0, Rk=4.4, C_corr=1.3,
        q_bpd=1000, duration_d=260, dt_d=2,
        fracture_on=False,
    ),
    "Matzen UVIM+Frac (SPE-215083)": dict(
        model="vertical_uvim_fracture",
        k_md=550, h_m=16.5 * FT2M, hf_m=16.5 * FT2M,
        phi=0.26, Sw=0.87, Sorw=0.13,
        krwmax=0.30, mu_w_cp=0.6, rho_r=2650, rho_w=1020,
        boundary_pressure_psi=1552, re_ft=984, rw_m=0.30 * FT2M,
        Cp_ppm=1000, powerlaw_K_cp=7.0, powerlaw_n=0.85, mu0_cp=7.8,
        D_ug_g=0.0, IPV=0.0, Rk=1.0, C_corr=1.0,
        q_bpd=2400, duration_d=53, dt_d=1,
        uvim_muinf_cp=0.6, uvim_mu0_cp=7.8, uvim_n1=0.4,
        uvim_lambda_s=0.65, uvim_mu_max_cp=33.0,
        uvim_lambda2=0.01, uvim_n2=2.0, tau_s=0.6,
        sigma_hmin_psi=2650, E_psi=1060000, nu_pr=0.365,
        KIC_psi_sqrt_in=0.0, Cturb=1.7, fracture_model="PKN",
    ),
    "Wara UVIM+Frac (SPE-224998)": dict(
        model="vertical_uvim_fracture",
        k_md=650, h_m=43 * FT2M, hf_m=43 * FT2M,
        phi=0.18, Sw=0.70, Sorw=0.30,
        krwmax=0.40, mu_w_cp=1.1, rho_r=2650, rho_w=1020,
        boundary_pressure_psi=1581, re_ft=400, rw_m=0.30 * FT2M,
        Cp_ppm=1800, powerlaw_K_cp=16, powerlaw_n=0.80, mu0_cp=10.0,
        D_ug_g=0.0, IPV=0.0, Rk=1.17, C_corr=1.0,
        q_bpd=3200, duration_d=7, dt_d=0.25,
        uvim_muinf_cp=1.1, uvim_mu0_cp=10.0, uvim_n1=0.80,
        uvim_lambda_s=0.1, uvim_mu_max_cp=65.0,
        uvim_lambda2=0.01, uvim_n2=3.3, tau_s=0.35,
        sigma_hmin_psi=2300, E_psi=7000000, nu_pr=0.21,
        KIC_psi_sqrt_in=2500, Cturb=1.7, fracture_model="PKN",
    ),
    "Pelican Lake Horiz (SPE-195513)": dict(
        model="horizontal_aitkulov",
        k_md=610, h_m=4.5, phi=0.29, Sw=0.30, Sorw=0.35,
        krwmax=0.05, mu_w_cp=1.0, rho_r=2650, rho_w=1020,
        boundary_pressure_psi=580, re_ft=600, rw_m=0.10,
        Cp_ppm=2000, powerlaw_K_cp=72, powerlaw_n=0.6, mu0_cp=20.0,
        D_ug_g=25, IPV=0.0, Rk=2.0, C_corr=1.0,
        q_bpd=800, duration_d=400, dt_d=2,
        L_horiz_m=244, well_dist_m=300,
        outer_zone_mu_cp=10.0, outer_zone_kr=0.30,
        fracture_on=False,
    ),
}

# Backward-compatible alias
PRESETS = FIELD_PRESETS

MODEL_LABELS = {
    "vertical_delamaide":     "Delamaide/Lake (shear-thin only)",
    "vertical_uvim":          "UVIM (full viscoelastic, no frac)",
    "vertical_uvim_fracture": "UVIM + Fracture (PKN/KGD)",
    "horizontal_aitkulov":    "Horizontal Well (Aitkulov)",
}
