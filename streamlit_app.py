#!/usr/bin/env python3
"""
streamlit_app.py
================
Web app version of the Polymer Injectivity Calculator — eppok EURL.

Run locally:
    streamlit run streamlit_app.py

Deploy free on Streamlit Community Cloud (see DEPLOY_WEB.md).

Password protection:
    Set a password in .streamlit/secrets.toml (local) or in the
    Streamlit Cloud "Secrets" panel:

        app_password = "your-password-here"

    If no password is configured, the app runs open (useful for local dev).
"""

import io
import hmac
import datetime

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from polymer_injectivity.engine import (
    make_model,
    load_rate_schedule_from_excel,
    rate_schedule_summary,
    BASE_DEFAULTS,
)
from polymer_injectivity.presets import (
    PARAM_GROUPS,
    INT_PARAMS,
    PRESETS,
    MODEL_LABELS,
)

# ──────────────────────────────────────────────────────────────────────────────
#  NULLABLE PARAMS — these must stay None so the engine computes cross-defaults
#  (e.g. uvim_muinf_cp inherits from mu_w_cp when None).
#  We show them as optional text inputs, not number inputs.
# ──────────────────────────────────────────────────────────────────────────────
NULLABLE_PARAMS = {
    "uvim_muinf_cp",
    "uvim_mu0_cp",
    "uvim_n1",
    "hf_m",
    "outer_zone_length_m",
    "fracture_initiation_pressure_psi",
    "fracture_bhp_target_psi",
}

# ──────────────────────────────────────────────────────────────────────────────
#  PAGE CONFIG + THEME
# ──────────────────────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="Polymer Injectivity Calculator — eppok",
    page_icon="🛢️",
    layout="wide",
    initial_sidebar_state="expanded",
)

NAVY  = "#1F3A5F"
NAVY_DARK = "#13243B"
ORANGE = "#E87722"
LIGHT  = "#F5F7FA"

st.markdown(f"""
<style>
    .stApp {{ background: {LIGHT}; }}
    .eppok-header {{
        background: linear-gradient(100deg, {NAVY_DARK} 0%, {NAVY} 70%);
        padding: 1.4rem 1.8rem; border-radius: 12px; margin-bottom: 1.2rem;
        border-left: 6px solid {ORANGE};
    }}
    .eppok-header h1 {{ color:#FFFFFF; font-size:1.55rem; margin:0; font-weight:700; }}
    .eppok-header p  {{ color:#B9C6D6; font-size:0.92rem; margin:0.25rem 0 0 0; }}
    .eppok-header .accent {{ color:{ORANGE}; }}
    div[data-testid="stMetric"] {{
        background:#FFFFFF; border:1px solid #E2E8F0; border-radius:10px;
        padding:0.9rem 1.1rem; box-shadow:0 1px 3px rgba(19,36,59,0.06);
    }}
    .stButton > button[kind="primary"] {{
        background:{ORANGE}; border:none; font-weight:700; font-size:1.02rem;
        padding:0.55rem 0; border-radius:8px;
    }}
    .stButton > button[kind="primary"]:hover {{ background:#C9631A; }}
    section[data-testid="stSidebar"] {{ background:#FFFFFF; border-right:1px solid #E2E8F0; }}
    details summary {{ font-weight:600; color:{NAVY}; }}
    div[data-testid="stNumberInput"] label p {{ font-size:0.82rem; }}
    .footer-note {{
        color:#8595A8; font-size:0.8rem; text-align:center;
        margin-top:2rem; padding-top:1rem; border-top:1px solid #E2E8F0;
    }}
    .nullable-hint {{ font-size:0.75rem; color:#8595A8; margin-top:-0.5rem; }}
</style>
""", unsafe_allow_html=True)


# ──────────────────────────────────────────────────────────────────────────────
#  PASSWORD GATE
# ──────────────────────────────────────────────────────────────────────────────

def check_password() -> bool:
    configured = None
    try:
        configured = st.secrets.get("app_password", None)
    except Exception:
        configured = None

    if not configured:
        return True
    if st.session_state.get("password_ok", False):
        return True

    _, mid, _ = st.columns([1, 1.3, 1])
    with mid:
        st.markdown(f"""
        <div class="eppok-header" style="text-align:center; margin-top:3rem;">
            <h1>🛢️ Polymer Injectivity Calculator</h1>
            <p><span class="accent">eppok EURL</span> · secure access</p>
        </div>
        """, unsafe_allow_html=True)
        pw = st.text_input("Password", type="password", placeholder="Enter access password")
        if st.button("Enter", type="primary", use_container_width=True):
            if hmac.compare_digest(pw, str(configured)):
                st.session_state["password_ok"] = True
                st.rerun()
            else:
                st.error("Incorrect password.")
        st.caption("Contact eppok EURL for access credentials.")
    return False


if not check_password():
    st.stop()


# ──────────────────────────────────────────────────────────────────────────────
#  HEADER
# ──────────────────────────────────────────────────────────────────────────────

st.markdown(f"""
<div class="eppok-header">
    <h1>🛢️ Polymer Injectivity Calculator</h1>
    <p><span class="accent">eppok EURL</span> · Analytical toolkit for polymer flooding injectivity
    · Delamaide · UVIM · UVIM+Fracture · Aitkulov</p>
</div>
""", unsafe_allow_html=True)


# ──────────────────────────────────────────────────────────────────────────────
#  SESSION STATE
# ──────────────────────────────────────────────────────────────────────────────

def _initial_value(key):
    """Return the right initial value for a parameter — None for nullable fields."""
    if key in NULLABLE_PARAMS:
        return None
    v = BASE_DEFAULTS.get(key)
    if v is not None and isinstance(v, (int, float)):
        return float(v)
    return 0.0


if "params" not in st.session_state:
    st.session_state.params = {
        key: _initial_value(key)
        for _, plist in PARAM_GROUPS
        for _, key, _, _ in plist
    }

if "excel_schedule" not in st.session_state:
    st.session_state.excel_schedule = None


def apply_preset(preset_name: str):
    preset = PRESETS[preset_name]
    # Reset nullable params to None first so engine uses cross-defaults
    for key in NULLABLE_PARAMS:
        st.session_state.params[key] = None
    # Apply preset values
    for _, plist in PARAM_GROUPS:
        for _, key, _, _ in plist:
            if key in preset:
                val = preset[key]
                st.session_state.params[key] = float(val) if val is not None else None
    st.session_state.preset_model     = preset.get("model", "vertical_uvim_fracture")
    st.session_state.preset_frac_model = preset.get("fracture_model", "PKN")
    st.session_state.preset_frac_on   = preset.get("fracture_on", True)
    st.session_state.excel_schedule   = None


# ──────────────────────────────────────────────────────────────────────────────
#  SIDEBAR
# ──────────────────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown("## ⚙️ Configuration")

    preset_choice = st.selectbox(
        "📋 Load a field preset",
        ["— none —"] + list(PRESETS.keys()),
        help="Pre-filled parameters from published field cases.",
    )
    if st.button("Apply preset", use_container_width=True) and preset_choice != "— none —":
        apply_preset(preset_choice)
        st.success(f"Loaded: {preset_choice}")

    st.divider()

    model_keys = list(MODEL_LABELS.keys())
    default_model_idx = (
        model_keys.index(st.session_state["preset_model"])
        if st.session_state.get("preset_model") in model_keys else 2
    )
    model = st.selectbox(
        "🧮 Model",
        model_keys,
        index=default_model_idx,
        format_func=lambda k: MODEL_LABELS[k],
    )

    st.divider()
    st.markdown("### Options")
    frac_on      = st.checkbox("Enable fracture",
                               value=st.session_state.get("preset_frac_on", True))
    frac_model   = st.radio("Fracture model", ["PKN", "KGD"],
                            index=0 if st.session_state.get("preset_frac_model","PKN")=="PKN" else 1,
                            horizontal=True)
    frac_growth  = st.radio("Growth mode", ["criterion", "bhp_target"], horizontal=True)
    frac_persist = st.checkbox("Fracture persistence (monotonic growth)", value=True)
    tertiary     = st.checkbox("Tertiary flood (post-waterflood Sw)", value=True,
                               help="Uncheck for secondary flood.")
    use_krw_wf   = st.checkbox("Use krw in waterflood zone", value=True)

    st.divider()
    st.caption("eppok EURL · polymer injectivity v1.0")


# ──────────────────────────────────────────────────────────────────────────────
#  PARAMETER INPUTS
# ──────────────────────────────────────────────────────────────────────────────

RELEVANT_GROUPS = {
    "vertical_delamaide":     {"Reservoir","Polymer (Power-Law)","Injection & Time","Densities & Other"},
    "vertical_uvim":          {"Reservoir","Polymer (Power-Law)","UVIM (Viscoelastic)","Injection & Time","Densities & Other"},
    "vertical_uvim_fracture": {"Reservoir","Polymer (Power-Law)","UVIM (Viscoelastic)","Injection & Time","Fracture (SPE-215083)","Densities & Other"},
    "horizontal_aitkulov":    {"Reservoir","Polymer (Power-Law)","Injection & Time","Horizontal Well","Densities & Other"},
}
active_groups = RELEVANT_GROUPS.get(model, {g for g, _ in PARAM_GROUPS})

st.markdown("### 📊 Parameters")
group_names = [g for g, _ in PARAM_GROUPS if g in active_groups]
tabs = st.tabs(group_names)

for tab, gname in zip(tabs, group_names):
    with tab:
        plist = next(pl for g, pl in PARAM_GROUPS if g == gname)
        cols  = st.columns(3)
        for i, (label, key, unit, tip) in enumerate(plist):
            with cols[i % 3]:
                unit_label = f"{label}" + (f"  ({unit})" if unit != "—" else "")
                current    = st.session_state.params.get(key)

                if key in NULLABLE_PARAMS:
                    # Show as text input: empty string = None (engine uses cross-default)
                    placeholder = "auto"
                    raw = st.text_input(
                        unit_label,
                        value="" if current is None else str(current),
                        key=f"input_{key}",
                        help=tip + " · Leave blank to use automatic default.",
                        placeholder=placeholder,
                    )
                    st.markdown('<p class="nullable-hint">blank = auto</p>',
                                unsafe_allow_html=True)
                    if raw.strip() == "":
                        st.session_state.params[key] = None
                    else:
                        try:
                            st.session_state.params[key] = float(raw)
                        except ValueError:
                            st.session_state.params[key] = None
                else:
                    is_int = key in INT_PARAMS
                    val    = float(current) if current is not None else 0.0
                    new_val = st.number_input(
                        unit_label,
                        value=val,
                        key=f"input_{key}",
                        help=tip,
                        format="%.4g",
                    )
                    st.session_state.params[key] = int(new_val) if is_int else new_val


# ──────────────────────────────────────────────────────────────────────────────
#  RATE SCHEDULE
# ──────────────────────────────────────────────────────────────────────────────

st.markdown("### 📅 Injection rate schedule (optional)")
rs_col1, rs_col2 = st.columns([1, 1])

with rs_col1:
    st.markdown("**Option A — Upload Excel**")
    xl_c1, xl_c2 = st.columns(2)
    with xl_c1:
        date_col = st.text_input("Date column name", value="date")
    with xl_c2:
        rate_col = st.text_input("Rate column name", value="q_bpd")

    uploaded = st.file_uploader(
        "Upload .xlsx with date + rate columns",
        type=["xlsx", "xls"],
        help="Two columns: a date and an injection rate (bbl/day).",
    )
    if uploaded is not None:
        try:
            schedule = load_rate_schedule_from_excel(
                uploaded,
                date_col=date_col.strip() or "date",
                rate_col=rate_col.strip() or "q_bpd",
            )
            st.session_state.excel_schedule = schedule
            st.success(f"✅ {rate_schedule_summary(schedule)}")
        except Exception as e:
            st.session_state.excel_schedule = None
            st.error(f"Could not load file: {e}")

with rs_col2:
    st.markdown("**Option B — Type manually**")
    manual_text = st.text_area(
        "One `day, rate` pair per line",
        placeholder="0, 800\n30, 1000\n90, 1200",
        height=150,
    )


def build_template_bytes() -> bytes:
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    wb = Workbook(); ws = wb.active; ws.title = "Rate Schedule"
    fill = PatternFill("solid", fgColor="1F3A5F")
    for c, h in enumerate(["date","q_bpd","notes"], start=1):
        cell = ws.cell(row=1, column=c, value=h)
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = fill
        cell.alignment = Alignment(horizontal="center")
    ws.column_dimensions["A"].width = 16
    ws.column_dimensions["B"].width = 12
    ws.column_dimensions["C"].width = 30
    for r, (d, q, n) in enumerate([
        (datetime.date(2024,1,1), 800,  "Injection start"),
        (datetime.date(2024,2,1), 1000, "Step increase"),
        (datetime.date(2024,4,1), 1200, "Full rate"),
        (datetime.date(2024,7,1), 1500, "Plateau"),
    ], start=2):
        ws.cell(row=r, column=1, value=d).number_format = "YYYY-MM-DD"
        ws.cell(row=r, column=2, value=q)
        ws.cell(row=r, column=3, value=n)
    buf = io.BytesIO(); wb.save(buf); return buf.getvalue()

st.download_button(
    "⬇️ Download Excel template",
    data=build_template_bytes(),
    file_name="rate_schedule_template.xlsx",
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
)


# ──────────────────────────────────────────────────────────────────────────────
#  RUN
# ──────────────────────────────────────────────────────────────────────────────

st.markdown("")
run = st.button("▶  RUN SIMULATION", type="primary", use_container_width=True)


def assemble_config() -> dict:
    cfg = {}
    for _, plist in PARAM_GROUPS:
        for _, key, _, _ in plist:
            val = st.session_state.params.get(key)
            if val is not None:          # only pass values the user actually set
                cfg[key] = val
            # None → omit → engine uses its own cross-default

    cfg["model"]              = model
    cfg["fracture_on"]        = frac_on
    cfg["fracture_model"]     = frac_model
    cfg["fracture_growth_mode"] = frac_growth
    cfg["fracture_persistence"] = frac_persist
    cfg["tertiary_flood"]     = tertiary
    cfg["use_krw_waterflood"] = use_krw_wf

    if st.session_state.excel_schedule:
        cfg["rate_schedule"] = st.session_state.excel_schedule
    elif manual_text.strip():
        sched = []
        for line in manual_text.strip().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(";",",").split(",")
            if len(parts) >= 2:
                try:
                    sched.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    pass
        if sched:
            cfg["rate_schedule"] = sched
    return cfg


def make_figure(rows):
    days = [r["day"]                          for r in rows]
    bhp  = [r["bhp_psi"]                      for r in rows]
    rp   = [r.get("rp_ft", 0.0)               for r in rows]
    ii   = [r.get("injectivity_bpd_per_psi", 0.0) for r in rows]
    Lf   = [r.get("Lf_ft", 0.0)               for r in rows]

    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), facecolor="white")
    ax = axes.ravel()
    for a, (series, title, ylab, color) in zip(ax, [
        (bhp, "Bottom-hole pressure",  "psi",     NAVY),
        (rp,  "Polymer front extent",  "ft",      ORANGE),
        (ii,  "Injectivity index",     "bpd/psi", "#2E7D5B"),
        (Lf,  "Fracture half-length",  "ft",      "#9C2B2B"),
    ]):
        a.plot(days, series, color=color, linewidth=2)
        a.set_title(title, fontsize=11, fontweight="bold", color=NAVY)
        a.set_xlabel("day", fontsize=9)
        a.set_ylabel(ylab, fontsize=9)
        a.grid(True, alpha=0.25)
        a.set_facecolor("#FCFCFD")
    fig.tight_layout()
    return fig


if run:
    try:
        cfg = assemble_config()
        with st.spinner("Running simulation..."):
            model_obj = make_model(cfg)
            rows = model_obj.run()
        st.session_state.results     = rows
        st.session_state.results_cfg = cfg
    except Exception as e:
        st.error(f"Simulation error: {e}")
        st.session_state.results = None


# ──────────────────────────────────────────────────────────────────────────────
#  RESULTS
# ──────────────────────────────────────────────────────────────────────────────

if st.session_state.get("results"):
    rows = st.session_state.results
    first, last = rows[0], rows[-1]

    st.markdown("### 📈 Results")

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Final BHP",
              f"{last['bhp_psi']:,.0f} psi",
              f"{last['bhp_psi']-first['bhp_psi']:+,.0f} psi")
    m2.metric("Final injectivity",
              f"{last['injectivity_bpd_per_psi']:.3f} bpd/psi")
    m3.metric("Polymer front",
              f"{last.get('rp_ft',0):,.1f} ft")
    max_lf = max(r.get("Lf_ft",0.0) for r in rows)
    m4.metric("Max fracture Lf",
              f"{max_lf:,.1f} ft" if max_lf > 0 else "— no frac")

    fig = make_figure(rows)
    st.pyplot(fig, use_container_width=True)

    df = pd.DataFrame(rows)
    with st.expander("📋 Full data table"):
        st.dataframe(df, use_container_width=True, height=320)

    dl1, dl2 = st.columns(2)
    with dl1:
        st.download_button("⬇️ Download CSV",
                           data=df.to_csv(index=False).encode("utf-8"),
                           file_name="polymer_injectivity_results.csv",
                           mime="text/csv", use_container_width=True)
    with dl2:
        png_buf = io.BytesIO()
        fig.savefig(png_buf, format="png", dpi=150, bbox_inches="tight")
        st.download_button("⬇️ Download plot (PNG)",
                           data=png_buf.getvalue(),
                           file_name="polymer_injectivity_plot.png",
                           mime="image/png", use_container_width=True)
    plt.close(fig)
else:
    st.info("Set your parameters above (or load a preset), then click **Run simulation**.")


st.markdown(
    '<div class="footer-note">eppok EURL · Polymer injectivity analytical toolkit · '
    'Models: SPE-195513, SPE-215083, Abdullah et al. 2023, Aitkulov et al. 2021</div>',
    unsafe_allow_html=True,
)
