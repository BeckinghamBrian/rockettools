"""
╔══════════════════════════════════════════════════════════════════════════════╗
║        HALL EFFECT THRUSTER  ·  ANALYSIS SYSTEM  ·  v3.0                  ║
║        Electric Propulsion Performance Simulator                             ║
╚══════════════════════════════════════════════════════════════════════════════╝

Monatomic propellants : Xenon (Xe), Krypton (Kr), Argon (Ar)
Molecular propellants : Nitrogen (N₂), Oxygen (O₂), Water (H₂O)

Physics model
─────────────
  • Ion acceleration through crossed E×B discharge
  • Voltage-dependent ionization efficiency with dissociation penalties
  • Beam divergence correction as function of B-field strength
  • Bohm anomalous electron diffusion across field lines
  • Double-stage acceleration mode
  • Three magnetic field profile shapes (uniform / Gaussian / double-peak)

Install
───────
    pip install numpy scipy matplotlib

Run
───
    python hall_thruster_analyzer.py
"""

import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np

# ── auto-select a working GUI backend ────────────────────────────────────────
import matplotlib
_BACKENDS = ["TkAgg", "Qt5Agg", "QtAgg", "WXAgg", "GTK3Agg", "Agg"]
for _b in _BACKENDS:
    try:
        matplotlib.use(_b)
        break
    except Exception:
        continue

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.widgets import Slider, RadioButtons, CheckButtons, Button
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch
from matplotlib.collections import LineCollection
from dataclasses import dataclass
from typing import Optional


# ══════════════════════════════════════════════════════════════════════════════
# PHYSICAL CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

E_C  = 1.602176634e-19   # C  — elementary charge
ME   = 9.10938e-31       # kg — electron mass
KB   = 1.380649e-23      # J/K
G0   = 9.80665           # m/s²
NA   = 6.02214076e23     # mol⁻¹


# ══════════════════════════════════════════════════════════════════════════════
# PROPELLANT DATABASE
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class Propellant:
    name:                str
    symbol:              str
    molar_mass:          float           # kg/mol
    ionization_energy:   float           # eV  (first ionization potential)
    prop_type:           str             # 'monatomic' | 'molecular'
    color:               str             # matplotlib hex color
    dissociation_energy: Optional[float] = None   # eV  (None for monatomic)
    degrees_of_freedom:  int  = 3        # internal DoF (3 mono, 5-6 mol)
    max_ioniz_eff:       float = 0.95    # ceiling on ionization fraction

    @property
    def ion_mass_kg(self) -> float:
        """Single-ion mass [kg] with molecular fragmentation correction."""
        m_mol = self.molar_mass / NA
        if self.prop_type == "monatomic":
            return m_mol
        # Partial dissociation at thruster temperatures:
        # Water fragments into 3 pieces, diatomics into 2.
        frag  = 0.50 if self.symbol == "H₂O" else 0.60
        frags = 3    if self.symbol == "H₂O" else 2
        return m_mol / (1.0 + frag * (frags - 1.0))


PROPELLANTS = {
    # ── Monatomic noble gases ─────────────────────────────────────────────────
    "Xenon":    Propellant("Xenon",    "Xe",  0.13129, 12.13, "monatomic", "#4EA8DE", None,  3, 0.95),
    "Krypton":  Propellant("Krypton",  "Kr",  0.08380, 14.00, "monatomic", "#52B788", None,  3, 0.92),
    "Argon":    Propellant("Argon",    "Ar",  0.03995, 15.76, "monatomic", "#F4A261", None,  3, 0.88),
    # ── Molecular propellants ─────────────────────────────────────────────────
    "Nitrogen": Propellant("Nitrogen", "N₂",  0.02802, 15.58, "molecular", "#C77DFF", 9.76,  5, 0.75),
    "Oxygen":   Propellant("Oxygen",   "O₂",  0.03200, 12.07, "molecular", "#E07A5F", 5.12,  5, 0.72),
    "Water":    Propellant("Water",    "H₂O", 0.01802, 12.62, "molecular", "#70C1B3", 9.51,  6, 0.65),
}

PROP_NAMES = list(PROPELLANTS.keys())


# ══════════════════════════════════════════════════════════════════════════════
# PHYSICS ENGINE
# ══════════════════════════════════════════════════════════════════════════════

class HETPhysics:
    """
    Semi-analytical Hall Effect Thruster performance model.

    Parameters are stored as plain attributes so the dashboard can
    modify them freely and call any method to get updated results.
    """

    def __init__(self):
        self.prop_name:       str   = "Xenon"
        self.voltage:         float = 300.0   # V
        self.power:           float = 1500.0  # W
        self.B_field:         float = 150.0   # Gauss
        self.channel_length:  float = 25.0    # mm
        self.channel_width:   float = 20.0    # mm
        self.mean_radius:     float = 50.0    # mm
        self.B_profile:       int   = 0       # 0=uniform 1=peaked 2=double
        self.double_stage:    bool  = False
        self.stage2_voltage:  float = 100.0   # V

    # ── convenience accessors ─────────────────────────────────────────────────

    @property
    def prop(self) -> Propellant:
        return PROPELLANTS[self.prop_name]

    @property
    def channel_area_m2(self) -> float:
        """Annular discharge channel cross-section [m²]."""
        r_out = (self.mean_radius + self.channel_width / 2) * 1e-3
        r_in  = (self.mean_radius - self.channel_width / 2) * 1e-3
        return np.pi * (r_out**2 - r_in**2)

    # ── ionization efficiency ─────────────────────────────────────────────────

    def ioniz_eff(self, prop_name: Optional[str] = None,
                  voltage: Optional[float] = None) -> float:
        """
        Fraction of propellant atoms successfully ionized [0–1].

        Falls for:
          • propellants with high ionization energy relative to Vd
          • molecular propellants that must first dissociate
        """
        p  = PROPELLANTS[prop_name] if prop_name else self.prop
        Vd = voltage if voltage is not None else self.voltage

        volt_margin = np.clip((Vd - p.ionization_energy) / max(Vd, 1.0), 0.05, 1.0)
        eff = p.max_ioniz_eff * (1.0 - np.exp(-3.0 * volt_margin))

        if p.prop_type == "molecular" and p.dissociation_energy:
            diss_penalty = p.dissociation_energy / max(Vd, 1.0)
            eff *= max(0.3, 1.0 - 0.8 * diss_penalty)

        return float(np.clip(eff, 0.05, p.max_ioniz_eff))

    # ── effective acceleration voltage ────────────────────────────────────────

    def accel_voltage(self, prop_name: Optional[str] = None,
                      voltage: Optional[float] = None) -> float:
        """
        Net ion acceleration voltage after plasma potential losses [V].
        """
        p  = PROPELLANTS[prop_name] if prop_name else self.prop
        Vd = voltage if voltage is not None else self.voltage

        plasma_loss = 8.0                        # sheath + plasma potential drop
        ion_cost    = p.ionization_energy * 0.15 # energy diverted to ionization per ion
        Veff = Vd - plasma_loss - ion_cost
        if self.double_stage:
            Veff += self.stage2_voltage * 0.85
        return max(Veff, 10.0)

    # ── exhaust velocity ──────────────────────────────────────────────────────

    def exhaust_velocity(self, prop_name: Optional[str] = None,
                         voltage: Optional[float] = None) -> float:
        """Ion exhaust velocity [m/s] = sqrt(2 e Veff / mi)."""
        p    = PROPELLANTS[prop_name] if prop_name else self.prop
        Veff = self.accel_voltage(prop_name, voltage)
        return float(np.sqrt(2.0 * E_C * Veff / p.ion_mass_kg))

    # ── beam divergence correction ────────────────────────────────────────────

    def divergence_factor(self, B_field: Optional[float] = None) -> float:
        """
        cos²(θ) beam divergence correction.
        Higher B → tighter beam → less divergence loss.
        """
        B_T  = (B_field if B_field is not None else self.B_field) * 1e-4  # Gauss → T
        theta_deg = 12.0 + 20.0 * (1.0 - np.clip(B_T / 0.02, 0.0, 1.0))
        return float(np.cos(np.radians(theta_deg))**2)

    # ── specific impulse ──────────────────────────────────────────────────────

    def isp(self, prop_name: Optional[str] = None,
            voltage: Optional[float] = None,
            B_field: Optional[float] = None) -> float:
        """Specific impulse [s]."""
        ve  = self.exhaust_velocity(prop_name, voltage)
        div = self.divergence_factor(B_field)
        return float(ve * div / G0)

    # ── mass flow rate ────────────────────────────────────────────────────────

    def mass_flow(self, prop_name: Optional[str] = None,
                  voltage: Optional[float] = None,
                  power: Optional[float]   = None) -> float:
        """Propellant mass flow rate [kg/s]."""
        p   = PROPELLANTS[prop_name] if prop_name else self.prop
        Vd  = voltage if voltage is not None else self.voltage
        Pd  = power   if power   is not None else self.power
        Id  = Pd / max(Vd, 1.0)
        eta = self.ioniz_eff(prop_name, voltage)
        return float(Id * p.ion_mass_kg / (E_C * eta))

    # ── thrust ────────────────────────────────────────────────────────────────

    def thrust(self, prop_name: Optional[str] = None,
               voltage: Optional[float] = None,
               power:   Optional[float] = None,
               B_field: Optional[float] = None) -> float:
        """Thrust [mN]."""
        mdot = self.mass_flow(prop_name, voltage, power)
        ve   = self.exhaust_velocity(prop_name, voltage)
        eta  = self.ioniz_eff(prop_name, voltage)
        div  = self.divergence_factor(B_field)
        T    = mdot * eta * ve * div
        return float(T * 1e3)

    # ── thrust efficiency ─────────────────────────────────────────────────────

    def thrust_efficiency(self, prop_name: Optional[str] = None,
                          voltage: Optional[float] = None,
                          power:   Optional[float] = None,
                          B_field: Optional[float] = None) -> float:
        """η_T = T² / (2 ṁ P)  [0–1]."""
        T    = self.thrust(prop_name, voltage, power, B_field) * 1e-3
        mdot = self.mass_flow(prop_name, voltage, power)
        Pd   = power if power is not None else self.power
        if Pd <= 0 or mdot <= 0:
            return 0.0
        return float(np.clip(T**2 / (2.0 * mdot * Pd), 0.0, 1.0))

    # ── thrust-to-power ───────────────────────────────────────────────────────

    def thrust_to_power(self, **kw) -> float:
        """T/P [mN/W]."""
        Pd = kw.get("power", self.power)
        return self.thrust(**kw) / max(Pd, 1.0)

    # ── electron Larmor radius ────────────────────────────────────────────────

    def larmor_radius(self, Te_eV: float = 15.0,
                      B_field: Optional[float] = None) -> float:
        """Electron Larmor radius [mm]."""
        B_T   = (B_field if B_field is not None else self.B_field) * 1e-4
        Te    = Te_eV * E_C / KB
        ve_th = np.sqrt(KB * Te / ME)
        return float(ME * ve_th / (E_C * max(B_T, 1e-6)) * 1e3)

    # ── ionization zone length ────────────────────────────────────────────────

    def ioniz_length(self) -> float:
        """Characteristic ionization zone length [mm]."""
        L   = self.channel_length
        B_T = self.B_field * 1e-4
        return float(np.clip(L * 0.3 / max(B_T / 0.01, 0.1), 0.5, L))

    # ── discharge current ─────────────────────────────────────────────────────

    def discharge_current(self) -> float:
        return self.power / max(self.voltage, 1.0)

    # ── power density ─────────────────────────────────────────────────────────

    def power_density(self) -> float:
        """W/cm²."""
        return self.power / max(self.channel_area_m2 * 1e4, 1e-6)

    # ── B-field axial profile ─────────────────────────────────────────────────

    def B_axial(self, n: int = 200) -> tuple:
        """
        Returns (x_mm, B_gauss) arrays for the chosen profile shape.
        """
        L   = self.channel_length
        B0  = self.B_field
        x   = np.linspace(0, L, n)
        mid = L * 0.45

        if self.B_profile == 0:
            B = np.full(n, B0)
        elif self.B_profile == 1:
            sigma = L * 0.20
            B = np.clip(B0 * np.exp(-((x - mid)**2) / (2 * sigma**2)), B0 * 0.08, B0)
        else:
            s1, s2   = L * 0.20, L * 0.25
            p1, p2   = L * 0.35, L * 0.65
            B = np.clip(
                0.70 * B0 * np.exp(-((x - p1)**2) / (2 * s1**2)) +
                0.50 * B0 * np.exp(-((x - p2)**2) / (2 * s2**2)),
                B0 * 0.04, B0
            )
        return x, B

    # ── voltage sweep ─────────────────────────────────────────────────────────

    def sweep_voltage(self, prop_name: Optional[str] = None,
                      V_range=(100, 800), n=100):
        pn  = prop_name or self.prop_name
        Vs  = np.linspace(*V_range, n)
        return Vs, {
            "isp":        np.array([self.isp(pn, v)              for v in Vs]),
            "thrust":     np.array([self.thrust(pn, v)           for v in Vs]),
            "efficiency": np.array([self.thrust_efficiency(pn, v)*100 for v in Vs]),
            "ve":         np.array([self.exhaust_velocity(pn, v)/1e3  for v in Vs]),
            "tp":         np.array([self.thrust_to_power(prop_name=pn, voltage=v) for v in Vs]),
        }

    # ── power sweep ───────────────────────────────────────────────────────────

    def sweep_power(self, prop_name: Optional[str] = None,
                    P_range=(200, 10000), n=100):
        pn  = prop_name or self.prop_name
        Ps  = np.linspace(*P_range, n)
        return Ps, {
            "isp":      np.array([self.isp(pn)                              for _ in Ps]),
            "thrust":   np.array([self.thrust(pn, power=p)                  for p in Ps]),
            "mflow":    np.array([self.mass_flow(pn, power=p)*1e6           for p in Ps]),
            "efficiency": np.array([self.thrust_efficiency(pn, power=p)*100 for p in Ps]),
        }

    # ── B-field sweep ─────────────────────────────────────────────────────────

    def sweep_B(self, prop_name: Optional[str] = None,
                B_range=(50, 400), n=100):
        pn = prop_name or self.prop_name
        Bs = np.linspace(*B_range, n)
        return Bs, {
            "isp":        np.array([self.isp(pn, B_field=b)               for b in Bs]),
            "thrust":     np.array([self.thrust(pn, B_field=b)            for b in Bs]),
            "efficiency": np.array([self.thrust_efficiency(pn, B_field=b)*100 for b in Bs]),
            "larmor":     np.array([self.larmor_radius(B_field=b)         for b in Bs]),
        }

    # ── efficiency heatmap ────────────────────────────────────────────────────

    def efficiency_map(self, prop_name: Optional[str] = None,
                       V_range=(100, 700), B_range=(50, 400), n=40):
        pn  = prop_name or self.prop_name
        Vs  = np.linspace(*V_range, n)
        Bs  = np.linspace(*B_range, n)
        VV, BB = np.meshgrid(Vs, Bs)
        ZZ = np.vectorize(
            lambda v, b: self.thrust_efficiency(pn, voltage=v, B_field=b) * 100
        )(VV, BB)
        return VV, BB, ZZ

    # ── full performance summary ──────────────────────────────────────────────

    def summary(self) -> dict:
        p   = self.prop
        mdot = self.mass_flow()
        return {
            "Propellant":         f"{p.name}  ({p.prop_type})",
            "Discharge Voltage":  f"{self.voltage:.0f} V",
            "Discharge Power":    f"{self.power:.0f} W",
            "Discharge Current":  f"{self.discharge_current():.2f} A",
            "B-Field":            f"{self.B_field:.0f} G",
            "Thrust":             f"{self.thrust():.3f} mN",
            "Specific Impulse":   f"{self.isp():.0f} s",
            "Exhaust Velocity":   f"{self.exhaust_velocity()/1e3:.2f} km/s",
            "Thrust Efficiency":  f"{self.thrust_efficiency()*100:.1f} %",
            "Ioniz. Efficiency":  f"{self.ioniz_eff()*100:.1f} %",
            "Mass Flow Rate":     f"{mdot*1e6:.4f} mg/s",
            "Thrust / Power":     f"{self.thrust_to_power():.5f} mN/W",
            "Power Density":      f"{self.power_density():.1f} W/cm²",
            "Accel. Voltage":     f"{self.accel_voltage():.1f} V",
            "e⁻ Larmor Radius":   f"{self.larmor_radius():.2f} mm",
            "Ioniz. Length":      f"{self.ioniz_length():.1f} mm",
            "Channel Area":       f"{self.channel_area_m2*1e4:.2f} cm²",
            "Double-Stage":       "ON" if self.double_stage else "OFF",
        }


# ══════════════════════════════════════════════════════════════════════════════
# COLOUR THEME
# ══════════════════════════════════════════════════════════════════════════════

BG       = "#0b0f1e"
PANEL    = "#131929"
BORDER   = "#1e2d44"
ACCENT   = "#4EA8DE"
ACCENT2  = "#F4A261"
SUCCESS  = "#52B788"
WARNING  = "#F4D35E"
DANGER   = "#E07A5F"
TEXT     = "#cdd6f4"
MUTED    = "#5a7a9a"
GRID_C   = "#1a2a3d"

def apply_theme():
    plt.rcParams.update({
        "figure.facecolor":       BG,
        "axes.facecolor":         PANEL,
        "axes.edgecolor":         BORDER,
        "axes.labelcolor":        MUTED,
        "axes.titlecolor":        ACCENT,
        "axes.titlesize":         9,
        "axes.labelsize":         8,
        "axes.titlepad":          5,
        "xtick.color":            MUTED,
        "ytick.color":            MUTED,
        "xtick.labelsize":        7,
        "ytick.labelsize":        7,
        "grid.color":             GRID_C,
        "grid.alpha":             1.0,
        "grid.linewidth":         0.5,
        "lines.linewidth":        2.0,
        "lines.solid_capstyle":   "round",
        "text.color":             TEXT,
        "font.family":            "monospace",
        "figure.dpi":             110,
        "savefig.facecolor":      BG,
    })


# ══════════════════════════════════════════════════════════════════════════════
# DASHBOARD
# ══════════════════════════════════════════════════════════════════════════════

class Dashboard:
    """
    Interactive matplotlib dashboard.

    Layout
    ──────
    Left strip  : controls (propellant, sliders, options)
    Right area  : 3 × 3 plot grid  +  schematic  +  summary
    """

    _SLIDER_KW = dict(color=ACCENT, track_color=BORDER, handle_style={"facecolor": ACCENT})

    def __init__(self):
        apply_theme()
        self.phys = HETPhysics()
        self._overlay_all  = False
        self._active_tab   = "sweep"   # 'sweep' | 'compare' | 'map' | 'schematic'
        self._build_layout()
        self._build_controls()
        self._connect_events()
        self._full_redraw()

    # ── layout ────────────────────────────────────────────────────────────────

    def _build_layout(self):
        self.fig = plt.figure(figsize=(22, 13))
        self.fig.canvas.manager.set_window_title(
            "Hall Effect Thruster  ·  Analysis System"
        )

        # Banner
        self.fig.text(0.5, 0.985,
                      "⚡  HALL EFFECT THRUSTER  ·  ANALYSIS SYSTEM",
                      ha="center", va="top", fontsize=14,
                      color=ACCENT, fontweight="bold")
        self.fig.text(0.5, 0.960,
                      "Electric Propulsion Performance Simulator  "
                      "·  Xe · Kr · Ar · N₂ · O₂ · H₂O",
                      ha="center", va="top", fontsize=8.5, color=MUTED)

        outer = gridspec.GridSpec(
            1, 2, figure=self.fig,
            left=0.01, right=0.995, top=0.950, bottom=0.025,
            wspace=0.025, width_ratios=[0.225, 0.775],
        )

        # ── Left panel: controls ──────────────────────────────────────────────
        ctrl = gridspec.GridSpecFromSubplotSpec(
            16, 1, subplot_spec=outer[0], hspace=0.45
        )
        self.ax_prop      = self.fig.add_subplot(ctrl[0:3])
        self.ax_b_prof    = self.fig.add_subplot(ctrl[3:5])
        self.ax_sl_v      = self.fig.add_subplot(ctrl[5])
        self.ax_sl_p      = self.fig.add_subplot(ctrl[6])
        self.ax_sl_b      = self.fig.add_subplot(ctrl[7])
        self.ax_sl_cl     = self.fig.add_subplot(ctrl[8])
        self.ax_sl_cw     = self.fig.add_subplot(ctrl[9])
        self.ax_opts      = self.fig.add_subplot(ctrl[10:12])
        self.ax_gauges    = self.fig.add_subplot(ctrl[12:16])

        for ax in [self.ax_prop, self.ax_b_prof, self.ax_sl_v, self.ax_sl_p,
                   self.ax_sl_b, self.ax_sl_cl, self.ax_sl_cw,
                   self.ax_opts, self.ax_gauges]:
            ax.set_facecolor(PANEL)
            for sp in ax.spines.values():
                sp.set_edgecolor(BORDER)

        # ── Right panel: plot grid ────────────────────────────────────────────
        plots = gridspec.GridSpecFromSubplotSpec(
            3, 3, subplot_spec=outer[1],
            hspace=0.52, wspace=0.32,
        )
        self.axs = np.empty((3, 3), dtype=object)
        for r in range(3):
            for c in range(3):
                ax = self.fig.add_subplot(plots[r, c])
                ax.set_facecolor(PANEL)
                ax.grid(True, lw=0.4, alpha=0.9)
                for sp in ax.spines.values():
                    sp.set_edgecolor(BORDER)
                self.axs[r, c] = ax

    # ── controls ──────────────────────────────────────────────────────────────

    def _build_controls(self):
        p = self.phys

        # ── Propellant radio ──────────────────────────────────────────────────
        self.ax_prop.set_title("PROPELLANT", fontsize=8, color=ACCENT, pad=3)
        self._radio_prop = RadioButtons(
            self.ax_prop, PROP_NAMES, active=0, activecolor=ACCENT
        )
        for lbl in self._radio_prop.labels:
            name = lbl.get_text()
            lbl.set_color(PROPELLANTS[name].color)
            lbl.set_fontsize(8.5)

        # ── B-field profile radio ─────────────────────────────────────────────
        self.ax_b_prof.set_title("B-FIELD PROFILE", fontsize=8, color=ACCENT, pad=3)
        self._radio_bp = RadioButtons(
            self.ax_b_prof,
            ["Uniform", "Peaked (Gaussian)", "Double-Peak"],
            active=0, activecolor=ACCENT2
        )
        for lbl in self._radio_bp.labels:
            lbl.set_color(TEXT); lbl.set_fontsize(8)

        # ── Sliders ───────────────────────────────────────────────────────────
        kw = self._SLIDER_KW
        self._sl_v  = Slider(self.ax_sl_v,  "Voltage  (V)",  100, 800,  valinit=p.voltage,        valstep=10,  **kw)
        self._sl_p  = Slider(self.ax_sl_p,  "Power    (W)",  200, 10000,valinit=p.power,           valstep=100, **kw)
        self._sl_b  = Slider(self.ax_sl_b,  "B-Field  (G)",   50, 400,  valinit=p.B_field,         valstep=5,   **kw)
        self._sl_cl = Slider(self.ax_sl_cl, "Ch. Len  (mm)",   5,  80,  valinit=p.channel_length,  valstep=1,   **kw)
        self._sl_cw = Slider(self.ax_sl_cw, "Ch. Wid  (mm)",   5,  60,  valinit=p.channel_width,   valstep=1,   **kw)

        for sl in [self._sl_v, self._sl_p, self._sl_b, self._sl_cl, self._sl_cw]:
            sl.label.set_color(TEXT);    sl.label.set_fontsize(8)
            sl.valtext.set_color(ACCENT);sl.valtext.set_fontsize(8)

        # ── Options check ─────────────────────────────────────────────────────
        self.ax_opts.set_title("OPTIONS", fontsize=8, color=ACCENT, pad=3)
        self._chk = CheckButtons(
            self.ax_opts,
            ["Double-Stage Accel.", "Overlay All Propellants"],
            actives=[False, False]
        )
        for lbl in self._chk.labels:
            lbl.set_color(TEXT); lbl.set_fontsize(8)

        # ── Performance gauges ────────────────────────────────────────────────
        self.ax_gauges.set_title("PERFORMANCE GAUGES", fontsize=8, color=ACCENT, pad=3)
        self.ax_gauges.axis("off")

    # ── event wiring ──────────────────────────────────────────────────────────

    def _connect_events(self):
        self._radio_prop.on_clicked(self._on_prop)
        self._radio_bp.on_clicked(self._on_bp)
        self._sl_v.on_changed(self._on_slider)
        self._sl_p.on_changed(self._on_slider)
        self._sl_b.on_changed(self._on_slider)
        self._sl_cl.on_changed(self._on_slider)
        self._sl_cw.on_changed(self._on_slider)
        self._chk.on_clicked(self._on_check)

    def _on_prop(self, label):
        self.phys.prop_name = label
        self._full_redraw()

    def _on_bp(self, label):
        mapping = {"Uniform": 0, "Peaked (Gaussian)": 1, "Double-Peak": 2}
        self.phys.B_profile = mapping[label]
        self._full_redraw()

    def _on_slider(self, _val):
        p = self.phys
        p.voltage        = self._sl_v.val
        p.power          = self._sl_p.val
        p.B_field        = self._sl_b.val
        p.channel_length = self._sl_cl.val
        p.channel_width  = self._sl_cw.val
        self._full_redraw()

    def _on_check(self, label):
        statuses = self._chk.get_status()
        self.phys.double_stage = statuses[0]
        self._overlay_all      = statuses[1]
        self._full_redraw()

    # ══════════════════════════════════════════════════════════════════════════
    # DRAWING HELPERS
    # ══════════════════════════════════════════════════════════════════════════

    def _clear_ax(self, ax):
        ax.cla()
        ax.set_facecolor(PANEL)
        ax.grid(True, lw=0.4, alpha=0.9)
        for sp in ax.spines.values():
            sp.set_edgecolor(BORDER)

    def _title(self, ax, t, xl, yl):
        ax.set_title(t,  fontsize=8,  color=ACCENT,  pad=4)
        ax.set_xlabel(xl, fontsize=7, color=MUTED)
        ax.set_ylabel(yl, fontsize=7, color=MUTED)

    # ══════════════════════════════════════════════════════════════════════════
    # PLOTS — ROW 0  (Voltage sweeps)
    # ══════════════════════════════════════════════════════════════════════════

    def _plot_isp_v(self):
        """[0,0] Isp vs Discharge Voltage."""
        ax = self.axs[0, 0]
        self._clear_ax(ax)
        names = PROP_NAMES if self._overlay_all else [self.phys.prop_name]
        for name in names:
            Vs, res = self.phys.sweep_voltage(name)
            ax.plot(Vs, res["isp"], color=PROPELLANTS[name].color,
                    lw=2.0 if name == self.phys.prop_name else 1.0,
                    alpha=1.0 if name == self.phys.prop_name else 0.55,
                    label=PROPELLANTS[name].symbol)
        ax.axvline(self.phys.voltage, color=SUCCESS, lw=1.0, ls="--", alpha=0.7)
        if self._overlay_all:
            ax.legend(fontsize=7, facecolor=PANEL, edgecolor=BORDER, loc="upper left")
        self._title(ax, f"Isp vs Voltage  [{PROPELLANTS[self.phys.prop_name].symbol}]",
                    "Discharge Voltage (V)", "Isp (s)")

    def _plot_thrust_v(self):
        """[0,1] Thrust vs Discharge Voltage."""
        ax = self.axs[0, 1]
        self._clear_ax(ax)
        names = PROP_NAMES if self._overlay_all else [self.phys.prop_name]
        for name in names:
            Vs, res = self.phys.sweep_voltage(name)
            ax.plot(Vs, res["thrust"], color=PROPELLANTS[name].color,
                    lw=2.0 if name == self.phys.prop_name else 1.0,
                    alpha=1.0 if name == self.phys.prop_name else 0.55,
                    label=PROPELLANTS[name].symbol)
        ax.axvline(self.phys.voltage, color=SUCCESS, lw=1.0, ls="--", alpha=0.7)
        if self._overlay_all:
            ax.legend(fontsize=7, facecolor=PANEL, edgecolor=BORDER)
        self._title(ax, "Thrust vs Voltage", "Discharge Voltage (V)", "Thrust (mN)")

    def _plot_eff_v(self):
        """[0,2] Thrust Efficiency vs Discharge Voltage."""
        ax = self.axs[0, 2]
        self._clear_ax(ax)
        names = PROP_NAMES if self._overlay_all else [self.phys.prop_name]
        for name in names:
            Vs, res = self.phys.sweep_voltage(name)
            ax.plot(Vs, res["efficiency"], color=PROPELLANTS[name].color,
                    lw=2.0 if name == self.phys.prop_name else 1.0,
                    alpha=1.0 if name == self.phys.prop_name else 0.55,
                    label=PROPELLANTS[name].symbol)
        ax.axvline(self.phys.voltage, color=SUCCESS, lw=1.0, ls="--", alpha=0.7)
        ax.set_ylim(0, 80)
        if self._overlay_all:
            ax.legend(fontsize=7, facecolor=PANEL, edgecolor=BORDER)
        self._title(ax, "Thrust Efficiency vs Voltage",
                    "Discharge Voltage (V)", "η_T (%)")

    # ══════════════════════════════════════════════════════════════════════════
    # PLOTS — ROW 1
    # ══════════════════════════════════════════════════════════════════════════

    def _plot_thrust_power(self):
        """[1,0] Thrust vs Power."""
        ax = self.axs[1, 0]
        self._clear_ax(ax)
        names = PROP_NAMES if self._overlay_all else [self.phys.prop_name]
        for name in names:
            Ps, res = self.phys.sweep_power(name)
            ax.plot(Ps / 1000, res["thrust"], color=PROPELLANTS[name].color,
                    lw=2.0 if name == self.phys.prop_name else 1.0,
                    alpha=1.0 if name == self.phys.prop_name else 0.55,
                    label=PROPELLANTS[name].symbol)
        ax.axvline(self.phys.power / 1000, color=SUCCESS, lw=1.0, ls="--", alpha=0.7)
        if self._overlay_all:
            ax.legend(fontsize=7, facecolor=PANEL, edgecolor=BORDER)
        self._title(ax, "Thrust vs Power", "Power (kW)", "Thrust (mN)")

    def _plot_eff_B(self):
        """[1,1] Efficiency vs B-field  + Larmor radius twin-axis."""
        ax = self.axs[1, 1]
        self._clear_ax(ax)
        Bs, res = self.phys.sweep_B()
        ax.plot(Bs, res["efficiency"], color=ACCENT, lw=2.0, label="η_T (%)")
        ax.axvline(self.phys.B_field, color=SUCCESS, lw=1.0, ls="--", alpha=0.7)

        ax2 = ax.twinx()
        ax2.set_facecolor(PANEL)
        ax2.plot(Bs, res["larmor"], color=WARNING, lw=1.5, ls="--", label="r_L (mm)")
        ax2.set_ylabel("e⁻ Larmor radius (mm)", fontsize=7, color=WARNING)
        ax2.tick_params(colors=MUTED, labelsize=7)
        ax2.set_facecolor("none")

        self._title(ax, "Efficiency + Larmor Radius vs B-Field",
                    "B-Field (Gauss)", "η_T (%)")
        lines  = ax.get_lines() + ax2.get_lines()
        labels = [l.get_label() for l in lines]
        ax.legend(lines, labels, fontsize=7, facecolor=PANEL, edgecolor=BORDER)

    def _plot_mflow_power(self):
        """[1,2] Mass flow + exhaust velocity vs Power."""
        ax = self.axs[1, 2]
        self._clear_ax(ax)
        Ps, res = self.phys.sweep_power()
        ax.plot(Ps / 1000, res["mflow"], color=ACCENT2, lw=2.0, label="ṁ (mg/s)")
        ax.axvline(self.phys.power / 1000, color=SUCCESS, lw=1.0, ls="--", alpha=0.7)

        ax2 = ax.twinx()
        ax2.set_facecolor("none")
        ax2.tick_params(colors=MUTED, labelsize=7)
        ax2.plot(Ps / 1000, res["efficiency"], color=ACCENT, lw=1.5, ls="--",
                 label="η_T (%)")
        ax2.set_ylabel("Thrust Efficiency (%)", fontsize=7, color=ACCENT)

        lines  = ax.get_lines() + ax2.get_lines()
        labels = [l.get_label() for l in lines]
        ax.legend(lines, labels, fontsize=7, facecolor=PANEL, edgecolor=BORDER)
        self._title(ax, "Mass Flow + Efficiency vs Power",
                    "Power (kW)", "ṁ (mg/s)")

    # ══════════════════════════════════════════════════════════════════════════
    # PLOTS — ROW 2
    # ══════════════════════════════════════════════════════════════════════════

    def _plot_propellant_compare(self):
        """[2,0] Normalised propellant comparison bar chart."""
        ax = self.axs[2, 0]
        self._clear_ax(ax)

        metrics = {}
        for name in PROP_NAMES:
            p = self.phys
            metrics[name] = {
                "isp":    p.isp(name),
                "thrust": p.thrust(name),
                "eta":    p.thrust_efficiency(name) * 100,
            }

        max_isp    = max(v["isp"]    for v in metrics.values()) or 1
        max_thrust = max(v["thrust"] for v in metrics.values()) or 1
        max_eta    = max(v["eta"]    for v in metrics.values()) or 1

        x  = np.arange(len(PROP_NAMES))
        w  = 0.26
        for i, (key, label, norm) in enumerate([
            ("isp",    "Isp",    max_isp),
            ("thrust", "Thrust", max_thrust),
            ("eta",    "η_T",    max_eta),
        ]):
            vals   = [metrics[n][key] / norm * 100 for n in PROP_NAMES]
            colors = [PROPELLANTS[n].color for n in PROP_NAMES]
            alpha  = 1.0 - i * 0.25
            ax.bar(x + (i - 1) * w, vals, w,
                   color=colors, alpha=alpha,
                   edgecolor=BORDER, linewidth=0.5, label=label)

        ax.axvline(2.5, color=MUTED, lw=0.8, ls=":", alpha=0.7)
        ax.text(1.0, 105, "monatomic", ha="center", fontsize=7, color=MUTED)
        ax.text(4.0, 105, "molecular", ha="center", fontsize=7, color=MUTED)

        ax.set_xticks(x)
        ax.set_xticklabels([PROPELLANTS[n].symbol for n in PROP_NAMES], fontsize=8)
        ax.set_ylim(0, 115)
        # Highlight active propellant
        idx = PROP_NAMES.index(self.phys.prop_name)
        ax.axvspan(idx - 0.45, idx + 0.45, alpha=0.10, color=SUCCESS, zorder=0)
        ax.legend(fontsize=7, facecolor=PANEL, edgecolor=BORDER)
        self._title(ax, "Propellant Comparison (Normalised)",
                    "Propellant", "Relative Performance (%)")

    def _plot_isp_eff_scatter(self):
        """[2,1] Isp vs efficiency scatter across all propellants (V sweep)."""
        ax = self.axs[2, 1]
        self._clear_ax(ax)
        Vs = np.linspace(150, 700, 25)
        for name in PROP_NAMES:
            isps  = [self.phys.isp(name, v)               for v in Vs]
            effs  = [self.phys.thrust_efficiency(name, v) * 100 for v in Vs]
            marker = "o" if PROPELLANTS[name].prop_type == "monatomic" else "s"
            ax.scatter(isps, effs, c=PROPELLANTS[name].color,
                       s=18, alpha=0.55, marker=marker,
                       label=PROPELLANTS[name].symbol)
        # Operating point star
        ax.scatter([self.phys.isp()], [self.phys.thrust_efficiency() * 100],
                   c=SUCCESS, s=140, marker="*", zorder=6, label="Op. Pt.")
        ax.legend(fontsize=7, facecolor=PANEL, edgecolor=BORDER, ncol=2)
        ax.text(0.98, 0.04, "○=mono  □=mol.", transform=ax.transAxes,
                fontsize=6.5, color=MUTED, ha="right")
        self._title(ax, "Isp vs Efficiency  (all propellants, V sweep)",
                    "Specific Impulse (s)", "Thrust Efficiency (%)")

    def _plot_efficiency_map(self):
        """[2,2] 2-D efficiency heatmap over voltage × B-field."""
        ax = self.axs[2, 2]
        self._clear_ax(ax)
        VV, BB, ZZ = self.phys.efficiency_map(n=35)
        cf = ax.contourf(VV, BB, ZZ, levels=22,
                         cmap=plt.cm.plasma, alpha=0.88)
        self.fig.colorbar(cf, ax=ax, fraction=0.038, pad=0.02,
                          label="η_T (%)").ax.tick_params(labelsize=6, colors=MUTED)
        # Mark current operating point
        ax.scatter([self.phys.voltage], [self.phys.B_field],
                   s=120, c=SUCCESS, marker="*", zorder=5)
        ax.text(self.phys.voltage + 12, self.phys.B_field + 8,
                "▲ OP", fontsize=7, color=SUCCESS)
        self._title(ax, "Efficiency Map  (V × B  sweep)",
                    "Discharge Voltage (V)", "B-Field (Gauss)")

    # ══════════════════════════════════════════════════════════════════════════
    # B-FIELD PROFILE SCHEMATIC  (replaces [0,0] when viewing schematic tab)
    # Here always drawn inline in gauges + dedicated axes row
    # ══════════════════════════════════════════════════════════════════════════

    def _plot_bfield_profile(self, ax):
        """Draw the axial B-field profile in a given axes."""
        ax.cla()
        ax.set_facecolor(PANEL)
        ax.grid(True, lw=0.4, alpha=0.9)
        x, B = self.phys.B_axial()
        ax.fill_between(x, 0, B, alpha=0.20, color=ACCENT)
        ax.plot(x, B, color=ACCENT, lw=1.8)
        iz = self.phys.ioniz_length()
        ax.axvspan(0, iz, alpha=0.15, color=ACCENT2, label=f"Ioniz. zone ~{iz:.1f}mm")
        ax.axvline(iz, color=ACCENT2, lw=1.0, ls=":")
        ax.set_ylim(0, self.phys.B_field * 1.2)
        ax.legend(fontsize=6.5, facecolor=PANEL, edgecolor=BORDER)
        self._title(ax, "B-Field Axial Profile", "Axial Position (mm)", "B (Gauss)")

    # ══════════════════════════════════════════════════════════════════════════
    # THRUSTER SCHEMATIC
    # ══════════════════════════════════════════════════════════════════════════

    def _draw_schematic(self, ax):
        """Render a labelled thruster cross-section in axes space [0,1]×[0,1]."""
        ax.cla()
        ax.set_facecolor("#060a14")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.set_title("Thruster Cross-Section Schematic", fontsize=9, color=ACCENT, pad=4)

        p   = self.phys
        Lf  = np.clip(p.channel_length / 80, 0.12, 0.55)   # fraction of axes width
        x0  = 0.12
        x1  = x0 + Lf
        cy  = 0.50
        hw  = np.clip(p.channel_width / 60, 0.08, 0.28)    # half-width fraction

        # ── walls ──────────────────────────────────────────────────────────
        wall_h = 0.06
        for y_edge, sign in [(cy + hw, 1), (cy - hw - wall_h, -1)]:
            ax.add_patch(FancyBboxPatch(
                (x0, y_edge), Lf, wall_h,
                boxstyle="square,pad=0", linewidth=0.8,
                edgecolor=ACCENT, facecolor="#1a3d5c"
            ))

        # ── channel ──────────────────────────────────────────────────────
        ax.add_patch(Rectangle(
            (x0, cy - hw), Lf, 2 * hw,
            linewidth=1.0, edgecolor=ACCENT, facecolor="#0d2a3d"
        ))

        # ── radial B-field lines ──────────────────────────────────────────
        for bx in np.linspace(x0 + Lf * 0.15, x0 + Lf * 0.90, 6):
            ys = np.linspace(cy - hw - wall_h, cy + hw + wall_h, 40)
            xs = bx + 0.008 * np.sin(np.linspace(0, np.pi, 40))
            ax.plot(xs, ys, color=WARNING, alpha=0.55, lw=0.9, ls=":")

        # ── plasma plume ─────────────────────────────────────────────────
        for i in range(35):
            frac  = i / 34
            px    = x1 + frac * 0.30
            spread = hw * (1 + frac * 1.8)
            alpha  = 0.22 * (1 - frac)**1.8
            ax.plot([px, px], [cy - spread, cy + spread],
                    color=ACCENT, alpha=alpha, lw=3.5)

        # ── ion trajectories ──────────────────────────────────────────────
        rng = np.random.default_rng(42)
        for y_start in np.linspace(cy - hw + 0.015, cy + hw - 0.015, 5):
            xs = np.linspace(x0 + Lf * 0.45, x1 + 0.22, 30)
            ys = y_start + 0.004 * rng.standard_normal(30).cumsum()
            ax.plot(xs, ys, color=ACCENT, alpha=0.35, lw=0.7)

        # ── labels ───────────────────────────────────────────────────────
        kw = dict(transform=ax.transAxes, fontsize=7.5)
        ax.text(x0 + Lf / 2, cy + hw + wall_h + 0.05, "OUTER WALL",
                ha="center", color=MUTED, **kw)
        ax.text(x0 + Lf / 2, cy - hw - wall_h - 0.05, "INNER WALL",
                ha="center", color=MUTED, **kw)
        ax.text(x0 - 0.03, cy, "← ANODE", ha="right",
                color=SUCCESS, **kw)
        ax.text(x1 + 0.16, cy, "PLUME →", ha="center",
                color=ACCENT2, **kw)
        ax.text(x0 + Lf / 2, cy + hw + wall_h + 0.12,
                f"B⃗ radial  ({p.B_field:.0f} G)  ·  {p.channel_length:.0f} mm channel",
                ha="center", color=WARNING, fontsize=7, **kw)
        ax.text(x0 + Lf / 2, 0.07,
                f"{PROPELLANTS[p.prop_name].symbol}  ·  "
                f"Vd={p.voltage:.0f}V  ·  P={p.power:.0f}W  ·  "
                f"Isp={p.isp():.0f}s  ·  T={p.thrust():.2f}mN",
                ha="center", color=MUTED, fontsize=7, **kw)

    # ══════════════════════════════════════════════════════════════════════════
    # PERFORMANCE GAUGES (left-panel)
    # ══════════════════════════════════════════════════════════════════════════

    def _draw_gauges(self):
        ax = self.ax_gauges
        ax.cla()
        ax.set_facecolor(PANEL)
        ax.axis("off")
        ax.set_title("PERFORMANCE GAUGES", fontsize=8, color=ACCENT, pad=3)

        phys = self.phys
        rows = [
            ("Thrust",  phys.thrust(),              500,  "mN"),
            ("Isp",     phys.isp(),                 3500, "s"),
            ("η_T",     phys.thrust_efficiency()*100, 70, "%"),
            ("η_i",     phys.ioniz_eff()*100,        95, "%"),
            ("T/P",     phys.thrust_to_power()*1000, 100, "μN/W"),
            ("Ve",      phys.exhaust_velocity()/1e3, 50,  "km/s"),
        ]

        n   = len(rows)
        bar_h = 0.09
        gap   = 0.14

        for i, (label, val, ref, unit) in enumerate(rows):
            frac = np.clip(val / max(ref, 1e-9), 0, 1)
            y    = 0.92 - i * gap
            # background track
            ax.add_patch(Rectangle((0.30, y), 0.65, bar_h,
                                   transform=ax.transAxes,
                                   color=BORDER, zorder=1))
            # filled bar
            bar_color = (SUCCESS if frac > 0.65 else
                         WARNING if frac > 0.35 else DANGER)
            ax.add_patch(Rectangle((0.30, y), 0.65 * frac, bar_h,
                                   transform=ax.transAxes,
                                   color=bar_color, alpha=0.85, zorder=2))
            # label
            ax.text(0.28, y + bar_h / 2, label,
                    transform=ax.transAxes, fontsize=7,
                    color=MUTED, ha="right", va="center")
            # value
            if unit == "mN":
                disp = f"{val:.2f} mN"
            elif unit == "s":
                disp = f"{val:.0f} s"
            elif unit == "%":
                disp = f"{val:.1f}%"
            elif unit == "km/s":
                disp = f"{val:.1f} km/s"
            elif unit == "μN/W":
                disp = f"{val:.2f} μN/W"
            else:
                disp = f"{val:.3g} {unit}"
            ax.text(0.97, y + bar_h / 2, disp,
                    transform=ax.transAxes, fontsize=7,
                    color=bar_color, ha="right", va="center", fontweight="bold")

    # ══════════════════════════════════════════════════════════════════════════
    # SUMMARY TABLE  (printed to terminal + drawn in last row if schematic tab)
    # ══════════════════════════════════════════════════════════════════════════

    def _draw_summary_table(self, ax):
        """Draw the summary key-value table in the given axes."""
        ax.cla()
        ax.set_facecolor(PANEL)
        ax.axis("off")
        ax.set_title("PERFORMANCE SUMMARY", fontsize=8, color=ACCENT, pad=3)

        summary = self.phys.summary()
        items   = list(summary.items())
        n       = len(items)
        for i, (k, v) in enumerate(items):
            y = 0.97 - i * (0.92 / n)
            ax.text(0.02, y, k + ":", fontsize=7, color=MUTED,
                    transform=ax.transAxes, va="top")
            ax.text(0.98, y, v, fontsize=7,
                    color=ACCENT if i < 3 else TEXT,
                    transform=ax.transAxes, va="top", ha="right",
                    fontweight="bold" if i < 3 else "normal")

    # ══════════════════════════════════════════════════════════════════════════
    # MASTER REDRAW
    # ══════════════════════════════════════════════════════════════════════════

    def _full_redraw(self):
        # Row 0
        self._plot_isp_v()
        self._plot_thrust_v()
        self._plot_eff_v()
        # Row 1
        self._plot_thrust_power()
        self._plot_eff_B()
        self._plot_mflow_power()
        # Row 2
        self._plot_propellant_compare()
        self._plot_isp_eff_scatter()
        self._plot_efficiency_map()

        # Schematic in a special axes: reuse axs[1,1] twin or draw over it
        # We dedicate axs[2,1] for scatter and embed schematic in the
        # middle of the row. Instead: append a separate axes under the grid.
        # For clean layout we place the schematic text in ax_gauges area
        # and draw the thruster in a fresh inset_axes of the figure.
        self._draw_bfield_inset()
        self._draw_gauges()
        self.fig.canvas.draw_idle()

    def _draw_bfield_inset(self):
        """Draw the B-field profile directly in axs[1,1] secondary view."""
        # We use axs[1,1] for eff-B sweep already.
        # Instead we insert a small thruster schematic as a figure-level inset.
        if hasattr(self, "_schem_ax") and self._schem_ax in self.fig.axes:
            self._schem_ax.remove()
        self._schem_ax = self.fig.add_axes([0.27, 0.032, 0.26, 0.20])
        self._schem_ax.set_facecolor("#060a14")
        self._draw_schematic(self._schem_ax)

    # ── public entry point ────────────────────────────────────────────────────

    def show(self):
        # Print summary to terminal
        print("\n" + "─" * 54)
        print(f"  Current Performance Summary — {self.phys.prop_name}")
        print("─" * 54)
        for k, v in self.phys.summary().items():
            print(f"  {k:<22}  {v}")
        print("─" * 54)
        print("  Use the sliders and radio buttons to explore.\n")
        plt.show()


# ══════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

def main():
    print("╔══════════════════════════════════════════════════════╗")
    print("║   HALL EFFECT THRUSTER  ·  ANALYSIS SYSTEM  v3.0   ║")
    print("╠══════════════════════════════════════════════════════╣")
    print("║  Monatomic : Xenon · Krypton · Argon                ║")
    print("║  Molecular : Nitrogen · Oxygen · Water              ║")
    print("╠══════════════════════════════════════════════════════╣")
    print("║  Controls  : propellant selector  ·  6 sliders      ║")
    print("║              B-field profile  ·  double-stage mode  ║")
    print("╚══════════════════════════════════════════════════════╝")
    print()

    dash = Dashboard()
    dash.show()


if __name__ == "__main__":
    main()
