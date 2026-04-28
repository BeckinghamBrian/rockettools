# =============================================================================
# mat_dict.py
# Thermomechanical property dictionary for solid materials.
#
# STRUCTURE
# ---------
# Top-level key  : material name string  (e.g. "Inconel_718")
# Second-level   : temperature in Kelvin (float)  — single state variable
# Value          : dict of property name → float or None
#
# INCLUDED PROPERTIES
# -------------------
# Required (used directly in Engine_Class.heatTransfer):
#   density_kg_m3              kg/m³
#   thermal_conductivity_W_mK  W/(m·K)
#   coef_thermal_expansion_1_K 1/K      (linear CTE)
#   youngs_modulus_GPa         GPa
#   yield_strength_MPa         MPa
#
# Suggested additions (not included — add per-project as needed):
#   ultimate_tensile_strength_MPa  — needed for burst/fracture margin calc
#   poissons_ratio                 — needed for full 3-D thermal stress tensor
#   specific_heat_J_kgK           — needed if transient (non-steady) heat soak is modelled
#   creep_rupture_stress_MPa       — critical for high-cycle fatigue / life prediction
#   fatigue_strength_MPa           — high-cycle fatigue limit at given T
#   oxidation_resistance_rating    — qualitative/quantitative for hot-gas wall selection
#   emissivity                     — needed for radiation boundary condition
#
# USAGE
# -----
#   from mat_dict import mat_dict, get_mat_props
#
#   # Exact O(1) lookup
#   props = mat_dict["Inconel_718"][973.15]
#   k     = props["thermal_conductivity_W_mK"]
#
#   # Nearest-temperature lookup with linear interpolation
#   props = get_mat_props("Inconel_718", T_K=900.0)
# =============================================================================

from __future__ import annotations
from typing import Dict, Any

_MatDict = Dict[str, Dict[float, Dict[str, Any]]]

# =============================================================================
# Main dictionary
# =============================================================================

mat_dict: _MatDict = {

    # ── Inconel 718 (UNS N07718) ──────────────────────────────────────────────
    # Primary nickel-based superalloy used for rocket engine chamber walls.
    # Data sourced from:
    #   • Special Metals Corp. data sheet (2007)
    #   • ASM Aerospace Specification Metals
    #   • NASA SP-8120 (material properties for liquid rocket engines)
    # All values are for wrought + precipitation-hardened condition unless noted.

    "Inconel_718": {

        # ── 25 °C (298.15 K) — room-temperature baseline ─────────────────────
        298.15: {
            "temperature_K":                298.15,
            "density_kg_m3":                8190.0,     # kg/m³
            "thermal_conductivity_W_mK":    11.4,       # W/(m·K)
            "coef_thermal_expansion_1_K":   1.28e-5,   # 1/K  (20–100 °C mean)
            "youngs_modulus_GPa":           200.0,      # GPa
            "yield_strength_MPa":           1100.0,     # MPa  (0.2% offset, annealed+aged)
            "notes": "Inconel 718 room-temperature baseline. "
                     "Special Metals data sheet Rev. Sept 2007.",
        },

        # ── 700 °C (973.15 K) — representative hot chamber-wall temperature ──
        973.15: {
            "temperature_K":                973.15,
            "density_kg_m3":               8060.0,     # kg/m³  (thermal expansion corrected)
            "thermal_conductivity_W_mK":   18.7,       # W/(m·K) — increases with T for IN718
            "coef_thermal_expansion_1_K":   1.46e-5,   # 1/K  (20–700 °C mean)
            "youngs_modulus_GPa":          156.0,      # GPa  — significant drop at elevated T
            "yield_strength_MPa":          860.0,      # MPa  (0.2% offset @ 700 °C)
            "notes": "Inconel 718 @ 700°C (973.15 K). "
                     "ASM Aerospace Spec Metals; NASA SP-8120.",
        },

        # ── 500 °C (773.15 K) ────────────────────────────────────────────────
        773.15: {
            "temperature_K":                773.15,
            "density_kg_m3":               8130.0,
            "thermal_conductivity_W_mK":   16.0,
            "coef_thermal_expansion_1_K":   1.40e-5,
            "youngs_modulus_GPa":          172.0,
            "yield_strength_MPa":          980.0,
            "notes": "Inconel 718 @ 500°C (773.15 K).",
        },

        # ── 800 °C (1073.15 K) ───────────────────────────────────────────────
        1073.15: {
            "temperature_K":                1073.15,
            "density_kg_m3":               8030.0,
            "thermal_conductivity_W_mK":   20.5,
            "coef_thermal_expansion_1_K":   1.52e-5,
            "youngs_modulus_GPa":          141.0,
            "yield_strength_MPa":          620.0,      # sharp drop approaching γ″ solvus
            "notes": "Inconel 718 @ 800°C (1073.15 K). Near precipitate dissolution range.",
        },
    },

    # ── Copper (C10100 / OFHC) ────────────────────────────────────────────────
    # Often used for inner chamber liner in regeneratively cooled engines.
    "Copper_OFHC": {

        298.15: {
            "temperature_K":                298.15,
            "density_kg_m3":               8960.0,
            "thermal_conductivity_W_mK":   391.0,
            "coef_thermal_expansion_1_K":   1.67e-5,
            "youngs_modulus_GPa":          117.0,
            "yield_strength_MPa":           70.0,      # annealed — very low, work-hardened >> higher
            "notes": "OFHC copper (C10100) room-temperature. CRC Handbook.",
        },

        573.15: {
            "temperature_K":                573.15,
            "density_kg_m3":               8900.0,
            "thermal_conductivity_W_mK":   374.0,
            "coef_thermal_expansion_1_K":   1.74e-5,
            "youngs_modulus_GPa":           98.0,
            "yield_strength_MPa":           50.0,
            "notes": "OFHC copper @ 300°C (573.15 K).",
        },
    },

    # ── 304 Stainless Steel ───────────────────────────────────────────────────
    "SS_304": {

        298.15: {
            "temperature_K":                298.15,
            "density_kg_m3":               8000.0,
            "thermal_conductivity_W_mK":   16.3,
            "coef_thermal_expansion_1_K":   1.72e-5,
            "youngs_modulus_GPa":          193.0,
            "yield_strength_MPa":          215.0,
            "notes": "304 SS room-temperature. AISI data.",
        },

        773.15: {
            "temperature_K":                773.15,
            "density_kg_m3":               7870.0,
            "thermal_conductivity_W_mK":   20.5,
            "coef_thermal_expansion_1_K":   1.84e-5,
            "youngs_modulus_GPa":          154.0,
            "yield_strength_MPa":          130.0,
            "notes": "304 SS @ 500°C (773.15 K).",
        },
    },
}

# =============================================================================
# Helper functions
# =============================================================================

def add_mat_state(material: str, T_K: float, props: dict) -> None:
    """
    Register a new (material, T) state point at runtime.

    Parameters
    ----------
    material : str     Material name key.
    T_K      : float   Temperature in Kelvin.
    props    : dict    Property name → value.

    Example
    -------
        add_mat_state("Inconel_718", 1173.15, {
            "thermal_conductivity_W_mK": 22.1,
            "yield_strength_MPa":        310.0,
            ...
        })
    """
    if material not in mat_dict:
        mat_dict[material] = {}
    mat_dict[material][T_K] = props


def get_mat_props(material: str, T_K: float,
                  interpolate: bool = True) -> dict | None:
    """
    Exact lookup with optional linear interpolation between registered
    temperature points.

    Parameters
    ----------
    material    : str    Material name key.
    T_K         : float  Desired temperature (K).
    interpolate : bool   If True and exact key absent, linearly interpolate
                         between the two nearest registered temperatures.
                         The numeric properties are interpolated; string fields
                         ('notes') are taken from the lower bracket.

    Returns None if the material is unknown or T_K is out of range.
    """
    states = mat_dict.get(material)
    if states is None:
        return None

    # O(1) exact hit
    exact = states.get(T_K)
    if exact is not None:
        return exact

    if not interpolate:
        return None

    # Linear interpolation between bracketing temperatures
    temps = sorted(states.keys())
    if T_K < temps[0] or T_K > temps[-1]:
        return None     # out of registered range

    # Find bracket  (O(n) sort already done; could cache sorted keys for O(log n))
    for i in range(len(temps) - 1):
        T_lo, T_hi = temps[i], temps[i + 1]
        if T_lo <= T_K <= T_hi:
            f  = (T_K - T_lo) / (T_hi - T_lo)   # interpolation fraction
            lo = states[T_lo]
            hi = states[T_hi]
            result = {"temperature_K": T_K}
            for key in lo:
                v_lo = lo[key]
                v_hi = hi.get(key, v_lo)
                if isinstance(v_lo, (int, float)) and isinstance(v_hi, (int, float)):
                    result[key] = v_lo + f * (v_hi - v_lo)
                else:
                    result[key] = v_lo   # non-numeric: keep lower-bracket value
            return result

    return None
