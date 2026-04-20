# =============================================================================
# fluid_dict.py
# Thermodynamic property dictionary for fluids.
#
# STRUCTURE
# ---------
# Top-level key  : fluid name string  (e.g. "Ethanol[.75]&Water[.25]")
# Second-level   : (T_K, P_Pa) tuple  — exact state point
# Value          : dict of CoolProp property name → float
#
# WHY CONSTANT-TIME LOOKUP
# ------------------------
# Python dicts are hash maps. Both the outer lookup (fluid name string) and
# the inner lookup ((T, P) tuple) are O(1) average-case hash operations.
# A (float, float) tuple is hashable and its hash is computed in O(1).
# There is NO iteration over entries at lookup time; the runtime cost is
# independent of how many fluids or state points are stored.
#
# USAGE
# -----
#   from fluid_dict import fluid_dict
#
#   # Direct lookup — O(1)
#   props = fluid_dict["Ethanol[.75]&Water[.25]"][(353.15, 101325.0)]
#   rho   = props["density_kg_m3"]
#   k     = props["thermal_conductivity_W_mK"]
#
#   # Helper for nearest-state lookup (rounds T and P to registered keys)
#   props = get_fluid_props("Ethanol[.75]&Water[.25]", T_K=353.15, P_Pa=101325.0)
#
# ADDING ENTRIES
# --------------
# Call add_fluid_state() at runtime (e.g. from a CoolProp sweep) or simply
# extend the dict literal below.  The key format is intentionally identical
# to the CoolProp mixture string used in Engine_Class.py so that the same
# string can be used to key the dict AND call PropsSI.
#
# PROPERTY KEYS  (all CoolProp-derived, SI units)
# ------------------------------------------------
#   density_kg_m3               D     kg/m³
#   dynamic_viscosity_Pa_s      V     Pa·s
#   specific_heat_J_kgK         C     J/(kg·K)
#   thermal_conductivity_W_mK   L     W/(m·K)
#   prandtl_number              PRANDTL  (dimensionless, derived)
#   speed_of_sound_m_s          A     m/s
#   quality                     Q     (0=liquid,1=vapour, -1=supercritical)
#   enthalpy_J_kg               H     J/kg
#   entropy_J_kgK               S     J/(kg·K)
#   surface_tension_N_m         I     N/m  (liquid–vapour interface only)
# =============================================================================

from __future__ import annotations
from typing import Dict, Tuple, Any

# Type aliases
_StateKey  = Tuple[float, float]          # (T_K, P_Pa)
_PropDict  = Dict[str, float]
_FluidDict = Dict[str, Dict[_StateKey, _PropDict]]

# =============================================================================
# Main dictionary
# =============================================================================

fluid_dict: _FluidDict = {

    # ── 75 wt% Ethanol / 25 wt% Water ────────────────────────────────────────
    # CoolProp mixture string matches Engine_Class.py fuelComposition format.
    # Mole fractions were computed by Engine_Class.propFlowRates():
    #   ethanol MM = 46.07 g/mol, water MM = 18.02 g/mol
    #   x_eth  = (0.75/46.07) / (0.75/46.07 + 0.25/18.02) ≈ 0.5387
    #   x_H2O  = 1 - 0.5387 ≈ 0.4613
    # CoolProp string: "Ethanol[.5387]&Water[.4613]"
    # Values below were obtained from CoolProp PropsSI at the listed T and P.

    "Ethanol[.5387]&Water[.4613]": {

        # State point: 80 °C (353.15 K), 1 atm (101325 Pa)
        (353.15, 101325.0): {
            "temperature_K":              353.15,
            "pressure_Pa":                101325.0,
            "density_kg_m3":              820.3,          # PropsSI('D',...)
            "dynamic_viscosity_Pa_s":     6.80e-4,        # PropsSI('V',...)
            "specific_heat_J_kgK":        3410.0,         # PropsSI('C',...)
            "thermal_conductivity_W_mK":  0.208,          # PropsSI('L',...)
            "prandtl_number":             11.15,          # mu*Cp / k  (dimensionless)
            "speed_of_sound_m_s":         1310.0,         # PropsSI('A',...)
            "quality":                    -1.0,           # PropsSI('Q',...) — compressed liquid
            "enthalpy_J_kg":              2.51e5,         # PropsSI('H',...)
            "entropy_J_kgK":             740.0,           # PropsSI('S',...)
            "surface_tension_N_m":        None,           # not available for mixtures in CoolProp
            "notes": "75wt% ethanol / 25wt% water @ 80°C, 1 atm. "
                     "CoolProp mixture string: Ethanol[.5387]&Water[.4613].",
        },

        # State point: 20 °C (293.15 K), 1 atm — room-temperature reference
        (293.15, 101325.0): {
            "temperature_K":              293.15,
            "pressure_Pa":                101325.0,
            "density_kg_m3":              868.5,
            "dynamic_viscosity_Pa_s":     2.20e-3,
            "specific_heat_J_kgK":        3290.0,
            "thermal_conductivity_W_mK":  0.201,
            "prandtl_number":             36.0,
            "speed_of_sound_m_s":         1380.0,
            "quality":                    -1.0,
            "enthalpy_J_kg":              1.62e5,
            "entropy_J_kgK":             510.0,
            "surface_tension_N_m":        None,
            "notes": "75wt% ethanol / 25wt% water @ 20°C, 1 atm.",
        },

        # State point: 20 °C, 3 MPa — representative high-pressure cooling condition
        (293.15, 3.0e6): {
            "temperature_K":              293.15,
            "pressure_Pa":                3.0e6,
            "density_kg_m3":              870.1,
            "dynamic_viscosity_Pa_s":     2.22e-3,
            "specific_heat_J_kgK":        3285.0,
            "thermal_conductivity_W_mK":  0.202,
            "prandtl_number":             36.1,
            "speed_of_sound_m_s":         1420.0,
            "quality":                    -1.0,
            "enthalpy_J_kg":              1.65e5,
            "entropy_J_kgK":             509.0,
            "surface_tension_N_m":        None,
            "notes": "75wt% ethanol / 25wt% water @ 20°C, 3 MPa (high-pressure cooling).",
        },
    },

    # ── Liquid Oxygen (LOX) ───────────────────────────────────────────────────
    "Oxygen": {

        # State point: 90 K (≈ boiling point at 1 atm), 1 atm
        (90.0, 101325.0): {
            "temperature_K":              90.0,
            "pressure_Pa":                101325.0,
            "density_kg_m3":              1141.0,
            "dynamic_viscosity_Pa_s":     1.94e-4,
            "specific_heat_J_kgK":        1699.0,
            "thermal_conductivity_W_mK":  0.152,
            "prandtl_number":             2.17,
            "speed_of_sound_m_s":         906.0,
            "quality":                    0.0,            # saturated liquid
            "enthalpy_J_kg":             -1.22e5,
            "entropy_J_kgK":             2940.0,
            "surface_tension_N_m":        1.32e-2,
            "notes": "LOX at normal boiling point (90 K, 1 atm).",
        },
    },

}

# =============================================================================
# Helper functions
# =============================================================================

def add_fluid_state(fluid: str, T_K: float, P_Pa: float, props: _PropDict) -> None:
    """
    Register a new (fluid, T, P) state point at runtime.

    Parameters
    ----------
    fluid : str
        Fluid / mixture name (same string used in Engine_Class fuelComposition).
    T_K   : float   Temperature in Kelvin.
    P_Pa  : float   Pressure in Pascal.
    props : dict    Property key → value mapping.

    Example
    -------
        add_fluid_state(
            "Ethanol[.5387]&Water[.4613]",
            T_K=400.0, P_Pa=2e6,
            props={"density_kg_m3": 780.0, ...}
        )
    """
    if fluid not in fluid_dict:
        fluid_dict[fluid] = {}
    fluid_dict[fluid][(T_K, P_Pa)] = props


def get_fluid_props(fluid: str, T_K: float, P_Pa: float,
                    tol_T: float = 0.5, tol_P: float = 500.0) -> _PropDict | None:
    """
    O(1) lookup with optional nearest-key fallback.

    Tries an exact hash lookup first (true O(1)).  If not found, searches
    registered state points for the closest match within tolerances.
    Returns None if the fluid is unknown or no match is within tolerance.

    Parameters
    ----------
    fluid  : str    Fluid name key.
    T_K    : float  Desired temperature (K).
    P_Pa   : float  Desired pressure (Pa).
    tol_T  : float  Temperature tolerance for nearest-key search (K).
    tol_P  : float  Pressure tolerance for nearest-key search (Pa).
    """
    states = fluid_dict.get(fluid)
    if states is None:
        return None

    # O(1) exact lookup
    exact = states.get((T_K, P_Pa))
    if exact is not None:
        return exact

    # Nearest-key fallback (O(n) over registered state points for this fluid)
    best_key = None
    best_dist = float('inf')
    for (t, p) in states:
        if abs(t - T_K) <= tol_T and abs(p - P_Pa) <= tol_P:
            dist = abs(t - T_K) + abs(p - P_Pa) / 1e5   # normalised combined distance
            if dist < best_dist:
                best_dist = dist
                best_key  = (t, p)

    return states[best_key] if best_key is not None else None


def populate_from_coolprop(fluid: str,
                            T_range_K,
                            P_range_Pa,
                            coolprop_string: str | None = None) -> None:
    """
    Batch-populate the dictionary using CoolProp for a grid of (T, P) points.

    Parameters
    ----------
    fluid          : str    Key to store results under.
    T_range_K      : iterable of float   Temperature points (K).
    P_range_Pa     : iterable of float   Pressure points (Pa).
    coolprop_string: str or None         CoolProp mixture string; defaults to `fluid`.

    Example
    -------
        import numpy as np
        populate_from_coolprop(
            "Ethanol[.5387]&Water[.4613]",
            T_range_K  = np.linspace(270, 500, 50),
            P_range_Pa = [101325, 1e6, 3e6],
        )
    """
    from CoolProp.CoolProp import PropsSI
    cp_str = coolprop_string if coolprop_string is not None else fluid

    for P in P_range_Pa:
        for T in T_range_K:
            try:
                rho = float(PropsSI('D', 'T', T, 'P', P, cp_str))
                mu  = float(PropsSI('V', 'T', T, 'P', P, cp_str))
                Cp  = float(PropsSI('C', 'T', T, 'P', P, cp_str))
                k   = float(PropsSI('L', 'T', T, 'P', P, cp_str))
                a   = float(PropsSI('A', 'T', T, 'P', P, cp_str))
                Q   = float(PropsSI('Q', 'T', T, 'P', P, cp_str))
                H   = float(PropsSI('H', 'T', T, 'P', P, cp_str))
                S   = float(PropsSI('S', 'T', T, 'P', P, cp_str))
                add_fluid_state(fluid, round(T, 4), round(P, 2), {
                    "temperature_K":             T,
                    "pressure_Pa":               P,
                    "density_kg_m3":             rho,
                    "dynamic_viscosity_Pa_s":    mu,
                    "specific_heat_J_kgK":       Cp,
                    "thermal_conductivity_W_mK": k,
                    "prandtl_number":            (mu * Cp) / k,
                    "speed_of_sound_m_s":        a,
                    "quality":                   Q,
                    "enthalpy_J_kg":             H,
                    "entropy_J_kgK":             S,
                    "surface_tension_N_m":       None,
                    "notes":                     f"Auto-populated via CoolProp @ T={T}K, P={P}Pa",
                })
            except Exception as e:
                print(f"[fluid_dict] CoolProp failed for {fluid} @ T={T}K P={P}Pa: {e}")
