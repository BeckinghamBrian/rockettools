# =============================================================================
# fluid_names.py
# Central fluid name translation table.
#
# Three fluid property sources each use different name strings for the same
# fluid, and none of them are compatible with each other out of the box:
#
#   Source        LOX example    Ethanol example          Mixture example
#   ----------    -----------    --------------------     ----------------------------
#   NASA CEA      "O2"           "C2H5OH(L)"              not applicable (card system)
#   CoolProp      "Oxygen"       "Ethanol"                "Ethanol[.75]&Water[.25]"
#   RocketProps   "LOX"          "Ethanol"                not supported
#   fluid_dict    "Oxygen"       "Ethanol"                "Ethanol[.75]&Water[.25]"
#
# The (L) suffix used in CEA names (e.g. "C2H5OH(L)") breaks CoolProp and
# RocketProps lookups.  The CoolProp mixture string format
# "Fluid1[x1]&Fluid2[x2]" is not valid for CEA or RocketProps.
# "O2" works for CEA but not CoolProp ("Oxygen") or RocketProps ("LOX").
#
# This module defines a single authoritative table (FLUID_NAME_MAP) that maps
# a human-readable "canonical" name to the correct string for each system.
# Engine_Class calls translate(name, target) to get the right string before
# passing it to any property source.
#
# CANONICAL NAME
# --------------
# The canonical name is the human-readable common name used in TOML input files
# and in fluid_dict keys.  It is chosen to match the CoolProp / fluid_dict
# string wherever possible (since that system is most widely used in the code).
# For pure fluids this is just the CoolProp name.
# For mixtures it is the CoolProp mixture string.
#
# ADDING NEW FLUIDS
# -----------------
# Add a row to FLUID_NAME_MAP.  Set the value to None for any system that
# does not support that fluid.  translate() will raise a clear error if a
# None entry is requested.
#
# USAGE
# -----
#   from fluid_names import translate, coolprop_name, cea_name, rocketprops_name
#
#   translate("Oxygen",   "coolprop")      # → "Oxygen"
#   translate("Oxygen",   "cea")           # → "O2"
#   translate("Oxygen",   "rocketprops")   # → "LOX"
#   translate("Ethanol",  "cea")           # → "C2H5OH(L)"
#   translate("H2O2",     "rocketprops")   # → "H2O2"
#
#   # Mixture: CoolProp string is the canonical key
#   translate("Ethanol[.5387]&Water[.4613]", "coolprop")   # → the same string
#   translate("Ethanol[.5387]&Water[.4613]", "rocketprops") # → None (not supported)
#
# =============================================================================

from __future__ import annotations
from typing import Literal, Optional
from CoolProp.CoolProp import PropsSI
from rocketprops.rocket_prop import get_prop
import sys, os, importlib, importlib.util

# Target system identifiers
Target = Literal["coolprop", "cea", "rocketprops", "fluid_dict"]

# =============================================================================
# Main translation table
# Format: "canonical_name": {"coolprop": str|None, "cea": str|None,
#                             "rocketprops": str|None, "fluid_dict": str|None}
#
# fluid_dict key = coolprop key for pure fluids.
# For mixtures the fluid_dict key IS the coolprop mixture string.
# "None" means the system does not support that fluid.
# =============================================================================

FLUID_NAME_MAP: dict[str, dict[str, Optional[str]]] = {

    # ── Oxidisers ─────────────────────────────────────────────────────────────

    "Oxygen": {
        "coolprop":     "Oxygen",
        "cea":          "O2",
        "rocketprops":  "LOX",
        "fluid_dict":   "Oxygen",
    },
    "N2O4": {
        "coolprop":     "N2O4",
        "cea":          "N2O4(L)",
        "rocketprops":  "N2O4",
        "fluid_dict":   "N2O4",
    },
    "NitrousOxide": {
        "coolprop":     None,
        "cea":          "N2O",
        "rocketprops":  "N2O",           # rocketprops does not carry N2O
        "fluid_dict":   "NitrousOxide",
    },
    "H2O2": {
        "coolprop":     "H2O2",
        "cea":          "H2O2",
        "rocketprops":  "H2O2",
        "fluid_dict":   "H2O2",
    },
    "IRFNA": {
        "coolprop":     None,           # CoolProp does not carry IRFNA
        "cea":          "IRFNA-III",
        "rocketprops":  "IRFNA",
        "fluid_dict":   "IRFNA",
    },
    "ParaHydrogen": {
        "coolprop":     "ParaHydrogen",
        "cea":          "H2(L)",
        "rocketprops":  "LH2",
        "fluid_dict":   "ParaHydrogen",
    },

    # ── Fuels ─────────────────────────────────────────────────────────────────

    "Ethanol": {
        "coolprop":     "Ethanol",
        "cea":          "C2H5OH(L)",
        "rocketprops":  "Ethanol",
        "fluid_dict":   "Ethanol",
    },
    "Methanol": {
        "coolprop":     "Methanol",
        "cea":          "CH3OH(L)",
        "rocketprops":  "Methanol",
        "fluid_dict":   "Methanol",
    },
    "Propanol": {
        "coolprop":     "Propanol",
        "cea":          "C3H7OH(L)",
        "rocketprops":  "IPA",          # rocketprops uses IPA for isopropanol
        "fluid_dict":   "Propanol",
    },
    "RP1": {
        "coolprop":     None,           # CoolProp does not carry RP-1
        "cea":          "RP-1",
        "rocketprops":  "RP1",
        "fluid_dict":   "RP1",
    },
    "Methane": {
        "coolprop":     "Methane",
        "cea":          "CH4(L)",
        "rocketprops":  "LCH4",
        "fluid_dict":   "Methane",
    },
    "Hydrazine": {
        "coolprop":     "Hydrazine",
        "cea":          "N2H4(L)",
        "rocketprops":  "N2H4",
        "fluid_dict":   "Hydrazine",
    },
    "MMH": {
        "coolprop":     None,           # CoolProp does not carry MMH
        "cea":          "CH6N2(L)",
        "rocketprops":  "MMH",
        "fluid_dict":   "MMH",
    },
    "UDMH": {
        "coolprop":     None,
        "cea":          "UDMH",
        "rocketprops":  "UDMH",
        "fluid_dict":   "UDMH",
    },
    "Ammonia": {
        "coolprop":     "Ammonia",
        "cea":          "NH3(L)",
        "rocketprops":  "NH3",
        "fluid_dict":   "Ammonia",
    },

    # ── Common pure fluids ────────────────────────────────────────────────────

    "Water": {
        "coolprop":     "Water",
        "cea":          "H2O(L)",
        "rocketprops":  None,
        "fluid_dict":   "Water",
    },
    "Nitrogen": {
        "coolprop":     "Nitrogen",
        "cea":          "N2(L)",
        "rocketprops":  "LN2",
        "fluid_dict":   "Nitrogen",
    },
    "Propane": {
        "coolprop":     "Propane",
        "cea":          "C3H8(L)",
        "rocketprops":  "Propane",
        "fluid_dict":   "Propane",
    },

    # ── Mixtures (CoolProp mixture string = canonical key) ────────────────────
    # Mixtures use the CoolProp bracket notation as the canonical name.
    # CEA handles mixtures via the card system — not a single string.
    # RocketProps does not support mixtures.
    # fluid_dict key = the same CoolProp mixture string.

    "Ethanol[.5387]&Water[.4613]": {
        "coolprop":     "Ethanol[.5387]&Water[.4613]",
        "cea":          None,           # mixtures use CEA card, not a name
        "rocketprops":  None,
        "fluid_dict":   "Ethanol[.5387]&Water[.4613]",
    },

    # ── Alternative / legacy aliases (map to canonical entries above) ─────────
    # These aliases allow users to use informal names in TOML files without
    # breaking anything.  The alias resolver runs before the main table lookup.

}

# Aliases: informal / legacy name → canonical name in FLUID_NAME_MAP
# All lookups resolve aliases before hitting the main table.
ALIASES: dict[str, str] = {
    # CEA-style names that users might type
    "O2":           "Oxygen",
    "O2(L)":        "Oxygen",
    "H2(L)":        "ParaHydrogen",
    "LH2":          "ParaHydrogen",
    "N2O4(L)":      "N2O4",
    "H2O2(L)":      "H2O2",
    "N2H4(L)":      "Hydrazine",
    "CH6N2(L)":     "MMH",
    "CH3OH(L)":     "Methanol",
    "C2H5OH(L)":    "Ethanol",
    "C3H7OH(L)":    "Propanol",
    "CH4(L)":       "Methane",
    "N2(L)":        "Nitrogen",
    "NH3(L)":       "Ammonia",
    "H2O(L)":       "Water",
    "C3H8(L)":      "Propane",
    "N2O":          "NitrousOxide",
    "RP-1":         "RP1",
    # RocketProps-style names
    "LOX":          "Oxygen",
    "LCH4":         "Methane",
    "LN2":          "Nitrogen",
    "IPA":          "Propanol",
    # Common shorthands
    "lox":          "Oxygen",
    "ethanol":      "Ethanol",
    "methanol":     "Methanol",
    "water":        "Water",
    "h2o2":         "H2O2",
    "n2o4":         "N2O4",
}


# =============================================================================
# Public API
# =============================================================================

def canonical(name: str) -> str:
    """
    Resolve any alias or variant name to its canonical form.
    Returns the name unchanged if it is already canonical or unknown.
    """
    return ALIASES.get(name, name)


def translate(name: str, target: Target) -> str:
    """
    Translate a fluid name to the string required by the target system.

    Parameters
    ----------
    name   : str     Fluid name — canonical, alias, or CoolProp mixture string.
    target : str     One of "coolprop", "cea", "rocketprops", "fluid_dict".

    Returns
    -------
    str  The name string appropriate for the target system.

    Raises
    ------
    KeyError   If the fluid is not in the translation table.
    ValueError If the target system does not support this fluid (mapped to None).
    """
    canon = canonical(name)

    # Pass through CoolProp mixture strings unchanged — they are already
    # in the correct format for CoolProp and fluid_dict, and cannot be
    # used with CEA or RocketProps.
    if '&' in canon:
        if target in ('coolprop', 'fluid_dict'):
            return canon
        raise ValueError(
            f"Fluid mixture '{name}' cannot be used with '{target}'. "
            "Mixture property lookup is only supported via CoolProp and fluid_dict."
        )

    entry = FLUID_NAME_MAP.get(canon)
    if entry is None:
        # Unknown fluid — return unchanged and let the caller handle the error.
        # This is intentional: it allows users to pass non-standard CoolProp
        # strings (e.g. custom backends) without breaking the translation layer.
        return name

    result = entry.get(target)
    if result is None:
        raise ValueError(
            f"Fluid '{name}' (canonical: '{canon}') is not supported by '{target}'. "
            f"Available targets for this fluid: "
            f"{[t for t, v in entry.items() if v is not None]}"
        )
    return result


def coolprop_name(name: str) -> str:
    """Shorthand: translate name to CoolProp format."""
    return translate(name, "coolprop")

def cea_name(name: str) -> str:
    """Shorthand: translate name to NASA CEA format."""
    return translate(name, "cea")

def rocketprops_name(name: str) -> str:
    """Shorthand: translate name to RocketProps format."""
    return translate(name, "rocketprops")

def fluid_dict_key(name: str) -> str:
    """Shorthand: translate name to fluid_dict key format (same as CoolProp)."""
    return translate(name, "fluid_dict")


def supported_targets(name: str) -> list[str]:
    """Return the list of target systems that support this fluid."""
    canon = canonical(name)
    entry = FLUID_NAME_MAP.get(canon, {})
    return [t for t, v in entry.items() if v is not None]


def list_fluids() -> list[str]:
    """Return all canonical fluid names in the translation table."""
    return sorted(FLUID_NAME_MAP.keys())


# =============================================================================
# Property source priority and dispatcher
# =============================================================================
# get_fluid_props() tries sources in this priority order:
#   1. fluid_dict  — O(1) in-memory lookup, zero I/O if cached
#   2. rocketprops — fast, covers most rocket propellants, no CoolProp dep
#   3. CoolProp    — broadest coverage, supports mixtures, slowest first call
#
# This function is the single entry point Engine_Class uses to get a fluid
# property value.  It handles all name translation internally.

def get_props(name: str, T_K: float, P_Pa: float,
              properties: list[str],
              prefer_source: str | None = None) -> dict:
    """
    Get fluid thermodynamic properties from the best available source.

    Tries sources in priority order: fluid_dict → rocketprops → CoolProp.
    Name translation is handled internally — pass any supported name format.

    Parameters
    ----------
    name        : str          Fluid name (any format).
    T_K         : float        Temperature (K).
    P_Pa        : float        Pressure (Pa).
    properties  : list[str]    Property keys to retrieve.
                               Recognised keys (SI units):
                                 "density_kg_m3"
                                 "dynamic_viscosity_Pa_s"
                                 "specific_heat_J_kgK"
                                 "thermal_conductivity_W_mK"
                                 "prandtl_number"
    prefer_source : str|None   Force a specific source: "fluid_dict",
                               "rocketprops", or "coolprop". None = auto.

    Returns
    -------
    dict  {property_key: float, "_source": str}
          "_source" records which system provided the values.
    """

    canon = canonical(name)

    # ── 1. fluid_dict ─────────────────────────────────────────────────────────
    if prefer_source in (None, 'fluid_dict'):
        try:
            _fd_key = fluid_dict_key(canon)
            _fd_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    'fluid_dict.py')
            if 'fluid_dict' not in sys.modules and os.path.exists(_fd_path):
                _spec = importlib.util.spec_from_file_location('fluid_dict', _fd_path)
                _mod  = importlib.util.module_from_spec(_spec)
                _spec.loader.exec_module(_mod)
                sys.modules['fluid_dict'] = _mod
            if 'fluid_dict' in sys.modules:
                _fd = sys.modules['fluid_dict'].fluid_dict
                _entry = _fd.get(_fd_key, {}).get((round(T_K, 4), round(P_Pa, 2)))
                if _entry is not None:
                    result = {k: _entry[k] for k in properties if k in _entry}
                    if len(result) == len(properties):
                        result['_source'] = 'fluid_dict'
                        return result
        except Exception:
            pass

    # ── 2. rocketprops ────────────────────────────────────────────────────────
    if prefer_source in (None, 'rocketprops'):
        try:
            rp_key = translate(canon, 'rocketprops')
            prop = get_prop(rp_key)

            # unit conversion constants
            BTU_lbmR_to_J_kgK   = 4186.8
            BTU_hrftR_to_W_mK   = 1.730735
            cP_to_Pa_s          = 1e-3
            SG_to_kg_m3         = 1000.0      # × water density at 4°C

            result = {}
            T_R = T_K * 1.8    # Rankine for rocketprops internal calls

            if 'density_kg_m3' in properties:
                P_atm = P_Pa / 101325.0
                result['density_kg_m3'] = prop.SG_compressed(T_K, P_atm) * SG_to_kg_m3

            if 'dynamic_viscosity_Pa_s' in properties:
                result['dynamic_viscosity_Pa_s'] = prop.viscosity(T_K) * cP_to_Pa_s

            if 'specific_heat_J_kgK' in properties:
                result['specific_heat_J_kgK'] = prop.specific_heat(T_K) * BTU_lbmR_to_J_kgK

            if 'thermal_conductivity_W_mK' in properties:
                result['thermal_conductivity_W_mK'] = (
                    prop.thermal_conductivity(T_K) * BTU_hrftR_to_W_mK)

            if 'prandtl_number' in properties:
                mu = result.get('dynamic_viscosity_Pa_s',
                                prop.viscosity(T_K) * cP_to_Pa_s)
                Cp = result.get('specific_heat_J_kgK',
                                prop.specific_heat(T_K) * BTU_lbmR_to_J_kgK)
                k  = result.get('thermal_conductivity_W_mK',
                                prop.thermal_conductivity(T_K) * BTU_hrftR_to_W_mK)
                result['prandtl_number'] = (mu * Cp) / k

            if len(result) == len(properties):
                result['_source'] = 'rocketprops'
                return result
        except Exception:
            pass

    # ── 3. CoolProp ───────────────────────────────────────────────────────────
    if prefer_source in (None, 'coolprop'):
        try:
            cp_key = coolprop_name(canon)
            result = {}

            CP_MAP = {
                'density_kg_m3':             'D',
                'dynamic_viscosity_Pa_s':    'V',
                'specific_heat_J_kgK':       'C',
                'thermal_conductivity_W_mK': 'L',
            }
            for prop_key in properties:
                if prop_key == 'prandtl_number':
                    continue
                cp_prop = CP_MAP.get(prop_key)
                if cp_prop:
                    result[prop_key] = float(PropsSI(cp_prop, 'T', T_K, 'P', P_Pa, cp_key))

            if 'prandtl_number' in properties:
                mu = result.get('dynamic_viscosity_Pa_s',
                                float(PropsSI('V', 'T', T_K, 'P', P_Pa, cp_key)))
                Cp = result.get('specific_heat_J_kgK',
                                float(PropsSI('C', 'T', T_K, 'P', P_Pa, cp_key)))
                k  = result.get('thermal_conductivity_W_mK',
                                float(PropsSI('L', 'T', T_K, 'P', P_Pa, cp_key)))
                result['prandtl_number'] = (mu * Cp) / k

            result['_source'] = 'coolprop'
            return result
        except Exception as e:
            raise RuntimeError(
                f"All property sources failed for fluid '{name}' at T={T_K}K, P={P_Pa}Pa.\n"
                f"Last error (CoolProp): {e}"
            )

    raise RuntimeError(
        f"No property source available for fluid '{name}' at T={T_K}K, P={P_Pa}Pa "
        f"with prefer_source='{prefer_source}'."
    )
