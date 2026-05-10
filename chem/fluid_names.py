"""
fluid_names.py
Central fluid name translation table.

Three fluid property sources each use different name strings for the same
fluid, and none of them are compatible with each other out of the box:

  Source        LOX example    Ethanol example          Mixture example
  ----------    -----------    --------------------     ----------------------------
  NASA CEA      "O2"           "C2H5OH(L)"              not applicable (card system)
  CoolProp      "Oxygen"       "Ethanol"                "Ethanol[.75]&Water[.25]"
  RocketProps   "LOX"          "Ethanol"                not supported
  fluid_dict    "Oxygen"       "Ethanol"                "Ethanol[.75]&Water[.25]"

The (L) suffix used in CEA names (e.g. "C2H5OH(L)") breaks CoolProp and
RocketProps lookups.  The CoolProp mixture string format
"Fluid1[x1]&Fluid2[x2]" is not valid for CEA or RocketProps.
"O2" works for CEA but not CoolProp ("Oxygen") or RocketProps ("LOX").

This module defines a single authoritative table (FLUID_NAME_MAP) that maps
a human-readable "canonical" name to the correct string for each system.
Engine_Class calls translate(name, target) to get the right string before
passing it to any property source.

CANONICAL NAME
--------------
The canonical name is the human-readable common name used in TOML input files
and in fluid_dict keys.  It is chosen to match the CoolProp / fluid_dict
string wherever possible (since that system is most widely used in the code).
For pure fluids this is just the CoolProp name.
For mixtures it is the CoolProp mixture string.

ADDING NEW FLUIDS
-----------------
Add a row to FLUID_NAME_MAP.  Set the value to None for any system that
does not support that fluid.  translate() will raise a clear error if a
None entry is requested.

USAGE
-----
  from fluid_names import translate, coolprop_name, cea_name, rocketprops_name

  translate("Oxygen",   "coolprop")      # → "Oxygen"
  translate("Oxygen",   "cea")           # → "O2"
  translate("Oxygen",   "rocketprops")   # → "LOX"
  translate("Ethanol",  "cea")           # → "C2H5OH(L)"
  translate("H2O2",     "rocketprops")   # → "H2O2"

  # Mixture: CoolProp string is the canonical key
  translate("Ethanol[.5387]&Water[.4613]", "coolprop")   # → the same string
  translate("Ethanol[.5387]&Water[.4613]", "rocketprops") # → None (not supported)
"""
from __future__ import annotations
from typing import Literal, Optional

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
    "NTO": {
        "coolprop":     None,           # CoolProp does not carry N2O4
        "cea":          "N2O4(L)",
        "rocketprops":  "N2O4",
        "fluid_dict":   "N2O4",
    },
    "MON3": {
        "coolprop":     None,           # CoolProp does not carry MON3
        "cea":          "MON3",
        "rocketprops":  None,           # RocketProps does not carry MON3
        "fluid_dict":   "MON3",
    },
    "MON15": {
        "coolprop":     None,           # CoolProp does not carry MON15
        "cea":          "MON15",
        "rocketprops":  None,           # RocketProps does not carry MON15
        "fluid_dict":   "MON15",
    },
    "MON25": {
        "coolprop":     None,           # CoolProp does not carry MON25
        "cea":          "MON25",
        "rocketprops":  None,           # RocketProps does not carry MON25
        "fluid_dict":   "MON25",
    },
    "N2F4": {
        "coolprop":     None,           # CoolProp does not carry N2F4
        "cea":          "N2F4",
        "rocketprops":  None,           # RocketProps does not carry N2F4
        "fluid_dict":   "N2F4",
    },
    "NitrousOxide": {
        "coolprop":     "NitrousOxide",
        "cea":          "N2O",
        "rocketprops":  "N20",           
        "fluid_dict":   "NitrousOxide",
    },
    "Peroxide": {
        "coolprop":     None,           # CoolProp does not carry peroxide
        "cea":          "H2O2(L)",
        "rocketprops":  "H2O2",
        "fluid_dict":   "H2O2",
    },
    "RFNA": {
        "coolprop":     None,           # CoolProp does not carry IRFNA
        "cea":          "HNO3(L)",
        "rocketprops":  "IRFNA",
        "fluid_dict":   "IRFNA",
    },
    "WFNA": {
        "coolprop":     None,           # CoolProp does not carry WFNA
        "cea":          "HNO3(L)",
        "rocketprops":  None,           # RocketProps does not carry WFNA
        "fluid_dict":   "HNO3",
    },
    "Fluorine": {
        "coolprop":     "Fluorine",
        "cea":          "F2(L)",
        "rocketprops":  "F2",
        "fluid_dict":   "Fluorine",
    },
    "ClF3": {
        "coolprop":     None,           # CoolProp does not carry CLF3
        "cea":          "ClF3(L)",
        "rocketprops":  None,           # RocketProps does not carry CLF3
        "fluid_dict":   "ClF3",
    },
    "CLF5": {
        "coolprop":     None,           # CoolProp does not carry CLF5
        "cea":          "ClF5",
        "rocketprops":  "CLF5",
        "fluid_dict":   "CLF5",
    },
    "FLOX70": {
        "coolprop":     "Fluorine[.6637]&Oxygen[.3363]",
        "cea":          "FLOX70",
        "rocketprops":  None,           # RocketProps does not carry FLOX
        "fluid_dict":   "Fluorine[.6637]&Oxygen[.3363]",
    },
    "FLOX82": {
        "coolprop":     "Fluorine[.798]&Oxygen[.202]",
        "cea":          None,   # use newOx card: F2(L) wt%=82.5 + O2 wt%=17.5
        "rocketprops":  None,           # RocketProps does not carry FLOX
        "fluid_dict":   "Fluorine[.798]&Oxygen[.202]",
    },
    "FLOX60": {
        "coolprop":     "Fluorine[.559]&Oxygen[.441]",
        "cea":          "FLOX60",
        "rocketprops":  None,           # RocketProps does not carry FLOX
        "fluid_dict":   "Fluorine[.559]&Oxygen[.441]",
    },
    "FLOX80": {
        "coolprop":     "Fluorine[.771]&Oxygen[.229]",
        "cea":          "FLOX80",
        "rocketprops":  None,           # RocketProps does not carry FLOX
        "fluid_dict":   "Fluorine[.771]&Oxygen[.229]",
    },
    "CarbonMonoxide": {
        "coolprop":     "CarbonMonoxide",
        "cea":          "CO",
        "rocketprops":  None,
        "fluid_dict":   "CarbonMonoxide",
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
    "IPA": {
        "coolprop":     None,           # CoolProp does not carry isopropanol or n-propanol
        "cea":          "C3H8O,1propanol",
        "rocketprops":  None,          # RocketProps does not carry isopropanol or n-propanol
        "fluid_dict":   "Propanol",
    },
    "RP1": {
        "coolprop":     None,           # CoolProp does not carry RP-1
        "cea":          "RP-1",
        "rocketprops":  "RP1",
        "fluid_dict":   "RP1",
    },
    "Hydrogen": {
        "coolprop":     "ParaHydrogen",
        "cea":          "H2(L)",
        "rocketprops":  "PH2",
        "fluid_dict":   "Hydrogen",
    },
    "Methane": {
        "coolprop":     "Methane",
        "cea":          "CH4(L)",
        "rocketprops":  "Methane",
        "fluid_dict":   "Methane",
    },
    "Hydrazine": {
        "coolprop":     None,           # CoolProp does not carry hydrazine
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
        "coolprop":     None,           # CoolProp does not carry UDMH
        "cea":          "C2H8N2(L),UDMH	",
        "rocketprops":  "UDMH",
        "fluid_dict":   "UDMH",
    },
    "Ammonia": {
        "coolprop":     "Ammonia",
        "cea":          "NH3(L)",
        "rocketprops":  "NH3",
        "fluid_dict":   "Ammonia",
    },
    "Propane": {
        "coolprop":     "Propane",
        "cea":          "C3H8(L)",
        "rocketprops":  "Propane",
        "fluid_dict":   "Propane",
    },
    "Gasoline": {
        "coolprop":     None,
        "cea":          "C8H18(L),n-octa",
        "rocketprops":  None,
        "fluid_dict":   "Gasoline",
    },
    "JetA": {
        "coolprop":     None,
        "cea":          "Jet-A(L)",
        "rocketprops":  None,
        "fluid_dict":   "JetA",
    },
    "M20": {
        "coolprop":     None,
        "cea":          "M20",
        "rocketprops":  None,
        "fluid_dict":   "M20",
    },
    "A50": {
        "coolprop":     None,
        "cea":          "A50",
        "rocketprops":  "A50",
        "fluid_dict":   "A50",
    },
    "Ethane": {
        "coolprop":     "Ethane",
        "cea":          "C2H6(L)",
        "rocketprops":  None,
        "fluid_dict":   "Ethane",
    },
    "Butane": {
        "coolprop":     "n-Butane",
        "cea":          "C4H10(L)",
        "rocketprops":  None,
        "fluid_dict":   "Butane",
    },
    "Ethylene": {
        "coolprop":     "Ethylene",
        "cea":          "C2H4(L)",
        "rocketprops":  None,
        "fluid_dict":   "Ethylene",
    },
    "Propylene": {
        "coolprop":     "Propylene",
        "cea":          "C3H6(L)",
        "rocketprops":  None,
        "fluid_dict":   "Propylene",
    },
    "Benzene": {
        "coolprop":     "Benzene",
        "cea":          "C6H6(L)",
        "rocketprops":  None,
        "fluid_dict":   "Benzene",
    },
    "Toluene": {
        "coolprop":     "Toluene",
        "cea":          "C7H8(L)",
        "rocketprops":  None,
        "fluid_dict":   "Toluene",
    },
    "Xylene": {
        "coolprop":     "o-Xylene",
        "cea":          "o-C8H10(L)",
        "rocketprops":  None,
        "fluid_dict":   "Xylene",
    },
    "Ethanol75": {
        "coolprop":     "Ethanol[.5387]&Water[.4613]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Ethanol[.5387]&Water[.4613]",
    },
    "Ethanol90": {
        "coolprop":     "Ethanol[.7787]&Water[.2213]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Ethanol[.7787]&Water[.2213]",
    },
    "Ethanol80": {
        "coolprop":     "Ethanol[.61]&Water[.39]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Ethanol[.61]&Water[.39]",
    },
    "Ethanol70": {
        "coolprop":     "Ethanol[.4771]&Water[.5229]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Ethanol[.4771]&Water[.5229]",
    },
    "Methanol90": {
        "coolprop":     "Methanol[.835]&Water[.165]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Methanol[.835]&Water[.165]",
    },
    "Methanol80": {
        "coolprop":     "Methanol[.6922]&Water[.3078]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Methanol[.6922]&Water[.3078]",
    },
    "Methanol70": {
        "coolprop":     "Methanol[.5675]&Water[.4325]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Methanol[.5675]&Water[.4325]",
    },
    "IPA90": {
        "coolprop":     "Propanol[.7296]&Water[.2704]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Propanol[.7296]&Water[.2704]",
    },
    "IPA80": {
        "coolprop":     "Propanol[.5453]&Water[.4547]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Propanol[.5453]&Water[.4547]",
    },
    "IPA70": {
        "coolprop":     "Propanol[.4116]&Water[.5884]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Propanol[.4116]&Water[.5884]",
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
}

# Aliases: informal / legacy name → canonical name in FLUID_NAME_MAP
# All lookups resolve aliases before hitting the main table.
ALIASES: dict[str, str] = {

    # ── CEA liquid-phase suffixes → canonical ─────────────────────────────────
    # CEA appends (L) to denote the liquid phase.  Users writing CEA names in
    # their TOML files will have them resolved to the canonical name first.
    "O2":           "Oxygen",
    "O2(L)":        "Oxygen",
    "F2(L)":        "Fluorine",
    "H2(L)":        "ParaHydrogen",
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
    "CO(L)":        "CarbonMonoxide",
    "C2H6(L)":      "Ethane",
    "C4H10(L)":     "Butane",
    "C2H4(L)":      "Ethylene",
    "C3H6(L)":      "Propylene",
    "C6H6(L)":      "Benzene",
    "C7H8(L)":      "Toluene",
    "o-C8H10(L)":   "Xylene",
    "ClF3(L)":      "ClF3",
    "ClF5(L)":      "ClF5",
    "HNO3(L)":      "HNO3",
    "N2O":          "NitrousOxide",
    "RP-1":         "RP1",
    "A-50":         "A50",

    # ── RocketProps short names → canonical ───────────────────────────────────
    "LOX":          "Oxygen",
    "LF2":          "Fluorine",
    "LCH4":         "Methane",
    "LN2":          "Nitrogen",
    "LH2":          "ParaHydrogen",
    "IPA":          "Propanol",
    "n-Butane":     "Butane",

    # ── Mixture weight-percent shorthand aliases ───────────────────────────────
    # These let users type a readable name in TOML instead of the full
    # CoolProp bracket string.  The canonical key IS the CoolProp string so
    # the alias maps directly to it.
    "Ethanol90":                        "Ethanol[.7787]&Water[.2213]",
    "Ethanol80":                        "Ethanol[.61]&Water[.39]",
    "Ethanol75":                        "Ethanol[.5387]&Water[.4613]",
    "Ethanol70":                        "Ethanol[.4771]&Water[.5229]",
    "Methanol90":                       "Methanol[.835]&Water[.165]",
    "Methanol80":                       "Methanol[.6922]&Water[.3078]",
    "Methanol70":                       "Methanol[.5675]&Water[.4325]",
    "IPA90":                            "Propanol[.7296]&Water[.2704]",
    "IPA80":                            "Propanol[.5453]&Water[.4547]",
    "IPA70":                            "Propanol[.4116]&Water[.5884]",
    "FLOX70":                           "Fluorine[.6637]&Oxygen[.3363]",
    "FLOX82":                           "Fluorine[.798]&Oxygen[.202]",

    # ── Common informal / lowercase shorthands ────────────────────────────────
    "lox":              "Oxygen",
    "ethanol":          "Ethanol",
    "methanol":         "Methanol",
    "water":            "Water",
    "h2o2":             "H2O2",
    "n2o4":             "N2O4",
    "fluorine":         "Fluorine",
    "clf3":             "ClF3",
    "clf5":             "ClF5",
    "ipa":              "Propanol",
    "a50":              "A50",
    "aerozine":         "A50",
    "aerozine-50":      "A50",
    "gasoline":         "Gasoline",
    "petrol":           "Gasoline",
    "benzene":          "Benzene",
    "toluene":          "Toluene",
    "xylene":           "Xylene",
    "ethane":           "Ethane",
    "butane":           "Butane",
    "ethylene":         "Ethylene",
    "propylene":        "Propylene",
    "nitric acid":      "HNO3",
    "fuming nitric":    "HNO3",
    "carbon monoxide":  "CarbonMonoxide",
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
# Module-level fluid_dict cache
# =============================================================================
# Loaded once on first call to get_props(); reused for the lifetime of the
# Python session via sys.modules.  Engine_Class never imports fluid_dict or
# CoolProp or rocketprops directly — all calls go through get_props() here.

_fluid_dict_loaded: bool = False          # sentinel so we only attempt load once

def _load_fluid_dict() -> None:
    """Import fluid_dict.py into sys.modules if not already loaded."""
    global _fluid_dict_loaded
    if _fluid_dict_loaded:
        return
    _fluid_dict_loaded = True             # set before the attempt so we don't retry on failure
    import sys, os, importlib.util
    _path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fluid_dict.py')
    if os.path.exists(_path) and 'fluid_dict' not in sys.modules:
        try:
            _spec = importlib.util.spec_from_file_location('fluid_dict', _path)
            _mod  = importlib.util.module_from_spec(_spec)
            _spec.loader.exec_module(_mod)
            sys.modules['fluid_dict'] = _mod
        except Exception as _e:
            print(f"[fluid_names] WARNING: could not load fluid_dict.py — {_e}")


def _write_back(fd_key: str, T_K: float, P_Pa: float, result: dict) -> None:
    """
    Write a property set returned by rocketprops or CoolProp back into the
    in-memory fluid_dict so subsequent calls are served from cache.
    Only writes keys that fluid_dict recognises (excludes '_source').
    """
    import sys
    if 'fluid_dict' not in sys.modules:
        return
    try:
        _add = sys.modules['fluid_dict'].add_fluid_state
        _fd  = sys.modules['fluid_dict'].fluid_dict
        entry = {k: v for k, v in result.items() if k != '_source'}
        entry['temperature_K'] = T_K
        entry['pressure_Pa']   = P_Pa
        entry.setdefault('notes', f"Auto-populated via {result.get('_source','?')}")
        _add(fd_key, round(T_K, 4), round(P_Pa, 2), entry)
    except Exception:
        pass   # write-back is best-effort; never raise from here


# =============================================================================
# Property source priority and dispatcher
# =============================================================================
# get_props() is the single entry point for ALL fluid property lookups.
# Engine_Class, and any other code in this repository, must call this
# function instead of importing CoolProp, rocketprops, or fluid_dict directly.
#
# Source priority order (fastest / most targeted first):
#   1. fluid_dict  — in-memory hash lookup; zero I/O after first session call
#   2. rocketprops — pre-computed curve fits; fast; covers common rocket propellants
#   3. CoolProp    — equation-of-state accuracy; supports mixtures; slowest first call
#
# On a cache miss from rocketprops or CoolProp the result is written back into
# the in-memory fluid_dict so the next call in the same session is a cache hit.

def get_props(name: str, T_K: float, P_Pa: float,
              properties: list[str],
              prefer_source: Optional[str] = None) -> dict:
    """
    Return thermodynamic properties for a fluid from the best available source.

    This is the **only** function that should be called for fluid properties.
    It handles name translation, source selection, unit conversion, and
    write-back to the in-memory cache — the caller needs none of that logic.

    Parameters
    ----------
    name        : str           Fluid name in any supported format (canonical,
                                CEA, CoolProp, RocketProps, alias, or mixture
                                string).  Examples: "Oxygen", "O2", "LOX",
                                "Ethanol90", "Ethanol[.5387]&Water[.4613]".
    T_K         : float         Temperature (K).
    P_Pa        : float         Pressure (Pa).
    properties  : list[str]     Property keys to retrieve.  All values are SI.
                                Supported keys:
                                  "density_kg_m3"              kg/m³
                                  "dynamic_viscosity_Pa_s"     Pa·s
                                  "specific_heat_J_kgK"        J/(kg·K)
                                  "thermal_conductivity_W_mK"  W/(m·K)
                                  "prandtl_number"             dimensionless
    prefer_source : str | None  Pin a specific source: "fluid_dict",
                                "rocketprops", or "coolprop".
                                None (default) = automatic priority order.

    Returns
    -------
    dict
        {property_key: float, ..., "_source": str}
        "_source" records which backend provided the values.

    Raises
    ------
    RuntimeError
        If all sources fail or the requested properties cannot be satisfied.
    """
    import sys

    # Ensure fluid_dict is loaded into the module cache
    _load_fluid_dict()

    canon  = canonical(name)
    fd_key = fluid_dict_key(canon)         # canonical = CoolProp key = fluid_dict key

    # ── 1. fluid_dict ─────────────────────────────────────────────────────────
    # Pure O(1) hash lookups — no I/O, no external library calls.
    if prefer_source in (None, 'fluid_dict'):
        if 'fluid_dict' in sys.modules:
            _fd    = sys.modules['fluid_dict'].fluid_dict
            _entry = _fd.get(fd_key, {}).get((round(T_K, 4), round(P_Pa, 2)))
            if _entry is not None:
                result = {k: _entry[k] for k in properties if k in _entry}
                if len(result) == len(properties):
                    result['_source'] = 'fluid_dict'
                    return result
                # Partial hit — fall through to get remaining properties from
                # a live source, then write the full set back to cache.

    # ── 2. rocketprops ────────────────────────────────────────────────────────
    # Fast pre-computed fits.  Covers LOX, LH2, LCH4, Ethanol, MMH, N2H4,
    # N2O4, RP-1, H2O2, A50, and other common propellants.
    # Does NOT support CoolProp mixture strings.
    if prefer_source in (None, 'rocketprops') and '&' not in canon:
        try:
            rp_key = translate(canon, 'rocketprops')   # raises ValueError if unsupported
            from rocketprops.rocket_prop import get_prop
            prop = get_prop(rp_key)

            BTU_lbmR_to_J_kgK  = 4186.8
            BTU_hrftR_to_W_mK  = 1.730735
            cP_to_Pa_s         = 1e-3
            SG_to_kg_m3        = 1000.0

            result = {}

            if 'density_kg_m3' in properties:
                result['density_kg_m3'] = prop.SG_compressed(T_K, P_Pa / 101325.0) * SG_to_kg_m3
            if 'dynamic_viscosity_Pa_s' in properties:
                result['dynamic_viscosity_Pa_s'] = prop.viscosity(T_K) * cP_to_Pa_s
            if 'specific_heat_J_kgK' in properties:
                result['specific_heat_J_kgK'] = prop.specific_heat(T_K) * BTU_lbmR_to_J_kgK
            if 'thermal_conductivity_W_mK' in properties:
                result['thermal_conductivity_W_mK'] = (
                    prop.thermal_conductivity(T_K) * BTU_hrftR_to_W_mK)
            if 'prandtl_number' in properties:
                _mu = result.get('dynamic_viscosity_Pa_s',
                                 prop.viscosity(T_K) * cP_to_Pa_s)
                _Cp = result.get('specific_heat_J_kgK',
                                 prop.specific_heat(T_K) * BTU_lbmR_to_J_kgK)
                _k  = result.get('thermal_conductivity_W_mK',
                                 prop.thermal_conductivity(T_K) * BTU_hrftR_to_W_mK)
                result['prandtl_number'] = (_mu * _Cp) / _k

            if len(result) == len(properties):
                result['_source'] = 'rocketprops'
                _write_back(fd_key, T_K, P_Pa, result)
                return result
        except (ValueError, ImportError, Exception):
            pass   # not supported by rocketprops → fall through to CoolProp

    # ── 3. CoolProp ───────────────────────────────────────────────────────────
    # Equation-of-state accuracy.  Supports pure fluids and mixtures.
    if prefer_source in (None, 'coolprop'):
        try:
            from CoolProp.CoolProp import PropsSI
            cp_key = coolprop_name(canon)   # raises ValueError if unsupported

            CP_PROPS = {
                'density_kg_m3':             'D',
                'dynamic_viscosity_Pa_s':    'V',
                'specific_heat_J_kgK':       'C',
                'thermal_conductivity_W_mK': 'L',
            }
            result = {}
            for prop_key in properties:
                if prop_key == 'prandtl_number':
                    continue
                cp_id = CP_PROPS.get(prop_key)
                if cp_id:
                    result[prop_key] = float(PropsSI(cp_id, 'T', T_K, 'P', P_Pa, cp_key))

            if 'prandtl_number' in properties:
                _mu = result.get('dynamic_viscosity_Pa_s',
                                 float(PropsSI('V', 'T', T_K, 'P', P_Pa, cp_key)))
                _Cp = result.get('specific_heat_J_kgK',
                                 float(PropsSI('C', 'T', T_K, 'P', P_Pa, cp_key)))
                _k  = result.get('thermal_conductivity_W_mK',
                                 float(PropsSI('L', 'T', T_K, 'P', P_Pa, cp_key)))
                result['prandtl_number'] = (_mu * _Cp) / _k

            result['_source'] = 'coolprop'
            _write_back(fd_key, T_K, P_Pa, result)
            return result

        except Exception as _cp_err:
            raise RuntimeError(
                f"[fluid_names] All property sources failed for '{name}' "
                f"(canonical='{canon}') at T={T_K} K, P={P_Pa} Pa.\n"
                f"  fluid_dict: {'loaded' if 'fluid_dict' in sys.modules else 'not found'}\n"
                f"  rocketprops: not available for this fluid\n"
                f"  CoolProp error: {_cp_err}"
            ) from _cp_err

    raise RuntimeError(
        f"[fluid_names] No property source available for '{name}' "
        f"at T={T_K} K, P={P_Pa} Pa with prefer_source={prefer_source!r}."
    )


def get_density(name: str, T_K: float, P_Pa: float) -> float:
    """Convenience wrapper: return density (kg/m³) only."""
    return get_props(name, T_K, P_Pa, ['density_kg_m3'])['density_kg_m3']


def get_viscosity(name: str, T_K: float, P_Pa: float) -> float:
    """Convenience wrapper: return dynamic viscosity (Pa·s) only."""
    return get_props(name, T_K, P_Pa, ['dynamic_viscosity_Pa_s'])['dynamic_viscosity_Pa_s']


def get_specific_heat(name: str, T_K: float, P_Pa: float) -> float:
    """Convenience wrapper: return specific heat Cp (J/kg·K) only."""
    return get_props(name, T_K, P_Pa, ['specific_heat_J_kgK'])['specific_heat_J_kgK']


def get_thermal_conductivity(name: str, T_K: float, P_Pa: float) -> float:
    """Convenience wrapper: return thermal conductivity (W/m·K) only."""
    return get_props(name, T_K, P_Pa, ['thermal_conductivity_W_mK'])['thermal_conductivity_W_mK']


def get_prandtl(name: str, T_K: float, P_Pa: float) -> float:
    """Convenience wrapper: return Prandtl number only."""
    return get_props(name, T_K, P_Pa, ['prandtl_number'])['prandtl_number']
