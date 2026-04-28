# =============================================================================
# fluid_names.py
# Central fluid name translation table.
#
# PROBLEM
# -------
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
# SOLUTION
# --------
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
        "coolprop":     "NitrousOxide",
        "cea":          "N2O",
        "rocketprops":  None,           # rocketprops does not carry N2O
        "fluid_dict":   "NitrousOxide",
    },
    "H2O2": {
        "coolprop":     "H2O2",
        "cea":          "H2O2(L)",
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

    # ── Oxidisers (halogenated / hypergolic) ──────────────────────────────────

    "Fluorine": {
        # Liquid fluorine (LF2). Extremely reactive; used in high-Isp oxidiser
        # research and FLOX mixtures.  CoolProp and rocketprops both carry it.
        "coolprop":     "Fluorine",
        "cea":          "F2(L)",
        "rocketprops":  "LF2",
        "fluid_dict":   "Fluorine",
    },
    "ClF3": {
        # Chlorine trifluoride. Hypergolic with almost everything including
        # asbestos and concrete.  Neither CoolProp nor rocketprops carry it;
        # fluid properties must come from fluid_dict or literature.
        "coolprop":     None,
        "cea":          "ClF3(L)",
        "rocketprops":  None,
        "fluid_dict":   "ClF3",
    },
    "ClF5": {
        # Chlorine pentafluoride.  Higher Isp than ClF3; equally corrosive.
        # No CoolProp or rocketprops entry.
        "coolprop":     None,
        "cea":          "ClF5(L)",
        "rocketprops":  None,
        "fluid_dict":   "ClF5",
    },
    "HNO3": {
        # Nitric acid (concentrated / fuming).  Used in IRFNA blends and
        # early hypergolic systems.  No CoolProp pure-fluid entry.
        # For RFNA/WFNA/IRFNA mixtures use the IRFNA entry or a custom card.
        "coolprop":     None,
        "cea":          "HNO3(L)",
        "rocketprops":  None,
        "fluid_dict":   "HNO3",
    },

    # ── Fuels (hydrocarbons and miscellaneous) ────────────────────────────────

    "Gasoline": {
        # Automotive gasoline has no standard CEA thermochemical entry because
        # its composition varies.  Users who need CEA analysis should define a
        # custom propellant card approximating their blend, or use RP-1 as the
        # accepted engineering surrogate.  CoolProp and rocketprops do not
        # carry gasoline either.
        "coolprop":     None,
        "cea":          None,           # no standard CEA name — define custom card
        "rocketprops":  None,
        "fluid_dict":   "Gasoline",
    },
    "A50": {
        # Aerozine-50: 50 wt% hydrazine + 50 wt% UDMH.  CEA has a built-in
        # "A-50" entry in its thermodynamic library.  rocketprops carries it.
        # CoolProp does not have a mixture model for this blend.
        "coolprop":     None,
        "cea":          "A-50",
        "rocketprops":  "A50",
        "fluid_dict":   "A50",
    },
    "CarbonMonoxide": {
        # Liquid carbon monoxide.  Exotic cryogenic fuel; low density.
        "coolprop":     "CarbonMonoxide",
        "cea":          "CO(L)",
        "rocketprops":  None,
        "fluid_dict":   "CarbonMonoxide",
    },
    "Ethane": {
        # Liquid ethane.  Proposed as a storable methane-surrogate for some
        # Mars ISRU concepts.
        "coolprop":     "Ethane",
        "cea":          "C2H6(L)",
        "rocketprops":  None,
        "fluid_dict":   "Ethane",
    },
    "Butane": {
        # n-Butane.  CoolProp uses "n-Butane"; CEA uses C4H10(L) (n-butane
        # implied).  Aliases cover "n-Butane" and "isobutane" variants.
        "coolprop":     "n-Butane",
        "cea":          "C4H10(L)",
        "rocketprops":  None,
        "fluid_dict":   "Butane",
    },
    "Ethylene": {
        # Liquid ethylene (ethene, C2H4).  Used in some hypergolic research.
        "coolprop":     "Ethylene",
        "cea":          "C2H4(L)",
        "rocketprops":  None,
        "fluid_dict":   "Ethylene",
    },
    "Propylene": {
        # Liquid propylene (propene, C3H6).
        "coolprop":     "Propylene",
        "cea":          "C3H6(L)",
        "rocketprops":  None,
        "fluid_dict":   "Propylene",
    },
    "Benzene": {
        # Liquid benzene (C6H6).  High energy density aromatic fuel.
        "coolprop":     "Benzene",
        "cea":          "C6H6(L)",
        "rocketprops":  None,
        "fluid_dict":   "Benzene",
    },
    "Toluene": {
        # Liquid toluene (methylbenzene, C7H8).
        "coolprop":     "Toluene",
        "cea":          "C7H8(L)",
        "rocketprops":  None,
        "fluid_dict":   "Toluene",
    },
    "Xylene": {
        # Liquid xylene (dimethylbenzene, C8H10).  CoolProp routes "Xylene"
        # to o-Xylene; CEA uses o-C8H10(L) as the default isomer.  For other
        # isomers use "m-Xylene" or "p-Xylene" directly in CoolProp, and
        # m-C8H10(L) or p-C8H10(L) in CEA.
        "coolprop":     "Xylene",
        "cea":          "o-C8H10(L)",
        "rocketprops":  None,
        "fluid_dict":   "Xylene",
    },

    # ── Mixtures (CoolProp mixture string = canonical key) ────────────────────
    # Mixtures use the CoolProp bracket notation as the canonical name.
    # CEA handles mixtures via the card system — not a single string.
    # RocketProps does not support mixtures.
    # fluid_dict key = the same CoolProp mixture string.

    # ── Existing 75 wt% ethanol blend ────────────────────────────────────────
    # Mole fracs: Ethanol 0.5387, Water 0.4613
    "Ethanol[.5387]&Water[.4613]": {
        "coolprop":     "Ethanol[.5387]&Water[.4613]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Ethanol[.5387]&Water[.4613]",
    },

    # ── Ethanol / water blends ────────────────────────────────────────────────
    # Mole fractions derived from weight percents using MM(ethanol)=46.068,
    # MM(water)=18.015.  CoolProp bracket notation strips leading zero.
    # CEA entry is None — use newFuel card with C2H5OH(L) + H2O(L) wt%.

    # 90 wt% ethanol / 10 wt% water
    # x_eth = (0.90/46.068) / (0.90/46.068 + 0.10/18.015) = 0.7787
    "Ethanol[.7787]&Water[.2213]": {
        "coolprop":     "Ethanol[.7787]&Water[.2213]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Ethanol[.7787]&Water[.2213]",
    },

    # 80 wt% ethanol / 20 wt% water
    # x_eth = (0.80/46.068) / (0.80/46.068 + 0.20/18.015) = 0.6100
    "Ethanol[.61]&Water[.39]": {
        "coolprop":     "Ethanol[.61]&Water[.39]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Ethanol[.61]&Water[.39]",
    },

    # 70 wt% ethanol / 30 wt% water
    # x_eth = (0.70/46.068) / (0.70/46.068 + 0.30/18.015) = 0.4771
    "Ethanol[.4771]&Water[.5229]": {
        "coolprop":     "Ethanol[.4771]&Water[.5229]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Ethanol[.4771]&Water[.5229]",
    },

    # ── Methanol / water blends ───────────────────────────────────────────────
    # MM(methanol) = 32.042

    # 90 wt% methanol / 10 wt% water
    # x_meth = (0.90/32.042) / (0.90/32.042 + 0.10/18.015) = 0.8350
    "Methanol[.835]&Water[.165]": {
        "coolprop":     "Methanol[.835]&Water[.165]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Methanol[.835]&Water[.165]",
    },

    # 80 wt% methanol / 20 wt% water
    # x_meth = (0.80/32.042) / (0.80/32.042 + 0.20/18.015) = 0.6922
    "Methanol[.6922]&Water[.3078]": {
        "coolprop":     "Methanol[.6922]&Water[.3078]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Methanol[.6922]&Water[.3078]",
    },

    # 70 wt% methanol / 30 wt% water
    # x_meth = (0.70/32.042) / (0.70/32.042 + 0.30/18.015) = 0.5675
    "Methanol[.5675]&Water[.4325]": {
        "coolprop":     "Methanol[.5675]&Water[.4325]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Methanol[.5675]&Water[.4325]",
    },

    # ── IPA (isopropanol) / water blends ─────────────────────────────────────
    # CoolProp name for IPA: "Propanol".  MM(IPA) = 60.096.

    # 90 wt% IPA / 10 wt% water
    # x_ipa = (0.90/60.096) / (0.90/60.096 + 0.10/18.015) = 0.7296
    "Propanol[.7296]&Water[.2704]": {
        "coolprop":     "Propanol[.7296]&Water[.2704]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Propanol[.7296]&Water[.2704]",
    },

    # 80 wt% IPA / 20 wt% water
    # x_ipa = (0.80/60.096) / (0.80/60.096 + 0.20/18.015) = 0.5453
    "Propanol[.5453]&Water[.4547]": {
        "coolprop":     "Propanol[.5453]&Water[.4547]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Propanol[.5453]&Water[.4547]",
    },

    # 70 wt% IPA / 30 wt% water
    # x_ipa = (0.70/60.096) / (0.70/60.096 + 0.30/18.015) = 0.4116
    "Propanol[.4116]&Water[.5884]": {
        "coolprop":     "Propanol[.4116]&Water[.5884]",
        "cea":          None,
        "rocketprops":  None,
        "fluid_dict":   "Propanol[.4116]&Water[.5884]",
    },

    # ── FLOX mixtures (liquid fluorine + LOX) ────────────────────────────────
    # FLOX is a high-performance oxidiser blending liquid F2 and LOX.
    # CEA handles it via two oxidiser component cards, not a single name.
    # CoolProp can compute mixture properties using its mixture API.
    #
    # Mole fractions from weight percents using MM(F2)=38.00, MM(O2)=32.00.
    #
    # FLOX-70: 70 wt% F2 / 30 wt% LOX
    # x_F2 = (0.70/38.00) / (0.70/38.00 + 0.30/32.00) = 0.6637
    "Fluorine[.6637]&Oxygen[.3363]": {
        "coolprop":     "Fluorine[.6637]&Oxygen[.3363]",
        "cea":          None,   # use newOx card: F2(L) wt%=70 + O2 wt%=30
        "rocketprops":  None,
        "fluid_dict":   "Fluorine[.6637]&Oxygen[.3363]",
    },

    # FLOX-82.5: 82.5 wt% F2 / 17.5 wt% LOX  (highest Isp FLOX variant)
    # x_F2 = (0.825/38.00) / (0.825/38.00 + 0.175/32.00) = 0.7980
    "Fluorine[.798]&Oxygen[.202]": {
        "coolprop":     "Fluorine[.798]&Oxygen[.202]",
        "cea":          None,   # use newOx card: F2(L) wt%=82.5 + O2 wt%=17.5
        "rocketprops":  None,
        "fluid_dict":   "Fluorine[.798]&Oxygen[.202]",
    },

    # ── Alternative / legacy aliases (map to canonical entries above) ─────────
    # These aliases allow users to use informal names in TOML files without
    # breaking anything.  The alias resolver runs before the main table lookup.

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
    import sys, os, importlib, importlib.util

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
            from rocketprops.rocket_prop import get_prop
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
            from CoolProp.CoolProp import PropsSI
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
