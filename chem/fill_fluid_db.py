#!/usr/bin/env python3
# =============================================================================
# fill_fluid_db.py
# Populate fluid_props.db with thermodynamic properties from CoolProp and
# rocketprops.
#
# RUN THIS ONCE before using Engine_Class. After it completes, all property
# lookups in get_props() will be served directly from the database with zero
# CoolProp or rocketprops calls until a new (fluid, T, P) point is requested.
#
# USAGE
# -----
#   # Fill the default set of fluids for Engine_Class:
#   python fill_fluid_db.py
#
#   # Fill a specific fluid only:
#   python fill_fluid_db.py --fluid Oxygen
#   python fill_fluid_db.py --fluid "Ethanol[.5387]&Water[.4613]"
#
#   # Fill with a custom T/P grid:
#   python fill_fluid_db.py --fluid Ethanol --T-min 270 --T-max 500 --T-step 1 --P 2068430
#
#   # Show what is already in the database:
#   python fill_fluid_db.py --status
#
#   # Delete all rows for a fluid and re-fill from scratch:
#   python fill_fluid_db.py --fluid Oxygen --refill
#
# HOW IT WORKS
# ------------
# For each (fluid, T, P) point:
#   1. Check if the row already exists in the DB — skip if so (unless --refill).
#   2. Try rocketprops first (faster, pre-computed fits).
#   3. Fall back to CoolProp if rocketprops doesn't carry the fluid.
#   4. INSERT the result into fluid_properties with source="coolprop"|"rocketprops".
#
# The script prints progress and a summary at the end.
# Interrupted runs are safe — already-inserted rows are skipped on re-run.
#
# DEFAULT FILL PLAN
# -----------------
# Covers the fluids and conditions most commonly used in Engine_Class:
#
#   Fluid                         T range (K)    P (Pa)        Notes
#   ──────────────────────────────────────────────────────────────────────
#   Oxygen                        70 – 150       1 atm         LOX tank
#   Oxygen                        70 – 150       <Pc>          regen cooling
#   Ethanol                       250 – 500      1 atm         fuel tank
#   Ethanol                       250 – 500      <Pc>          regen cooling
#   Ethanol[.5387]&Water[.4613]   250 – 500      1 atm         75 wt% blend
#   Ethanol[.5387]&Water[.4613]   250 – 500      <Pc>
#   Methane                       100 – 200      1 atm         LCH4 tank
#   H2O2                          270 – 420      1 atm
#   N2O4                          260 – 350      1 atm
#   N2H4                          270 – 400      1 atm
#   ParaHydrogen                  15 – 35        1 atm         LH2 tank
#   RP1                           250 – 450      1 atm
#
# <Pc> is 2 068 000 Pa (~300 psi) by default; pass --Pc to change it.
# =============================================================================

import argparse
import sys
import os
import time

# Make sure fluid_db and fluid_names are importable from the same directory
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from fluid_db import get_db, _PROP_COLUMNS


# ---------------------------------------------------------------------------
# Conversion constants for rocketprops (imperial → SI)
# ---------------------------------------------------------------------------
BTU_lbmR_to_J_kgK  = 4186.8
BTU_hrftR_to_W_mK  = 1.730735
cP_to_Pa_s         = 1e-3
SG_to_kg_m3        = 1000.0


# ---------------------------------------------------------------------------
# Default fill plan
# Each entry: (canonical_name, T_min, T_max, T_step, [P_Pa, ...])
# ---------------------------------------------------------------------------

def build_default_plan(Pc_Pa: float) -> list:
    import numpy as np
    atm = 101325.0

    plan = [
        # Pure oxidisers / cryogens
        ("Oxygen",       70,  155,  1.0, [atm, Pc_Pa]),
        ("ParaHydrogen", 14,   34,  0.5, [atm]),
        ("N2O4",        260,  350,  1.0, [atm]),
        ("H2O2",        270,  420,  1.0, [atm, Pc_Pa]),
        ("N2H4",        270,  400,  1.0, [atm, Pc_Pa]),
        ("Fluorine",     54,   90,  1.0, [atm]),

        # Pure fuels
        ("Ethanol",     240,  510,  1.0, [atm, Pc_Pa]),
        ("Methanol",    180,  510,  1.0, [atm, Pc_Pa]),
        ("Methane",      95,  200,  1.0, [atm]),
        ("RP1",         230,  450,  1.0, [atm, Pc_Pa]),
        ("Propane",     180,  370,  1.0, [atm]),
        ("Ammonia",     200,  405,  1.0, [atm, Pc_Pa]),
        ("Hydrazine",   275,  400,  1.0, [atm, Pc_Pa]),
        ("MMH",         220,  360,  1.0, [atm, Pc_Pa]),

        # Ethanol / water blends
        ("Ethanol[.7787]&Water[.2213]", 240, 510, 1.0, [atm, Pc_Pa]),  # 90 wt%
        ("Ethanol[.61]&Water[.39]",     240, 510, 1.0, [atm, Pc_Pa]),  # 80 wt%
        ("Ethanol[.5387]&Water[.4613]", 240, 510, 1.0, [atm, Pc_Pa]),  # 75 wt%
        ("Ethanol[.4771]&Water[.5229]", 240, 510, 1.0, [atm, Pc_Pa]),  # 70 wt%

        # Methanol / water blends
        ("Methanol[.835]&Water[.165]",    180, 510, 1.0, [atm, Pc_Pa]),
        ("Methanol[.6922]&Water[.3078]",  180, 510, 1.0, [atm, Pc_Pa]),
        ("Methanol[.5675]&Water[.4325]",  180, 510, 1.0, [atm, Pc_Pa]),

        # IPA / water blends
        ("Propanol[.7296]&Water[.2704]",  230, 510, 1.0, [atm, Pc_Pa]),
        ("Propanol[.5453]&Water[.4547]",  230, 510, 1.0, [atm, Pc_Pa]),
        ("Propanol[.4116]&Water[.5884]",  230, 510, 1.0, [atm, Pc_Pa]),
    ]
    return plan


# ---------------------------------------------------------------------------
# Property fetchers
# ---------------------------------------------------------------------------

def fetch_from_rocketprops(rp_name: str, T_K: float, P_Pa: float) -> tuple[dict | None, str | None]:
    """
    Fetch all available properties from rocketprops.
    Returns (props_dict, None) on success, (None, error_string) on failure.
    """
    try:
        from rocketprops.rocket_prop import get_prop
        prop = get_prop(rp_name)
        P_atm = P_Pa / 101325.0
        rho  = prop.SG_compressed(T_K, P_atm) * SG_to_kg_m3
        mu   = prop.viscosity(T_K)             * cP_to_Pa_s
        Cp   = prop.specific_heat(T_K)         * BTU_lbmR_to_J_kgK
        k    = prop.thermal_conductivity(T_K)  * BTU_hrftR_to_W_mK
        Pr   = (mu * Cp) / k if k else None
        return {
            "density_kg_m3":             rho,
            "dynamic_viscosity_Pa_s":    mu,
            "specific_heat_J_kgK":       Cp,
            "thermal_conductivity_W_mK": k,
            "prandtl_number":            Pr,
            "speed_of_sound_m_s":        None,
            "quality":                   None,
            "enthalpy_J_kg":             None,
            "entropy_J_kgK":             None,
        }, None
    except Exception as e:
        return None, f"rocketprops: {type(e).__name__}: {e}" 


def fetch_from_coolprop(cp_name: str, T_K: float, P_Pa: float) -> tuple[dict | None, str | None]:
    """
    Fetch all available properties from CoolProp PropsSI.
    Returns (props_dict, None) on success, (None, error_string) on failure.
    """
    try:
        from CoolProp.CoolProp import PropsSI
        rho = float(PropsSI('D', 'T', T_K, 'P', P_Pa, cp_name))
        mu  = float(PropsSI('V', 'T', T_K, 'P', P_Pa, cp_name))
        Cp  = float(PropsSI('C', 'T', T_K, 'P', P_Pa, cp_name))
        k   = float(PropsSI('L', 'T', T_K, 'P', P_Pa, cp_name))
        try:
            a = float(PropsSI('A', 'T', T_K, 'P', P_Pa, cp_name))
        except Exception:
            a = None
        try:
            Q = float(PropsSI('Q', 'T', T_K, 'P', P_Pa, cp_name))
        except Exception:
            Q = None
        try:
            H = float(PropsSI('H', 'T', T_K, 'P', P_Pa, cp_name))
        except Exception:
            H = None
        try:
            S = float(PropsSI('S', 'T', T_K, 'P', P_Pa, cp_name))
        except Exception:
            S = None
        Pr = (mu * Cp) / k if k else None
        return {
            "density_kg_m3":             rho,
            "dynamic_viscosity_Pa_s":    mu,
            "specific_heat_J_kgK":       Cp,
            "thermal_conductivity_W_mK": k,
            "prandtl_number":            Pr,
            "speed_of_sound_m_s":        a,
            "quality":                   Q,
            "enthalpy_J_kg":             H,
            "entropy_J_kgK":             S,
        }, None
    except Exception as e:
        return None, f"coolprop: {type(e).__name__}: {e}" 


# ---------------------------------------------------------------------------
# Core fill function
# ---------------------------------------------------------------------------

def fill_fluid(canonical: str, T_min: float, T_max: float, T_step: float,
               P_list: list[float], refill: bool = False,
               verbose: bool = True) -> dict:
    """
    Fill the DB for one fluid across a temperature grid and list of pressures.

    Returns a summary dict: {"inserted": int, "skipped": int, "failed": int}
    """
    import numpy as np

    # Resolve names for each backend
    try:
        from fluid_names import translate, canonical as canon_fn
        canon = canon_fn(canonical)
        try:
            rp_name = translate(canon, 'rocketprops')
        except ValueError:
            rp_name = None
        try:
            cp_name = translate(canon, 'coolprop')
        except ValueError:
            cp_name = None
    except ImportError:
        # fluid_names not available — use canonical as-is for CoolProp
        canon   = canonical
        rp_name = None
        cp_name = canonical

    db = get_db()

    if refill:
        deleted = db.delete_fluid(canon)
        if verbose and deleted:
            print(f"  Deleted {deleted} existing rows for '{canon}'")

    T_grid = np.round(np.arange(T_min, T_max + T_step * 0.5, T_step), 4)

    inserted = skipped = failed = 0
    total = len(T_grid) * len(P_list)

    # Track unique failure reasons across all (T, P) points.
    # Key: error string  →  Value: list of (T_K, P_Pa) that produced it.
    # Deduplicated so the same CoolProp or rocketprops error message doesn't
    # print hundreds of times for a fluid that is entirely unsupported.
    failure_log: dict[str, list[tuple]] = {}

    for P_Pa in P_list:
        P_r = round(P_Pa, 2)
        batch = []

        for T_K in T_grid:
            T_r = round(float(T_K), 4)

            # Skip if already in DB
            if not refill and db.lookup(canon, T_r, P_r) is not None:
                skipped += 1
                continue

            # Try rocketprops first
            props = source = None
            rp_err = cp_err = None

            if rp_name:
                props, rp_err = fetch_from_rocketprops(rp_name, T_r, P_r)
                if props:
                    source = "rocketprops"

            # Fall back to CoolProp on rocketprops miss
            if props is None and cp_name:
                props, cp_err = fetch_from_coolprop(cp_name, T_r, P_r)
                if props:
                    source = "coolprop"

            if props is None:
                failed += 1
                # Build a combined error key from whichever sources were tried
                parts = []
                if rp_err:
                    parts.append(rp_err)
                if cp_err:
                    parts.append(cp_err)
                if not parts:
                    parts.append("no source available (fluid not in rocketprops or CoolProp)")
                err_key = " | ".join(parts)
                failure_log.setdefault(err_key, []).append((T_r, P_r))
                continue

            props['temperature_K'] = T_r
            props['pressure_Pa']   = P_r
            batch.append((canon, T_r, P_r, props, source,
                          f"Filled by fill_fluid_db.py via {source}"))
            inserted += 1

        # Bulk insert this pressure slice and commit once
        if batch:
            db.insert_many(batch, commit=True)

        if verbose:
            status = f"  P={P_r/1e3:.1f} kPa : {inserted} new, {skipped} skipped"
            if failed:
                status += f", {failed} FAILED"
            print(status)

    # Print deduplicated failure report
    if failure_log and verbose:
        print(f"\n  ── Failure report for '{canon}' ──")
        for err_msg, points in failure_log.items():
            # Show first 3 failing points so the log doesn't flood the terminal
            sample = ", ".join(f"T={t}K P={p/1e3:.0f}kPa" for t, p in points[:3])
            extra  = f" ... and {len(points)-3} more" if len(points) > 3 else ""
            print(f"  [{len(points)} points] {err_msg}")
            print(f"    e.g. {sample}{extra}")
        print()

    return {
        "inserted":    inserted,
        "skipped":     skipped,
        "failed":      failed,
        "failure_log": failure_log,   # caller can inspect programmatically
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Populate fluid_props.db with thermodynamic properties."
    )
    parser.add_argument(
        "--fluid", default=None,
        help="Canonical fluid name to fill (default: all fluids in the plan)."
    )
    parser.add_argument(
        "--T-min", type=float, default=None,
        help="Minimum temperature (K). Overrides per-fluid default."
    )
    parser.add_argument(
        "--T-max", type=float, default=None,
        help="Maximum temperature (K). Overrides per-fluid default."
    )
    parser.add_argument(
        "--T-step", type=float, default=1.0,
        help="Temperature step size (K). Default: 1.0."
    )
    parser.add_argument(
        "--P", type=float, default=None,
        help="Single pressure (Pa) to fill. Overrides per-fluid default."
    )
    parser.add_argument(
        "--Pc", type=float, default=2068000.0,
        help="Chamber pressure (Pa) included in the default fill plan. "
             "Default: 2068000 Pa (~300 psi)."
    )
    parser.add_argument(
        "--refill", action="store_true",
        help="Delete existing rows and re-fill from scratch."
    )
    parser.add_argument(
        "--status", action="store_true",
        help="Print DB contents summary and exit."
    )
    args = parser.parse_args()

    db = get_db()

    # ── Status report ────────────────────────────────────────────────────────
    if args.status:
        print(f"\nDatabase: {db._path}")
        print(f"Total rows: {db.count()}\n")
        fluids = db.list_fluids()
        if not fluids:
            print("  (empty)")
        else:
            print(f"{'Fluid':<45} {'Rows':>6}  T range at 1 atm")
            print("-" * 75)
            for f in fluids:
                n   = db.count(f)
                rng = db.get_T_range(f, 101325.0)
                rng_str = f"{rng[0]:.1f} – {rng[1]:.1f} K" if rng else "(no 1-atm data)"
                print(f"  {f:<43} {n:>6}  {rng_str}")
        print()
        return

    # ── Build fill plan ──────────────────────────────────────────────────────
    plan = build_default_plan(args.Pc)

    if args.fluid:
        # Filter to only the requested fluid, or create a single-entry plan
        matched = [row for row in plan if row[0] == args.fluid]
        if not matched:
            # User specified a fluid not in the default plan — add it
            T_min  = args.T_min  or 250.0
            T_max  = args.T_max  or 500.0
            T_step = args.T_step or 1.0
            P_list = [args.P] if args.P else [101325.0]
            matched = [(args.fluid, T_min, T_max, T_step, P_list)]
        plan = matched

    # Apply T/P overrides if provided
    if args.T_min or args.T_max or args.P:
        new_plan = []
        for (name, T_min, T_max, T_step, P_list) in plan:
            T_min  = args.T_min  or T_min
            T_max  = args.T_max  or T_max
            T_step = args.T_step or T_step
            P_list = [args.P]    if args.P else P_list
            new_plan.append((name, T_min, T_max, T_step, P_list))
        plan = new_plan

    # ── Execute fill ─────────────────────────────────────────────────────────
    total_inserted = total_skipped = total_failed = 0
    all_failures: dict = {}   # fluid_name → failure_log from fill_fluid()
    t0 = time.time()

    for (name, T_min, T_max, T_step, P_list) in plan:
        print(f"\nFilling '{name}'  T={T_min}–{T_max} K  step={T_step} K  "
              f"pressures={[f'{p/1e3:.1f}kPa' for p in P_list]}")
        summary = fill_fluid(name, T_min, T_max, T_step, P_list,
                             refill=args.refill, verbose=True)
        total_inserted += summary["inserted"]
        total_skipped  += summary["skipped"]
        total_failed   += summary["failed"]
        if summary["failure_log"]:
            all_failures[name] = summary["failure_log"]

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"Done in {elapsed:.1f}s")
    print(f"  Inserted : {total_inserted}")
    print(f"  Skipped  : {total_skipped}  (already in DB)")
    print(f"  Failed   : {total_failed}")
    print(f"  DB size  : {os.path.getsize(db._path)/1024:.1f} KB")
    print(f"  DB path  : {db._path}")

    # Cross-fluid failure summary — useful when running the full default plan
    if all_failures:
        print(f"\n  ── Fluids with failures ──")
        for fluid_name, log in all_failures.items():
            total_pts = sum(len(v) for v in log.values())
            print(f"  {fluid_name}: {total_pts} failed points")
            for err_msg in log:
                print(f"    {err_msg}")
        print()
        print("  These fluids are not available in rocketprops or CoolProp.")
        print("  Add data manually via fluid_db.get_db().insert() or add")
        print("  entries to fluid_dict.py for use as a manual fallback.")


if __name__ == "__main__":
    main()
