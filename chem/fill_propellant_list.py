#!/usr/bin/env python3
# =============================================================================
# fill_propellant_list.py
# Populate fluid_props.db for the full rockettools propellant list.
#
# USAGE
# -----
#   python fill_propellant_list.py                  # fill everything
#   python fill_propellant_list.py --status         # show what is already in DB
#   python fill_propellant_list.py --dry-run        # print plan, insert nothing
#   python fill_propellant_list.py --fluid Ethanol  # single fluid only
#   python fill_propellant_list.py --refill         # delete and re-fill all
#   python fill_propellant_list.py --Pc 3000000     # use a different Pc (Pa)
#
# TEMPERATURE BOUNDS
# ------------------
# T_min is set just above the freeze point so only liquid-phase data is stored.
# T_max is set well above the normal boiling point — this is intentional.
# Rocket engine coolant is often pressurised above its saturation pressure,
# meaning the fluid remains liquid (or becomes supercritical) at temperatures
# that would cause it to boil at 1 atm. CoolProp and rocketprops both return
# valid properties above the normal boiling point as long as the pressure is
# above the saturation pressure at that temperature.
#
# Pressures stored per fluid:
#   1 atm (101 325 Pa)  — tank / ambient conditions
#   Pc (default 300 psi = 2 068 000 Pa) — representative chamber / cooling pressure
#
# FLUIDS THAT CANNOT BE AUTO-FILLED
# ----------------------------------
# Some fluids are not in CoolProp or rocketprops. fill_fluid_db.py will report
# them as failed points. Properties for these must be added manually via:
#   fluid_db.get_db().insert(canonical_name, T_K, P_Pa, props_dict)
#
# Fluids likely to fail (no CoolProp or rocketprops entry):
#   CLF3, CLF5, FLOX variants, Nitric Acid, MON3/15/25, N2F4, JetA, M20,
#   Gasoline — add manually from literature if needed.
# =============================================================================

import sys
import os
import time
import argparse
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from fill_fluid_db import fill_fluid
from fluid_db import get_db

# =============================================================================
# Full propellant list with temperature bounds
#
# Each entry:
#   (name, T_min_K, T_max_K, notes)
#
# T_min: 1 K above freeze point — keeps data in the liquid phase.
# T_max: set above the normal boiling point (bp) so the dataset covers
#        pressurised liquid and supercritical conditions relevant to
#        regenerative cooling. The value chosen is roughly:
#           pure cryogens:  T_crit + 20–50 K
#           storable fuels: T_bp  + 100–150 K (well below decomposition)
#           hypergolics:    T_bp  + 100 K     (below decomposition onset)
#           blends:         water bp + 110 K  (mixture stays below flash point)
#
# Normal boiling points are noted in comments so the above-bp range is clear.
# =============================================================================

PROPELLANT_LIST = [

    # ── Oxidisers ─────────────────────────────────────────────────────────────

    # Liquid oxygen: bp 90 K, supercrit 155 K.
    # Upper 275 K covers warm pressurised LOX and gaseous oxygen.
    ("Oxygen",          55,   275,  "LOX  bp=90K  supercrit=155K"),

    # NTO (pure N2O4): freeze 262 K, bp 294 K, decomp onset ~430 K.
    # Upper 425 K = bp + 131 K, well below decomp.
    ("NTO",            263,   425,  "N2O4  bp=294K  decomp>430K"),

    # MON variants (N2O4 + dissolved NO): NO depresses bp and freeze point.
    # MON3:  3% NO,  bp ~289 K.   MON15: 15% NO, bp ~278 K.   MON25: 25% NO, bp ~268 K.
    ("MON3",           259,   420,  "3%NO in NTO  bp~289K"),
    ("MON15",          249,   410,  "15%NO in NTO  bp~278K"),
    ("MON25",          239,   400,  "25%NO in NTO  bp~268K"),

    # Dinitrogen tetrafluoride (N2F4): bp 200 K, crit 309 K.
    # Exotic high-energy oxidiser; data from literature only.
    ("N2F4",           110,   305,  "N2F4  bp=200K  crit=309K"),

    # Nitrous oxide: freeze 182 K, bp 185 K, supercrit 310 K.
    # Upper = critical temperature.
    ("NitrousOxide",   183,   310,  "N2O  bp=185K  supercrit=310K"),

    # 90 wt% hydrogen peroxide: freeze 272 K, bp ~423 K, decomp ~570 K.
    # Upper 520 K = bp + 97 K, safe margin below decomp.
    ("Peroxide",       273,   520,  "90% H2O2  bp~423K  decomp>570K"),

    # Nitric acid (HNO3): freeze 231 K, bp 356 K, fumes heavily above ~400 K.
    # Upper 450 K = bp + 94 K.
    ("Nitric Acid",    232,   450,  "HNO3  bp=356K  fumes>400K"),

    # Liquid fluorine: freeze 54 K, bp 85 K, supercrit 144 K.
    # Upper 200 K covers supercritical use.
    ("Fluorine",        55,   200,  "LF2  bp=85K  supercrit=144K"),

    # Chlorine trifluoride: freeze 197 K, bp 285 K.
    # Upper 400 K = bp + 115 K.  CoolProp/RP likely fail — manual fill needed.
    ("CLF3",           198,   400,  "ClF3  bp=285K  manual fill likely needed"),

    # Chlorine pentafluoride: freeze 170 K, bp 260 K.
    # Upper 370 K = bp + 110 K.  Manual fill likely needed.
    ("CLF5",           171,   370,  "ClF5  bp=260K  manual fill likely needed"),

    # FLOX variants (liquid fluorine + LOX): bp ~86–88 K depending on ratio.
    # Upper 200 K covers supercritical region. CoolProp mixture — may fail.
    ("FLOX60",          63,   200,  "60%F2/40%LOX  bp~88K"),
    ("FLOX70",          61,   200,  "70%F2/30%LOX  bp~87K"),
    ("FLOX80",          59,   200,  "80%F2/20%LOX  bp~86K"),
    ("FLOX82",          58,   200,  "82.5%F2/17.5%LOX  bp~86K"),

    # ── Fuels ─────────────────────────────────────────────────────────────────

    # Carbon monoxide: freeze 68 K, bp 82 K, supercrit 133 K.
    ("CarbonMonoxide",  69,   200,  "LCO  bp=82K  supercrit=133K"),

    # Ethanol: freeze 159 K, bp 351 K, crit 514 K.
    # Upper 500 K = bp + 149 K, well below crit.
    ("Ethanol",        160,   500,  "EtOH  bp=351K  crit=514K"),

    # Methanol: freeze 175 K, bp 338 K, crit 513 K.
    ("Methanol",       176,   490,  "MeOH  bp=338K  crit=513K"),

    # IPA (isopropanol): freeze 185 K, bp 356 K, crit 508 K.
    ("IPA",            186,   490,  "IPA  bp=356K  crit=508K"),

    # RP-1: no sharp bp; distillation range ~450–510 K.
    # Upper 600 K covers film-cooling superheated conditions.
    ("RP1",            231,   600,  "RP-1  distil=450-510K  crit~675K"),

    # Liquid hydrogen: freeze 14 K, bp 20 K, supercrit 33 K.
    # Upper 80 K covers warm pressurised GH2.
    ("Hydrogen",        15,    80,  "LH2  bp=20K  supercrit=33K"),

    # Liquid methane: freeze 91 K, bp 112 K, supercrit 190 K.
    ("Methane",         92,   250,  "LCH4  bp=112K  supercrit=190K"),

    # Hydrazine: freeze 275 K, bp 387 K, decomp onset ~600 K.
    # Upper 500 K = bp + 113 K, safe below decomp.
    ("Hydrazine",      276,   500,  "N2H4  bp=387K  decomp>600K"),

    # MMH (monomethylhydrazine): freeze 221 K, bp 360 K, crit 567 K.
    ("MMH",            222,   480,  "MMH  bp=360K  crit=567K"),

    # UDMH: freeze 216 K, bp 336 K, crit 523 K.
    ("UDMH",           217,   450,  "UDMH  bp=336K  crit=523K"),

    # Ammonia: freeze 195 K, bp 240 K, crit 406 K.
    ("Ammonia",        196,   430,  "NH3  bp=240K  crit=406K"),

    # Propane: freeze 86 K, bp 231 K, crit 370 K.
    ("Propane",         87,   400,  "C3H8  bp=231K  crit=370K"),

    # Gasoline: no standard formula; approximated as C8H18.
    # CoolProp/RP will fail — manual fill from literature needed.
    ("Gasoline",       211,   450,  "~C8H18 surrogate  bp~330K  manual fill needed"),

    # Jet-A: distillation range ~430–570 K, freeze ~233 K.
    # No standard CoolProp entry; similar to RP-1.
    ("JetA",           234,   600,  "Jet-A  distil~430-570K  manual fill needed"),

    # M20 (20% methanol / 80% hydrazine blend): bp ~320 K.
    ("M20",            169,   430,  "20%MeOH/80%N2H4  bp~320K"),

    # Aerozine-50 (50% N2H4 / 50% UDMH): freeze 274 K, bp ~341 K.
    ("A50",            275,   480,  "A-50  bp~341K"),

    # Ethane: freeze 90 K, bp 185 K, supercrit 305 K.
    ("Ethane",          91,   310,  "C2H6  bp=185K  supercrit=305K"),

    # n-Butane: freeze 136 K, bp 273 K, crit 425 K.
    ("Butane",         137,   450,  "n-C4H10  bp=273K  crit=425K"),

    # Ethylene: freeze 104 K, bp 169 K, supercrit 282 K.
    ("Ethylene",       105,   300,  "C2H4  bp=169K  supercrit=282K"),

    # Propylene: freeze 88 K, bp 225 K, crit 365 K.
    ("Propylene",       89,   390,  "C3H6  bp=225K  crit=365K"),

    # Benzene: freeze 279 K, bp 353 K, crit 562 K.
    ("Benzene",        280,   580,  "C6H6  bp=353K  crit=562K"),

    # Toluene: freeze 178 K, bp 384 K, crit 592 K.
    ("Toluene",        179,   600,  "C7H8  bp=384K  crit=592K"),

    # Xylene (o-xylene default): freeze 248 K, bp 417 K, crit 630 K.
    ("Xylene",         249,   640,  "o-C8H10  bp=417K  crit=630K"),

    # ── Ethanol / water blends ────────────────────────────────────────────────
    # bp rises with increasing water content; upper bound same for all.
    # Lower bound set by ethanol freeze point (water raises it slightly but
    # the ethanol freeze point is the controlling constraint for the blend).

    # 90 wt% ethanol: bp ~352 K
    ("Ethanol90",      160,   490,  "90wt%EtOH/10wt%H2O  bp~352K"),
    # 80 wt% ethanol: bp ~354 K
    ("Ethanol80",      160,   490,  "80wt%EtOH/20wt%H2O  bp~354K"),
    # 75 wt% ethanol: bp ~356 K
    ("Ethanol75",      160,   490,  "75wt%EtOH/25wt%H2O  bp~356K"),
    # 70 wt% ethanol: bp ~360 K
    ("Ethanol70",      160,   490,  "70wt%EtOH/30wt%H2O  bp~360K"),

    # ── Methanol / water blends ───────────────────────────────────────────────

    # 90 wt% methanol: bp ~341 K
    ("Methanol90",     176,   490,  "90wt%MeOH/10wt%H2O  bp~341K"),
    # 80 wt% methanol: bp ~347 K
    ("Methanol80",     176,   490,  "80wt%MeOH/20wt%H2O  bp~347K"),
    # 70 wt% methanol: bp ~354 K
    ("Methanol70",     176,   490,  "70wt%MeOH/30wt%H2O  bp~354K"),

    # ── IPA / water blends ────────────────────────────────────────────────────

    # 90 wt% IPA: bp ~357 K
    ("IPA90",          186,   490,  "90wt%IPA/10wt%H2O  bp~357K"),
    # 80 wt% IPA: bp ~360 K
    ("IPA80",          186,   490,  "80wt%IPA/20wt%H2O  bp~360K"),
    # 70 wt% IPA: bp ~364 K
    ("IPA70",          186,   490,  "70wt%IPA/30wt%H2O  bp~364K"),

    # ── Common supporting fluids ──────────────────────────────────────────────

    # Water: freeze 273 K, bp 373 K, crit 647 K.
    # Upper 580 K covers high-pressure steam cooling and supercritical water.
    ("Water",          274,   580,  "H2O  bp=373K  crit=647K"),

    # Liquid nitrogen: freeze 63 K, bp 77 K, supercrit 126 K.
    # Used for pressurant cooldown and inert purge.
    ("Nitrogen",        64,   200,  "LN2  bp=77K  supercrit=126K"),
]


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Fill fluid_props.db for the full rockettools propellant list."
    )
    parser.add_argument("--fluid",   default=None,
                        help="Fill a single fluid by name instead of the full list.")
    parser.add_argument("--Pc",      type=float, default=2068000.0,
                        help="Chamber/cooling pressure in Pa. Default 2 068 000 Pa (~300 psi).")
    parser.add_argument("--T-step",  type=float, default=1.0,
                        help="Temperature grid step in K. Default 1.0 K.")
    parser.add_argument("--refill",  action="store_true",
                        help="Delete existing rows and re-fill from scratch.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the fill plan without inserting anything.")
    parser.add_argument("--status",  action="store_true",
                        help="Print DB contents summary and exit.")
    args = parser.parse_args()

    db   = get_db()
    P_atm = 101325.0
    Pc    = args.Pc
    P_range = np.arange(1e5, 10e6, 1000)  # Pressure from 1 bar to 100 bar

    # ── Status ────────────────────────────────────────────────────────────────
    if args.status:
        print(f"\nDatabase: {db._path}")
        print(f"Total rows: {db.count()}\n")
        in_db  = set(db.list_fluids())
        in_list = {canonical_of(name) for name, *_ in PROPELLANT_LIST}
        print(f"{'Fluid (canonical)':<45} {'In DB':>6}  {'Rows':>6}  T range @ 1 atm")
        print("-" * 80)
        for (name, T_min, T_max, notes) in PROPELLANT_LIST:
            c   = canonical_of(name)
            n   = db.count(c)
            rng = db.get_T_range(c, P_atm)
            rng_str = f"{rng[0]:.0f}–{rng[1]:.0f} K" if rng else "—"
            flag = "✓" if c in in_db else " "
            print(f"  {flag} {c:<43} {flag:>6}  {n:>6}  {rng_str}")
        unfilled = [canonical_of(n) for n, *_ in PROPELLANT_LIST
                    if canonical_of(n) not in in_db]
        if unfilled:
            print(f"\n  Not yet in DB ({len(unfilled)}): {unfilled}")
        print()
        return

    # ── Filter to single fluid if requested ───────────────────────────────────
    plan = PROPELLANT_LIST
    if args.fluid:
        target = canonical_of(args.fluid)
        plan   = [(n, lo, hi, note) for (n, lo, hi, note) in plan
                  if canonical_of(n) == target]
        if not plan:
            print(f"'{args.fluid}' not found in propellant list. "
                  f"Use --status to see available fluids.")
            sys.exit(1)

    # ── Dry run ───────────────────────────────────────────────────────────────
    if args.dry_run:
        print(f"\nDry run — T step = {args.T_step} K  "
              f"Pressures = [1 atm, {Pc/1e3:.0f} kPa]\n")
        print(f"  {'Fluid':<45} {'T_min':>7} {'T_max':>7} {'Points':>8}  Notes")
        print("  " + "-"*85)
        total_pts = 0
        for (name, T_min, T_max, notes) in plan:
            n_T    = int((T_max - T_min) / args.T_step) + 1
            n_pts  = n_T * 2   # two pressures
            total_pts += n_pts
            print(f"  {canonical_of(name):<45} {T_min:>7} {T_max:>7} {n_pts:>8}  {notes}")
        print(f"\n  Total points: {total_pts}  "
              f"(~{total_pts*200/1024:.0f} KB estimated DB size)")
        return

    # ── Fill ──────────────────────────────────────────────────────────────────
    total_inserted = total_skipped = total_failed = 0
    all_failures: dict = {}
    t0 = time.time()

    print(f"\nFilling {len(plan)} fluid(s)  "
          f"T step={args.T_step} K  "
          f"Pressures=[1 atm, {Pc/1e3:.0f} kPa]")
    if args.refill:
        print("  Mode: REFILL (existing rows will be deleted and replaced)")
    print()

    for (name, T_min, T_max, notes) in plan:
        canon = canonical_of(name)
        n_T   = int((T_max - T_min) / args.T_step) + 1
        print(f"{'─'*60}")
        print(f"  {canon}  [{T_min}–{T_max} K, {n_T} points × 2 pressures]")
        print(f"  ({notes})")

        summary = fill_fluid(
            canonical  = name,
            T_min      = T_min,
            T_max      = T_max,
            T_step     = args.T_step,
            P_list     = P_range,
            refill     = args.refill,
            verbose    = True,
        )
        total_inserted += summary["inserted"]
        total_skipped  += summary["skipped"]
        total_failed   += summary["failed"]
        if summary.get("failure_log"):
            all_failures[canon] = summary["failure_log"]

    elapsed = time.time() - t0

    print(f"\n{'='*60}")
    print(f"Complete in {elapsed:.1f}s")
    print(f"  Inserted : {total_inserted:>7}")
    print(f"  Skipped  : {total_skipped:>7}  (already in DB)")
    print(f"  Failed   : {total_failed:>7}")
    print(f"  DB rows  : {db.count():>7}")
    print(f"  DB size  : {os.path.getsize(db._path)/1024:.1f} KB")
    print(f"  DB path  : {db._path}")

    if all_failures:
        print(f"\n── Fluids that could not be auto-filled ──")
        for fluid_name, log in all_failures.items():
            total_pts = sum(len(v) for v in log.values())
            print(f"  {fluid_name}: {total_pts} failed points")
            for err_msg in list(log.keys())[:2]:
                print(f"    {err_msg[:90]}")
        print()
        print("  Add these manually via:")
        print("    from fluid_db import get_db")
        print("    db = get_db()")
        print("    db.insert('FluidName', T_K, P_Pa, props_dict, source='manual')")


def canonical_of(name: str) -> str:
    """Resolve a name to its canonical form using fluid_names if available."""
    try:
        from fluid_names import canonical
        return canonical(name)
    except Exception:
        return name


if __name__ == "__main__":
    main()
