"""
fast_fill_fluid_db.py
=====================
High-speed bulk population of fluid_props.db using CoolProp's AbstractState
interface and SQLite executemany() for maximum insert throughput, with
automatic failover through three property sources.

SOURCE PRIORITY ORDER (per fluid, per T/P point)
-------------------------------------------------
1. CoolProp AbstractState ("HEOS" backend)
     Fastest. One AbstractState object per fluid, updated per (T, P).
     Covers most pure fluids and mixtures with full EOS accuracy.
     Fails for: ClF3, ClF5, Gasoline, Jet-A, MON blends, custom fluids.

2. rocketprops
     Pre-computed polynomial curve fits. Covers LOX, RP-1, MMH, N2H4,
     N2O4, H2O2, A50, and other common rocket propellants.
     Called only if CoolProp AbstractState init fails for the fluid
     OR if CoolProp returned None for a specific property at a (T, P).

3. Custom equations file (fluid_equations.py)
     A file you create and populate with Python functions for fluids
     not covered by CoolProp or rocketprops — literature correlations,
     supplier data fits, or hand-coded equations of state.
     Called only if both CoolProp and rocketprops fail or return None.
     See the CUSTOM EQUATIONS section below for the required interface.

FAILOVER LOGIC
--------------
At the fluid level:
  - Try CoolProp AbstractState init.
  - If that fails, try rocketprops init.
  - If that fails, try fluid_equations.py.
  - If all three fail, skip the fluid and log it.

At the property level (per T/P point):
  - CoolProp fills all properties it can return.
  - Any property that came back as None is re-attempted from rocketprops.
  - Any property still None after rocketprops is re-attempted from equations.
  - The final source label written to the 'source' column reflects which
    source actually provided the majority of the data.

CUSTOM EQUATIONS FILE (fluid_equations.py)
------------------------------------------
Create fluid_equations.py in the same directory as this script.
It must define one function per fluid, named get_<canonical>_props,
that accepts (T_K, P_Pa) and returns a dict with any subset of the
standard property keys. Properties not returned are left as None (NULL).

Example fluid_equations.py:

    def get_ClF3_props(T_K: float, P_Pa: float) -> dict:
        \"\"\"
        Chlorine trifluoride — correlations from:
        Gmelin Handbook of Inorganic Chemistry, ClF3, 8th ed.
        Valid: 200–400 K, 0.1–5 MPa
        \"\"\"
        # Density: linear fit to experimental data
        rho = 1950.0 - 0.85 * (T_K - 200.0)
        # Viscosity: Andrade equation
        mu  = 1.2e-3 * (200.0 / T_K) ** 1.8
        # Specific heat: constant (weak T dependence)
        Cp  = 920.0
        # Thermal conductivity: linear fit
        k   = 0.195 - 0.0002 * (T_K - 200.0)
        Pr  = (mu * Cp) / k if k else None
        return {
            "density_kg_m3":             rho,
            "dynamic_viscosity_Pa_s":    mu,
            "specific_heat_J_kgK":       Cp,
            "thermal_conductivity_W_mK": k,
            "prandtl_number":            Pr,
        }

    def get_Gasoline_props(T_K: float, P_Pa: float) -> dict:
        \"\"\"Gasoline approximated as iso-octane (C8H18).\"\"\"
        rho = 720.0 - 0.6 * (T_K - 298.15)
        mu  = 5.5e-4 * (298.15 / T_K) ** 2.1
        Cp  = 2090.0 + 3.5 * (T_K - 298.15)
        k   = 0.132 - 0.00012 * (T_K - 298.15)
        Pr  = (mu * Cp) / k if k else None
        return {
            "density_kg_m3":             rho,
            "dynamic_viscosity_Pa_s":    mu,
            "specific_heat_J_kgK":       Cp,
            "thermal_conductivity_W_mK": k,
            "prandtl_number":            Pr,
        }

The function name convention is:  get_<canonical_name>_props
where <canonical_name> has spaces replaced with underscores and
hyphens replaced with underscores.  Examples:
    canonical "ClF3"        → get_ClF3_props
    canonical "Nitric Acid" → get_Nitric_Acid_props
    canonical "Jet-A"       → get_Jet_A_props

DATABASE COMPATIBILITY
----------------------
Writes to fluid_props.db in the schema that fluid_db.py creates.
Table: fluid_properties  PRIMARY KEY (fluid_canonical, T_K, P_Pa)

USAGE
-----
    python fast_fill_fluid_db.py
    python fast_fill_fluid_db.py --fluid ClF3
    python fast_fill_fluid_db.py --fluid Ethanol --T-min 160 --T-max 500
    python fast_fill_fluid_db.py --overwrite
    python fast_fill_fluid_db.py --status
    python fast_fill_fluid_db.py --dry-run
"""

import sqlite3
import os
import sys
import time
import argparse
import importlib.util
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Optional imports — each is used only if available
# ---------------------------------------------------------------------------

# CoolProp
try:
    import CoolProp.CoolProp as CP
    _COOLPROP_AVAILABLE = True
except ImportError:
    _COOLPROP_AVAILABLE = False

# rocketprops
try:
    from rocketprops.rocket_prop import get_prop as _rp_get_prop
    _ROCKETPROPS_AVAILABLE = True
    # Unit conversion constants
    _BTU_lbmR_to_J_kgK  = 4186.8
    _BTU_hrftR_to_W_mK  = 1.730735
    _cP_to_Pa_s         = 1e-3
    _SG_to_kg_m3        = 1000.0
except ImportError:
    _ROCKETPROPS_AVAILABLE = False

# Custom equations file
def _load_equations_module():
    """
    Load fluid_equations.py from the same directory as this script.
    Returns the module object or None if the file is not found.
    Called once at startup; result cached in _EQUATIONS_MOD.
    """
    eq_path = os.path.join(_HERE, 'fluid_equations.py')
    if not os.path.exists(eq_path):
        return None
    try:
        spec = importlib.util.spec_from_file_location('fluid_equations', eq_path)
        mod  = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod
    except Exception as e:
        print(f"[fast_fill] WARNING: could not load fluid_equations.py — {e}")
        return None

_EQUATIONS_MOD = _load_equations_module()


# ---------------------------------------------------------------------------
# fluid_names translation (optional but strongly recommended)
# ---------------------------------------------------------------------------

def _load_fluid_names():
    """Load fluid_names.py and return (translate_fn, canonical_fn, FLUID_NAME_MAP)."""
    fn_path = os.path.join(_HERE, 'fluid_names.py')
    if not os.path.exists(fn_path):
        return None, None, {}
    try:
        spec = importlib.util.spec_from_file_location('fluid_names', fn_path)
        mod  = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod.translate, mod.canonical, mod.FLUID_NAME_MAP
    except Exception as e:
        print(f"[fast_fill] WARNING: could not load fluid_names.py — {e}")
        return None, None, {}

_fn_translate, _fn_canonical, _FLUID_NAME_MAP = _load_fluid_names()

def _translate(name, target):
    if _fn_translate:
        try:    return _fn_translate(name, target)
        except: return None
    return name

def _canonical(name):
    if _fn_canonical:
        try:    return _fn_canonical(name)
        except: return name
    return name


# ---------------------------------------------------------------------------
# Database schema — must match fluid_db.py exactly
# ---------------------------------------------------------------------------

_CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS fluid_properties (
    fluid_canonical           TEXT    NOT NULL,
    T_K                       REAL    NOT NULL,
    P_Pa                      REAL    NOT NULL,
    density_kg_m3             REAL,
    dynamic_viscosity_Pa_s    REAL,
    specific_heat_J_kgK       REAL,
    thermal_conductivity_W_mK REAL,
    prandtl_number            REAL,
    speed_of_sound_m_s        REAL,
    quality                   REAL,
    enthalpy_J_kg             REAL,
    entropy_J_kgK             REAL,
    source                    TEXT,
    notes                     TEXT,
    PRIMARY KEY (fluid_canonical, T_K, P_Pa)
);
"""

_CREATE_INDEX = """
CREATE INDEX IF NOT EXISTS idx_fluid_T_P
ON fluid_properties (fluid_canonical, T_K, P_Pa);
"""

_CREATE_INDEX2 = """
CREATE INDEX IF NOT EXISTS idx_fluid_P_T
ON fluid_properties (fluid_canonical, P_Pa, T_K);
"""

_COLUMNS = [
    "fluid_canonical", "T_K", "P_Pa",
    "density_kg_m3", "dynamic_viscosity_Pa_s", "specific_heat_J_kgK",
    "thermal_conductivity_W_mK", "prandtl_number", "speed_of_sound_m_s",
    "quality", "enthalpy_J_kg", "entropy_J_kgK",
    "source", "notes",
]

_PROP_KEYS = [
    "density_kg_m3", "dynamic_viscosity_Pa_s", "specific_heat_J_kgK",
    "thermal_conductivity_W_mK", "prandtl_number", "speed_of_sound_m_s",
    "quality", "enthalpy_J_kg", "entropy_J_kgK",
]

def open_db(path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute("PRAGMA cache_size=-524288;")   # 512 MB page cache
    conn.execute("PRAGMA mmap_size=4294967296;") # 4 GB mmap
    conn.execute("PRAGMA temp_store=MEMORY;")
    conn.execute(_CREATE_TABLE)
    conn.execute(_CREATE_INDEX)
    conn.execute(_CREATE_INDEX2)
    conn.commit()
    return conn


# ---------------------------------------------------------------------------
# Per-source property fetchers
# ---------------------------------------------------------------------------

class CoolPropSource:
    """
    Wraps a single CoolProp AbstractState object for one fluid.
    Initialised once; state.update() called per (T, P) point.
    """
    def __init__(self, cp_name: str):
        if not _COOLPROP_AVAILABLE:
            raise ImportError("CoolProp is not installed.")
        self._state = CP.AbstractState("HEOS", cp_name)
        self._name  = cp_name

    def get(self, T_K: float, P_Pa: float) -> dict:
        """
        Return a property dict for (T_K, P_Pa).
        Properties that CoolProp cannot provide at this state point are None.
        Raises if state.update() itself fails (T/P outside EOS validity).
        """
        self._state.update(CP.PT_INPUTS, P_Pa, T_K)

        def _g(fn):
            try:    return fn()
            except: return None

        rho = _g(self._state.rhomass)
        mu  = _g(self._state.viscosity)
        Cp  = _g(self._state.cpmass)
        k   = _g(self._state.conductivity)
        a   = _g(self._state.speed_sound)
        Q   = _g(self._state.Q)
        H   = _g(self._state.hmass)
        S   = _g(self._state.smass)
        Pr  = (mu * Cp / k) if (mu and Cp and k) else None

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
        }


class RocketPropsSource:
    """
    Wraps a rocketprops PropObj for one fluid.
    Initialised once; property methods called per (T, P) point.
    rocketprops does not support mixtures or two-phase lookups.
    """
    def __init__(self, rp_name: str):
        if not _ROCKETPROPS_AVAILABLE:
            raise ImportError("rocketprops is not installed.")
        self._prop = _rp_get_prop(rp_name)
        self._name = rp_name

    def get(self, T_K: float, P_Pa: float) -> dict:
        """
        Return a property dict for (T_K, P_Pa).
        rocketprops only provides rho, mu, Cp, k — speed of sound,
        quality, enthalpy, entropy are returned as None.
        """
        P_atm = P_Pa / 101325.0

        def _g(fn, *args):
            try:    return fn(*args)
            except: return None

        rho = _g(self._prop.SG_compressed, T_K, P_atm)
        if rho is not None:
            rho *= _SG_to_kg_m3
        mu  = _g(self._prop.viscosity, T_K)
        if mu is not None:
            mu *= _cP_to_Pa_s
        Cp  = _g(self._prop.specific_heat, T_K)
        if Cp is not None:
            Cp *= _BTU_lbmR_to_J_kgK
        k   = _g(self._prop.thermal_conductivity, T_K)
        if k is not None:
            k *= _BTU_hrftR_to_W_mK
        Pr  = (mu * Cp / k) if (mu and Cp and k) else None

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
        }


class EquationsSource:
    """
    Calls a user-defined function from fluid_equations.py.

    The function must be named  get_<canonical>_props(T_K, P_Pa) -> dict
    where <canonical> is the fluid's canonical name with spaces and hyphens
    replaced by underscores.

    Examples:
        canonical "ClF3"        → get_ClF3_props
        canonical "Nitric Acid" → get_Nitric_Acid_props
        canonical "Jet-A"       → get_Jet_A_props
    """
    def __init__(self, canonical: str):
        if _EQUATIONS_MOD is None:
            raise FileNotFoundError(
                "fluid_equations.py not found in the script directory.")
        fn_name = "get_" + canonical.replace(" ", "_").replace("-", "_") + "_props"
        self._fn = getattr(_EQUATIONS_MOD, fn_name, None)
        if self._fn is None:
            raise AttributeError(
                f"fluid_equations.py has no function '{fn_name}' for '{canonical}'.")
        self._name = canonical

    def get(self, T_K: float, P_Pa: float) -> dict:
        """
        Call the user function and return its dict, filling missing keys with None.
        """
        try:
            result = self._fn(T_K, P_Pa)
        except Exception as e:
            raise RuntimeError(
                f"fluid_equations.py: {self._name} @ T={T_K}K P={P_Pa}Pa — {e}")
        # Normalise: ensure all standard keys are present
        return {k: result.get(k) for k in _PROP_KEYS}


# ---------------------------------------------------------------------------
# Source initialisation with failover
# ---------------------------------------------------------------------------

def _init_sources(canon: str) -> tuple[list, str]:
    """
    Attempt to initialise each source for a fluid in priority order.
    Returns (list_of_active_sources, description_string).

    Each source in the list is a (source_obj, label) tuple.
    Sources that fail to initialise are skipped with a printed reason.
    """
    sources  = []
    skipped  = []

    # ── 1. CoolProp ──────────────────────────────────────────────────────────
    cp_name = _translate(canon, 'coolprop')
    if cp_name:
        try:
            sources.append((CoolPropSource(cp_name), 'coolprop'))
        except Exception as e:
            skipped.append(f"CoolProp({cp_name}): {e}")
    else:
        skipped.append("CoolProp: no entry in fluid_names")

    # ── 2. rocketprops ───────────────────────────────────────────────────────
    rp_name = _translate(canon, 'rocketprops')
    if rp_name:
        try:
            sources.append((RocketPropsSource(rp_name), 'rocketprops'))
        except Exception as e:
            skipped.append(f"rocketprops({rp_name}): {e}")
    else:
        skipped.append("rocketprops: no entry in fluid_names")

    # ── 3. Custom equations ──────────────────────────────────────────────────
    try:
        sources.append((EquationsSource(canon), 'equations'))
    except (FileNotFoundError, AttributeError) as e:
        skipped.append(f"equations: {e}")
    except Exception as e:
        skipped.append(f"equations: unexpected error — {e}")

    desc = f"{len(sources)} source(s) active"
    if skipped:
        desc += f"  [{'; '.join(skipped)}]"
    return sources, desc


# ---------------------------------------------------------------------------
# Property merging — combine results from multiple sources
# ---------------------------------------------------------------------------

def _merge(primary: dict, fallback: dict) -> dict:
    """
    Fill any None values in `primary` with the corresponding value from
    `fallback`.  Returns a new dict; neither argument is modified.
    """
    return {k: (primary[k] if primary[k] is not None else fallback.get(k))
            for k in _PROP_KEYS}


# ---------------------------------------------------------------------------
# Core fill function
# ---------------------------------------------------------------------------

DB_PATH = os.environ.get(
    "FLUID_DB_PATH",
    os.path.join(_HERE, "fluid_props.db")
)

def _fill_fluid(conn: sqlite3.Connection,
                fluid_canonical: str,
                T_arr: np.ndarray,
                P_arr: np.ndarray,
                overwrite: bool = False,
                verbose: bool = True) -> dict:
    """
    Fill one fluid across the full (T, P) grid with automatic source failover.

    Source selection per (T, P) point:
      1. CoolProp AbstractState — primary; fastest; widest coverage.
      2. rocketprops            — fallback if CoolProp init fails or
                                  returns None for specific properties.
      3. fluid_equations.py     — final fallback for custom/literature data.

    The 'source' column written to the DB names whichever source provided
    the most data for that row.

    Batch strategy: accumulate one pressure slice in memory, then flush with
    executemany() + one COMMIT. Peak memory is O(len(T_arr)) rows.
    """
    sql = "INSERT OR REPLACE" if overwrite else "INSERT OR IGNORE"
    sql += (f" INTO fluid_properties ({', '.join(_COLUMNS)}) "
            f"VALUES ({', '.join(['?']*len(_COLUMNS))})")

    # Initialise all available sources
    sources, src_desc = _init_sources(fluid_canonical)
    if verbose:
        print(f"  Sources: {src_desc}")

    if not sources:
        if verbose:
            print(f"  [SKIP] No property sources available for '{fluid_canonical}'")
        return {"inserted": 0, "failed_init": 1, "failed_update": 0,
                "total_points": 0, "sources_used": []}

    notes_str = f"fast_fill_fluid_db.py — sources: {[s[1] for s in sources]}"
    n_inserted = n_failed = 0
    sources_used: set[str] = set()
    t0 = time.time()

    for P in P_arr:
        P_r    = round(float(P), 2)
        batch: list[tuple] = []

        for T in T_arr:
            T_r = round(float(T), 4)

            # Attempt each source in priority order, merging results
            props: dict | None = None
            source_label = "unknown"

            for src_obj, src_label in sources:
                try:
                    candidate = src_obj.get(T_r, P_r)
                except Exception:
                    # This source failed at this (T, P) — try next source
                    continue

                if props is None:
                    # First successful source — use its result as the base
                    props        = candidate
                    source_label = src_label
                    sources_used.add(src_label)
                else:
                    # Subsequent source — fill any None values from the base
                    before_nones = sum(1 for v in props.values() if v is None)
                    props = _merge(props, candidate)
                    after_nones  = sum(1 for v in props.values() if v is None)
                    if before_nones != after_nones:
                        # This source contributed at least one property
                        source_label += f"+{src_label}"
                        sources_used.add(src_label)

                # If all core properties are filled, no need to try more sources
                core = ["density_kg_m3", "dynamic_viscosity_Pa_s",
                        "specific_heat_J_kgK", "thermal_conductivity_W_mK",
                        "prandtl_number"]
                if all(props.get(k) is not None for k in core):
                    break

            if props is None:
                # All sources failed at this (T, P) — skip the row
                n_failed += 1
                continue

            batch.append((
                fluid_canonical, T_r, P_r,
                props.get("density_kg_m3"),
                props.get("dynamic_viscosity_Pa_s"),
                props.get("specific_heat_J_kgK"),
                props.get("thermal_conductivity_W_mK"),
                props.get("prandtl_number"),
                props.get("speed_of_sound_m_s"),
                props.get("quality"),
                props.get("enthalpy_J_kg"),
                props.get("entropy_J_kgK"),
                source_label,
                notes_str,
            ))

        # Flush this pressure slice
        if batch:
            conn.executemany(sql, batch)
            conn.commit()
            n_inserted += len(batch)

        if verbose:
            elapsed = time.time() - t0
            rate    = n_inserted / elapsed * 60 if elapsed > 0 else 0
            print(f"    P={P_r/1e5:6.2f} bar  "
                  f"rows={n_inserted:>8,}  "
                  f"failed={n_failed:>5}  "
                  f"rate={rate:>8,.0f} rows/min",
                  end='\r', flush=True)

    if verbose:
        print()

    return {
        "inserted":      n_inserted,
        "failed_init":   0,
        "failed_update": n_failed,
        "total_points":  len(T_arr) * len(P_arr),
        "sources_used":  sorted(sources_used),
    }


# ---------------------------------------------------------------------------
# Fluid list resolution
# ---------------------------------------------------------------------------

def _resolve_fluid_list(requested: list[str]) -> list[str]:
    """Resolve a list of names to canonical names."""
    return [_canonical(n) for n in requested]


def _all_fluids() -> list[str]:
    """
    Return all canonical fluid names from FLUID_NAME_MAP plus any fluid
    that has a function in fluid_equations.py but is not in fluid_names.
    """
    names = list(_FLUID_NAME_MAP.keys())

    # Also pick up any fluid defined only in fluid_equations.py
    if _EQUATIONS_MOD:
        import inspect
        for fn_name, fn in inspect.getmembers(_EQUATIONS_MOD, inspect.isfunction):
            if fn_name.startswith("get_") and fn_name.endswith("_props"):
                # Reverse the name mangling to get the canonical name
                middle = fn_name[4:-6]   # strip "get_" and "_props"
                # We can't perfectly reverse underscores back to spaces/hyphens
                # so just add it as-is if not already known
                if middle not in names and middle not in [_canonical(n) for n in names]:
                    names.append(middle)

    return names


# ---------------------------------------------------------------------------
# Status report
# ---------------------------------------------------------------------------

def _status(conn: sqlite3.Connection) -> None:
    total = conn.execute("SELECT COUNT(*) FROM fluid_properties").fetchone()[0]
    print(f"\nDatabase: {DB_PATH}")
    print(f"Total rows: {total:,}\n")

    rows = conn.execute(
        "SELECT fluid_canonical, COUNT(*), MIN(T_K), MAX(T_K), "
        "       MIN(P_Pa), MAX(P_Pa), "
        "       GROUP_CONCAT(DISTINCT source) "
        "FROM fluid_properties GROUP BY fluid_canonical ORDER BY fluid_canonical"
    ).fetchall()

    if not rows:
        print("  (empty)")
        return

    print(f"  {'Fluid':<45} {'Rows':>8}  {'T range (K)':^16}  "
          f"{'P range (bar)':^16}  Sources")
    print("  " + "─" * 105)
    for (fluid, n, T_min, T_max, P_min, P_max, srcs) in rows:
        print(f"  {fluid:<45} {n:>8,}  "
              f"{T_min:.0f}–{T_max:.0f} K{'':<6}  "
              f"{P_min/1e5:.1f}–{P_max/1e5:.1f} bar{'':<3}  "
              f"{srcs}")
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "High-speed bulk fill of fluid_props.db. "
            "Tries CoolProp → rocketprops → fluid_equations.py per fluid."
        )
    )
    parser.add_argument("--fluid",     nargs="+", default=None,
                        help="One or more fluid names. Default: all fluids in FLUID_NAME_MAP.")
    parser.add_argument("--T-min",     type=float, default=20.0)
    parser.add_argument("--T-max",     type=float, default=600.0)
    parser.add_argument("--T-step",    type=float, default=1.0)
    parser.add_argument("--P-min",     type=float, default=1e5)
    parser.add_argument("--P-max",     type=float, default=10e6)
    parser.add_argument("--P-step",    type=float, default=1000.0)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--status",    action="store_true")
    parser.add_argument("--dry-run",   action="store_true")
    args = parser.parse_args()

    conn = open_db(DB_PATH)

    if args.status:
        _status(conn)
        conn.close()
        return

    T_arr = np.round(np.arange(args.T_min,
                               args.T_max + args.T_step * 0.5,
                               args.T_step), 4)
    P_arr = np.round(np.arange(args.P_min,
                               args.P_max + args.P_step * 0.5,
                               args.P_step), 2)

    fluid_list = (_resolve_fluid_list(args.fluid)
                  if args.fluid else _all_fluids())

    n_T = len(T_arr); n_P = len(P_arr)

    print(f"\nFill plan")
    print(f"  Fluids      : {len(fluid_list)}")
    print(f"  T grid      : {T_arr[0]:.1f}–{T_arr[-1]:.1f} K  "
          f"step={args.T_step} K  ({n_T} pts)")
    print(f"  P grid      : {P_arr[0]/1e5:.1f}–{P_arr[-1]/1e5:.1f} bar  "
          f"step={args.P_step/1e3:.1f} kPa  ({n_P} pts)")
    print(f"  Total pts   : {n_T*n_P*len(fluid_list):,}")
    print(f"  Mode        : {'REPLACE' if args.overwrite else 'SKIP existing'}")
    print(f"  DB path     : {DB_PATH}")
    print(f"  CoolProp    : {'available' if _COOLPROP_AVAILABLE else 'NOT INSTALLED'}")
    print(f"  rocketprops : {'available' if _ROCKETPROPS_AVAILABLE else 'NOT INSTALLED'}")
    print(f"  equations   : {'fluid_equations.py loaded' if _EQUATIONS_MOD else 'fluid_equations.py not found'}")

    if args.dry_run:
        print("\nDry run — nothing inserted.\n")
        conn.close()
        return

    print()
    grand_inserted = grand_failed = 0
    grand_sources: dict[str, int] = {}
    t_global = time.time()

    for canon in fluid_list:
        print(f"{'─'*65}")
        print(f"  {canon}")
        print(f"  {n_T} T × {n_P} P = {n_T*n_P:,} points")

        t_fluid = time.time()
        result  = _fill_fluid(conn, canon, T_arr, P_arr,
                              overwrite=args.overwrite, verbose=True)
        elapsed = time.time() - t_fluid
        rate    = result["inserted"] / elapsed * 60 if elapsed > 0 else 0

        print(f"  Inserted : {result['inserted']:>10,}  "
              f"({rate:,.0f} rows/min)  "
              f"failed={result['failed_update']}  "
              f"sources={result['sources_used']}")

        grand_inserted += result["inserted"]
        grand_failed   += result["failed_update"]
        for s in result["sources_used"]:
            grand_sources[s] = grand_sources.get(s, 0) + 1

    total_elapsed = time.time() - t_global
    db_size = os.path.getsize(DB_PATH) / (1024 * 1024)

    print(f"\n{'='*65}")
    print(f"Complete in {total_elapsed:.1f}s  ({total_elapsed/60:.1f} min)")
    print(f"  Total inserted : {grand_inserted:>12,}")
    print(f"  Total failed   : {grand_failed:>12,}")
    print(f"  Sources used   : {dict(grand_sources)}")
    print(f"  DB size        : {db_size:.1f} MB")
    print(f"  DB path        : {DB_PATH}")

    conn.close()


if __name__ == "__main__":
    main()
