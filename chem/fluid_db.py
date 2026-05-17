# =============================================================================
# fluid_db.py
# SQLite database manager for fluid thermodynamic properties.
#
# REPLACES fluid_dict.py
# ----------------------
# fluid_dict.py stored all property data as a Python dict that lived entirely
# in process memory. For small datasets that is fine, but as the number of
# fluids and temperature/pressure points grows the dict can become large
# enough to meaningfully affect the memory footprint of Engine_Class.
#
# This module stores the same data in a SQLite file (fluid_props.db) that
# sits on disk next to the Python files. Only the rows actually needed by a
# given calculation are ever loaded into memory — typically a few hundred rows
# for one heat-transfer run, regardless of how many millions of rows the DB
# contains.
#
# SQLite is part of Python's standard library (sqlite3). No extra dependencies.
#
# DATABASE FILE
# -------------
# Default location: same directory as this file, named "fluid_props.db".
# Override by setting the environment variable FLUID_DB_PATH, or by calling
# FluidDB(path="...") directly.
#
# SCHEMA
# ------
# Table: fluid_properties
#   fluid_canonical          TEXT    canonical fluid name (e.g. "Oxygen")
#   T_K                      REAL    temperature (K), rounded to 4 dp
#   P_Pa                     REAL    pressure (Pa), rounded to 2 dp
#   density_kg_m3            REAL    kg/m³
#   dynamic_viscosity_Pa_s   REAL    Pa·s
#   specific_heat_J_kgK      REAL    J/(kg·K)
#   thermal_conductivity_W_mK REAL   W/(m·K)
#   prandtl_number           REAL    dimensionless
#   speed_of_sound_m_s       REAL    m/s        (nullable)
#   quality                  REAL    0=liq,1=vap,-1=super (nullable)
#   enthalpy_J_kg            REAL    J/kg       (nullable)
#   entropy_J_kgK            REAL    J/(kg·K)   (nullable)
#   source                   TEXT    "coolprop"|"rocketprops"|"manual"
#   notes                    TEXT    freeform annotation
#   PRIMARY KEY (fluid_canonical, T_K, P_Pa)
#
# The primary key is automatically indexed by SQLite as a B-tree, giving
# O(log n) lookup. For any realistic dataset (<500 k rows) this is
# indistinguishable from constant time.
#
# USAGE
# -----
#   from fluid_db import get_db
#
#   db = get_db()                          # singleton — same connection reused
#   row = db.lookup("Oxygen", 90.0, 101325.0)   # returns dict or None
#   db.insert("Oxygen", 90.0, 101325.0, props_dict, source="coolprop")
#   db.close()                             # call at program exit if desired
#
# THREAD SAFETY
# -------------
# SQLite connections are not thread-safe by default. If Engine_Class is ever
# run in parallel threads, pass check_same_thread=False to sqlite3.connect()
# and wrap writes in a threading.Lock(). For the current single-threaded use
# case the defaults are fine.
# =============================================================================

from __future__ import annotations
import sqlite3
import os
from typing import Optional

# ---------------------------------------------------------------------------
# Property columns stored in the database.
# The first three (fluid_canonical, T_K, P_Pa) form the primary key.
# ---------------------------------------------------------------------------
_PROP_COLUMNS = [
    "density_kg_m3",
    "dynamic_viscosity_Pa_s",
    "specific_heat_J_kgK",
    "thermal_conductivity_W_mK",
    "prandtl_number",
    "speed_of_sound_m_s",
    "quality",
    "enthalpy_J_kg",
    "entropy_J_kgK",
]

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

# A covering index on T_K and P_Pa speeds up range queries used by the
# nearest-neighbour fallback (e.g. "find the closest T_K to 312.5 K").
_CREATE_INDEX = """
CREATE INDEX IF NOT EXISTS idx_fluid_T_P
ON fluid_properties (fluid_canonical, T_K, P_Pa);
"""


# ---------------------------------------------------------------------------
# FluidDB class
# ---------------------------------------------------------------------------

class FluidDB:
    """
    Thin wrapper around a SQLite connection for fluid property storage.

    Instantiate via the module-level get_db() singleton rather than directly.
    """

    def __init__(self, path: str) -> None:
        self._path = path
        self._conn = sqlite3.connect(path)
        self._conn.row_factory = sqlite3.Row   # rows accessible as dicts
        self._conn.execute("PRAGMA journal_mode=WAL;")   # safe concurrent reads
        self._conn.execute("PRAGMA synchronous=NORMAL;") # faster writes, still safe
        self._conn.execute(_CREATE_TABLE)
        self._conn.execute(_CREATE_INDEX)
        self._conn.commit()

    # ── Read ─────────────────────────────────────────────────────────────────

    def lookup(self, fluid: str, T_K: float, P_Pa: float) -> Optional[dict]:
        """
        Exact lookup for (fluid, T, P). Returns a property dict or None.

        This is the primary lookup path — one indexed B-tree seek.
        T_K and P_Pa are rounded to 4 and 2 decimal places respectively,
        matching the rounding applied during insert().
        """
        row = self._conn.execute(
            "SELECT * FROM fluid_properties "
            "WHERE fluid_canonical=? AND T_K=? AND P_Pa=?",
            (fluid, round(T_K, 4), round(P_Pa, 2))
        ).fetchone()
        return dict(row) if row is not None else None

    def lookup_nearest(self, fluid: str, T_K: float, P_Pa: float,
                       tol_T: float = 1.0, tol_P: float = 5000.0) -> Optional[dict]:
        """
        Find the closest stored (T, P) point within the given tolerances.

        Used as a fallback when an exact match is not found and live
        CoolProp/rocketprops calls are not available. Returns None if no
        point is within tolerance.

        The query uses the index on (fluid_canonical, T_K, P_Pa) to restrict
        the candidate set to the tolerance window before sorting.
        """
        row = self._conn.execute(
            "SELECT * FROM fluid_properties "
            "WHERE fluid_canonical=? "
            "  AND T_K  BETWEEN ? AND ? "
            "  AND P_Pa BETWEEN ? AND ? "
            "ORDER BY (ABS(T_K - ?) + ABS(P_Pa - ?) / 1e5) ASC "
            "LIMIT 1",
            (fluid,
             T_K - tol_T,  T_K + tol_T,
             P_Pa - tol_P, P_Pa + tol_P,
             T_K, P_Pa)
        ).fetchone()
        return dict(row) if row is not None else None

    def get_T_range(self, fluid: str, P_Pa: float) -> tuple[float, float] | None:
        """
        Return (T_min, T_max) of stored temperature points for a fluid at
        a given pressure. Useful for checking whether the DB covers the
        operating range before committing to a fill run.
        """
        row = self._conn.execute(
            "SELECT MIN(T_K), MAX(T_K) FROM fluid_properties "
            "WHERE fluid_canonical=? AND P_Pa=?",
            (fluid, round(P_Pa, 2))
        ).fetchone()
        if row and row[0] is not None:
            return (row[0], row[1])
        return None

    def list_fluids(self) -> list[str]:
        """Return the list of fluid names that have at least one row in the DB."""
        rows = self._conn.execute(
            "SELECT DISTINCT fluid_canonical FROM fluid_properties "
            "ORDER BY fluid_canonical"
        ).fetchall()
        return [r[0] for r in rows]

    def count(self, fluid: Optional[str] = None) -> int:
        """Return the number of rows for a fluid (or total if fluid is None)."""
        if fluid:
            return self._conn.execute(
                "SELECT COUNT(*) FROM fluid_properties WHERE fluid_canonical=?",
                (fluid,)
            ).fetchone()[0]
        return self._conn.execute(
            "SELECT COUNT(*) FROM fluid_properties"
        ).fetchone()[0]

    # ── Write ────────────────────────────────────────────────────────────────

    def insert(self, fluid: str, T_K: float, P_Pa: float,
               props: dict, source: str = "manual",
               notes: str = "", overwrite: bool = False) -> None:
        """
        Insert or update a property row.

        Parameters
        ----------
        fluid     : canonical fluid name.
        T_K       : temperature (K).
        P_Pa      : pressure (Pa).
        props     : dict of property_key → value. Keys not in _PROP_COLUMNS
                    are silently ignored.
        source    : "coolprop" | "rocketprops" | "manual".
        notes     : freeform string stored alongside the data.
        overwrite : if True, replace an existing (fluid, T, P) row.
                    if False (default), skip silently if the row exists.
        """
        T_r = round(T_K, 4)
        P_r = round(P_Pa, 2)

        if not overwrite:
            exists = self._conn.execute(
                "SELECT 1 FROM fluid_properties "
                "WHERE fluid_canonical=? AND T_K=? AND P_Pa=?",
                (fluid, T_r, P_r)
            ).fetchone()
            if exists:
                return   # already have this point — skip

        values = {col: props.get(col) for col in _PROP_COLUMNS}
        self._conn.execute(
            f"""
            INSERT OR REPLACE INTO fluid_properties
                (fluid_canonical, T_K, P_Pa,
                 {', '.join(_PROP_COLUMNS)}, source, notes)
            VALUES
                (?, ?, ?,
                 {', '.join(['?']*len(_PROP_COLUMNS))}, ?, ?)
            """,
            [fluid, T_r, P_r] +
            [values[c] for c in _PROP_COLUMNS] +
            [source, notes]
        )

    def insert_many(self, rows: list[tuple], commit: bool = True) -> int:
        """
        Bulk insert for efficiency when filling large temperature grids.

        Each element of rows must be a tuple:
          (fluid, T_K, P_Pa, props_dict, source, notes)

        Returns the number of rows inserted.
        """
        inserted = 0
        for (fluid, T_K, P_Pa, props, source, notes) in rows:
            self.insert(fluid, T_K, P_Pa, props, source, notes)
            inserted += 1
        if commit:
            self._conn.commit()
        return inserted

    def commit(self) -> None:
        """Flush pending writes to disk."""
        self._conn.commit()

    def delete_fluid(self, fluid: str) -> int:
        """Delete all rows for a fluid. Returns the number of rows deleted."""
        cur = self._conn.execute(
            "DELETE FROM fluid_properties WHERE fluid_canonical=?", (fluid,))
        self._conn.commit()
        return cur.rowcount

    def close(self) -> None:
        """Close the database connection."""
        self._conn.close()

    def __repr__(self) -> str:
        return f"FluidDB(path={self._path!r}, rows={self.count()})"


# ---------------------------------------------------------------------------
# Module-level singleton
# ---------------------------------------------------------------------------
# A single FluidDB instance is shared across all callers within a Python
# session.  This avoids the overhead of opening and closing the SQLite file
# on every property lookup.  The connection is opened on first access and
# remains open until the process exits (or close() is called explicitly).

_DB_INSTANCE: Optional[FluidDB] = None


def get_db(path: Optional[str] = None) -> FluidDB:
    """
    Return the module-level FluidDB singleton, creating it if necessary.

    The database file is located by this priority order:
      1. The `path` argument if provided.
      2. The FLUID_DB_PATH environment variable.
      3. fluid_props.db in the same directory as this file (default).

    The connection is opened once and reused for the lifetime of the
    Python session.  Calling get_db() multiple times always returns the
    same instance unless close_db() has been called first.
    """
    global _DB_INSTANCE
    if _DB_INSTANCE is not None:
        return _DB_INSTANCE

    if path is None:
        path = os.environ.get(
            "FLUID_DB_PATH",
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "fluid_props.db")
        )

    _DB_INSTANCE = FluidDB(path)
    return _DB_INSTANCE


def close_db() -> None:
    """Close the singleton connection and allow it to be re-opened."""
    global _DB_INSTANCE
    if _DB_INSTANCE is not None:
        _DB_INSTANCE.close()
        _DB_INSTANCE = None
