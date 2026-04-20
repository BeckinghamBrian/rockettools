# rockettools

A Python toolkit for liquid rocket engine preliminary design, thermodynamic analysis, regenerative cooling, and automated report generation. Provide a plain-text `.toml` input file and get back a formatted Word document containing your inputs, chamber outputs, and selected plots — no scripting required beyond a single command.

---

## Table of Contents

- [What's in this repo](#whats-in-this-repo)
- [Quick Start](#quick-start)
- [Dependencies](#dependencies)
- [File Reference](#file-reference)
  - [Engine_Class.py](#engine_classpy)
  - [generate_report.js](#generate_reportjs)
  - [fluid_dict.py](#fluid_dictpy)
  - [mat_dict.py](#mat_dictpy)
- [Input File Format](#input-file-format)
  - [Required Parameters](#required-parameters)
  - [Optional Parameters](#optional-parameters)
  - [Plot Selection](#plot-selection)
  - [Custom Propellants](#custom-propellants)
- [Report Output](#report-output)
- [Using the Dictionaries](#using-the-dictionaries)
  - [fluid_dict](#fluid_dict)
  - [mat_dict](#mat_dict)
- [Python API](#python-api)
- [Performance Notes](#performance-notes)

---

## What's in this repo

| File | Language | Purpose |
|------|----------|---------|
| `Engine_Class.py` | Python | Core engine class — CEA wrap, geometry, heat transfer, report generation |
| `generate_report.js` | Node.js | Docx builder — assembles the Word report from a JSON data file |
| `fluid_dict.py` | Python | thermodynamic property lookup for fluids |
| `mat_dict.py` | Python | thermomechanic property lookup for solid materials |
| `example_engine.toml` | TOML | Annotated example input file — start here |

---

## Quick Start

**1. Clone and install dependencies** (see [Dependencies](#dependencies) below).

**2. Copy and edit the example input file:**

```bash
cp example.toml my_engine.toml
# Edit my_engine.toml with your thrust, chamber pressure, OF ratio, etc.
```

**3. Run from the command line:**

```bash
python Engine_Class.py my_engine.toml ./output
```

That single command runs all analyses, saves the selected plots as PNGs, and writes a formatted `.docx` report to `./output/`. 

**Or use the Python API directly:**

```python
from Engine_Class import Engine

eng = Engine.from_file('my_engine.toml')
eng.generate_report()
```

---

## Dependencies

### Python packages

```bash
pip install rocketcea CoolProp matplotlib numpy pandas molmass
# Python 3.11+: tomllib is built-in
# Python < 3.11: pip install tomli
```

| Package | Purpose |
|---------|---------|
| `rocketcea` | NASA CEA wrapper — combustion thermodynamics |
| `CoolProp` | Coolant fluid properties (density, viscosity, Cp, k) |
| `matplotlib` | Plot generation |
| `numpy` | Vectorised calculations |
| `pandas` | Nozzle geometry export to CSV |
| `molmass` | Molecular mass calculation for custom fuel blends |
| `tomllib` / `tomli` | TOML input file parsing |
| `MOC_nozzle` | Creates a Method of Characteristics nozzle profile |

### Node.js

Required only for `.docx` report generation. The `docx` package must be available globally or in a local `node_modules`:

```bash
npm install -g docx
# or locally:
npm install docx
```

---

## File Reference

### Engine_Class.py

The main class. Instantiate it with constructor arguments or a TOML file; call analysis methods individually or call `generate_report()` to do everything at once.

**Key methods:**

| Method | What it does |
|--------|-------------|
| `getCeaProperties()` | Calls NASA CEA via rocketcea to get Tc, gamma, Isp, C*, Mach at exit |
| `engineDimensions()` | Computes all chamber and nozzle geometry from CEA outputs |
| `engineContour()` | Builds the full engine wall profile using numpy-vectorised piecewise curves and MOC for the nozzle |
| `propFlowRates()` | Computes mass flow rates, volume flow rates, and propellant tank volumes required for burn duration |
| `heatTransfer()` | Bartz-equation regenerative cooling analysis along the entire engine contour |
| `printChamberInfo()` | Prints a full summary of all thermodynamic and geometric outputs to stdout |
| `generate_report()` | Runs all analyses, saves selected plots, and generates the `.docx` report |
| `from_file(path)` | Class method — construct an Engine entirely from a TOML file |

---

### generate_report.js

Called automatically by `Engine_Class.generate_report()` via `subprocess`. You do not need to run this manually.

It reads `report_data.json` from the output directory and produces a `.docx` with three sections:

1. **Inputs** — a parameter table plus a complete TOML block (copy-paste to create a new input file)
2. **Chamber Info & Dimensions** — all thermodynamic and geometric outputs
3. **Plots** — selected plots embedded at 9 × 4 inches, two per US Letter page

---

### fluid_dict.py

A thermodynamic property store for coolant fluids, keyed by `(fluid_name, T_K, P_Pa)`. Uses nested Python dicts (hash maps) so both the fluid name lookup and the state-point lookup are constant time regardless of how many entries are stored.

The key format is intentionally identical to the `fuelComposition` string that `Engine_Class` already passes to CoolProp, so the same string can be used for both dictionary lookup and live CoolProp calls.

Includes `populate_from_coolprop()` to batch-fill the dictionary over a temperature/pressure grid at startup, and `get_fluid_props()` for nearest-state fallback when an exact key is not found.

---

### mat_dict.py

A constant-time thermomechanical property store for solid materials, keyed by `(material_name, T_K)`. 

Includes `get_mat_props()` with linear interpolation between registered temperature points, so you can query an arbitrary wall temperature that falls between your data points.

Properties stored per entry:

| Property | Key | Units |
|----------|-----|-------|
| Density | `density_kg_m3` | kg/m³ |
| Thermal conductivity | `thermal_conductivity_W_mK` | W/(m·K) |
| Coefficient of thermal expansion | `coef_thermal_expansion_1_K` | 1/K |
| Young's modulus | `youngs_modulus_GPa` | GPa |
| Yield strength | `yield_strength_MPa` | MPa |

Pre-populated materials: **Inconel 718** (298 K, 773 K, 973 K, 1073 K), **OFHC Copper** (298 K, 573 K), **304 Stainless Steel** (298 K, 773 K).

---

## Input File Format

The input file is [TOML](https://toml.io) — a minimal key=value format that maps directly to Python types with no custom parsing. Python 3.11+ reads it with the built-in `tomllib`; earlier versions need `pip install tomli`.

See `example_engine.toml` for a fully annotated example. The structure is explained below.

### Required Parameters

These five keys must be present in every input file:

```toml
OF                     = 2.1      # oxidiser-to-fuel mass ratio
Pc_psi                 = 300.0    # chamber pressure (psi)
thrust_lbf             = 500.0    # thrust (lbf) — converted to N internally
height_of_optimization = 0.0      # altitude for nozzle area-ratio optimisation (m)
burnTime               = 30.0     # burn duration (s)
```

### Optional Parameters

Every other parameter has a class default and only needs to appear in the file if you want to override it:

```toml
engine_name           = "Demo"
Lstar                 = 1.0         # characteristic chamber length (m)
alpha                 = 30          # converging half-angle (degrees)
beta                  = 15          # initial diverging half-angle (degrees)
num_characteristics   = 300         # MOC characteristic lines for nozzle contour

# Regenerative cooling
Kwall                 = 12.0        # wall thermal conductivity W/(m·K)  — Inconel 718
Kc                    = 0.17        # coolant thermal conductivity W/(m·K)
chamberWallThickness  = 1.5         # wall thickness (mm)
channelHeight         = 3.0         # channel depth (mm)
channelWallThickness  = 0.8         # wall between channels (mm)
numChannels           = 60
coefThermEx           = 1.3e-5      # CTE (1/K)
youngMod              = 200.0       # Young's modulus (GPa)
```

### Plot Selection

Every plot defaults to `true`. Set any flag to `false` to exclude it from the report:

```toml
plot_engine_contour   = true
plot_nozzle_contour   = true
plot_isp_vs_of        = false   # skip the Isp vs OF sweep
plot_mach             = true
plot_heat_flux        = true
plot_wall_temp        = true
plot_thermal_stress   = false   # skip thermal stress
plot_coolant_velocity = true
plot_hx_coefficient   = true
plot_coolant_temp     = true
plot_channel_width    = false
```

Plots that require biprop mode are silently skipped when `monoMode = 1`.

### Custom Propellants

Custom fuel blends, oxidisers, and monopropellants follow the rocketCEA card format. Each component is a sub-array of `[CEA_name, elemental_composition, wt_percent]`. Weight percentages must sum to 100.

**Example — 75% ethanol / 25% water fuel blend:**

```toml
newFuelBool  = 1
newFuelName  = "75_ETH"
newFuel      = [
    ["C2H5OH(L)", "C 2 H 6 O 1", "75", "Ethanol"],
    ["H2O(L)",    "H 2 O 1",     "25", "Water"  ],
]
```

**Example — gelled NTO oxidiser:**

```toml
newOxBool  = 1
newOxName  = "GelN2O4"
newOx      = [
    ["N2O4(L)", "N 2 O 4",  "96.5"],
    ["SiO2",    "Si 1 O 2", "3.5" ],
]
```

**Example — monoprop mode (90% H2O2):**

```toml
monoMode     = 1
newMonoBool  = 1
newMonoName  = "90_HTP"
newMono      = [
    ["H2O2(L)", "H 2 O 2", "90"],
    ["H2O(L)",  "H 2 O 1", "10"],
]
```

---

## Report Output

Running `python Engine_Class.py my_engine.toml ./output` produces:

```
./output/
    my_engine_report.docx     ← Word document (open in Word or Google Docs)
    report_data.json          ← intermediate data (useful for debugging)
    plot_engine_contour.png
    plot_mach_number.png
    plot_heat_flux.png
    ...                       ← one PNG per enabled plot flag
```

The `.docx` contains three sections:

**Section 1 — Inputs:** every parameter in a formatted table, followed by a complete TOML block you can copy directly into a new input file. The TOML is syntactically valid — booleans render as `true`/`false`, strings are quoted, and list-of-lists (`newFuel`, `newOx`) are written in proper TOML array syntax.

**Section 2 — Chamber Info & Dimensions:** all outputs from `printChamberInfo()` in a labelled table with units — area ratio, chamber pressure, temperatures, mass flow rates, propellant volumes, Isp, C*, geometry radii.

**Section 3 — Plots:** enabled plots embedded at 9 × 4 inches, two per page, on US Letter paper with 1-inch margins. Each page break is placed automatically.

---

## Using the Dictionaries

### fluid_dict

```python
from fluid_dict import fluid_dict, get_fluid_props, add_fluid_state, populate_from_coolprop

# exact lookup — (fluid_string, T_K, P_Pa)
props = fluid_dict["Ethanol[.5387]&Water[.4613]"][(353.15, 101325.0)]
rho   = props["density_kg_m3"]          # 820.3 kg/m³
k     = props["thermal_conductivity_W_mK"]  # 0.208 W/(m·K)
Pr    = props["prandtl_number"]          # 11.15

# Nearest-state lookup with tolerance fallback
props = get_fluid_props("Ethanol[.5387]&Water[.4613]", T_K=354.0, P_Pa=101325.0)

# Add a single entry at runtime
add_fluid_state("Ethanol[.5387]&Water[.4613]", T_K=400.0, P_Pa=2e6, props={
    "density_kg_m3": 780.0,
    "dynamic_viscosity_Pa_s": 4.1e-4,
    ...
})

# Batch-populate over a grid using live CoolProp calls
import numpy as np
populate_from_coolprop(
    "Ethanol[.5387]&Water[.4613]",
    T_range_K  = np.linspace(270, 500, 50),
    P_range_Pa = [101325, 1e6, 3e6],
)
```

The fluid name string is intentionally the same format as `Engine_Class.fuelComposition`, so you can use it as a dictionary key and pass it to CoolProp without reformatting.

### mat_dict

```python
from mat_dict import mat_dict, get_mat_props, add_mat_state

# exact lookup at a registered temperature
props = mat_dict["Inconel_718"][973.15]     # 700 °C
k     = props["thermal_conductivity_W_mK"]  # 18.7 W/(m·K)
E     = props["youngs_modulus_GPa"]         # 156.0 GPa
sy    = props["yield_strength_MPa"]         # 860.0 MPa

# Linear interpolation between registered temperatures
props = get_mat_props("Inconel_718", T_K=850.0)   # interpolates 773 K ↔ 973 K
k     = props["thermal_conductivity_W_mK"]          # 17.35 W/(m·K)

# Add a new temperature point at runtime
add_mat_state("Inconel_718", T_K=1173.15, props={
    "temperature_K":               1173.15,
    "density_kg_m3":               7990.0,
    "thermal_conductivity_W_mK":   22.1,
    "coef_thermal_expansion_1_K":  1.59e-5,
    "youngs_modulus_GPa":          118.0,
    "yield_strength_MPa":          310.0,
})
```

Both dictionaries are plain `.py` files — import them directly. No database, no file I/O at lookup time.

---

## Python API

All three constructor styles produce an identical object:

```python
from Engine_Class import Engine

# 1. Constructor arguments (original style)
eng = Engine(thrust=500, Pc=300, OF=2.1, height_of_optimization=0, burnTime=30)

# 2. Constructor + file (file overrides matching args)
eng = Engine(thrust=500, Pc=300, OF=2.1, input_file='my_engine.toml')

# 3. File only (recommended for reproducible runs)
eng = Engine.from_file('my_engine.toml')
```

Override report settings programmatically after construction:

```python
eng = Engine.from_file('my_engine.toml')
eng.report_output_dir = '/data/results/run_07'
eng.report_dpi        = 300
eng.plot_thermal_stress = False   # disable one plot at runtime
eng.generate_report()
```

Call individual analysis stages without generating a full report:

```python
eng = Engine.from_file('my_engine.toml')
eng.getCeaProperties()
eng.engineDimensions()

print(f"Throat diameter: {eng.Dt*1000:.2f} mm")
print(f"Isp:             {eng.Isp:.1f} s")
print(f"C*:              {eng.Cstar:.1f} m/s")

eng.propFlowRates()
eng.heatTransfer()

import numpy as np
print(f"Peak heat flux:  {np.max(eng.qw)/1e6:.2f} MW/m²")
print(f"Max wall temp:   {np.max(eng.Tw1):.0f} K")
```

Save nozzle geometry for CAD import:

```python
eng = Engine.from_file('my_engine.toml')
eng.saveNozGeo('./cad_export')     # writes nozzle_geo1.txt (x y z columns)
eng.saveChamberGeo('./cad_export') # writes chamber_geo1.txt
```

---