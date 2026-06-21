# fluid_equations.py

def get_CLF3_props(T_K: float, P_Pa: float) -> dict:
    """Chlorine trifluoride — Gmelin Handbook correlations. Valid 200–400 K."""
    rho = 1950.0 - 0.85 * (T_K - 200.0)
    mu  = 1.2e-3 * (200.0 / T_K) ** 1.8
    Cp  = 920.0
    k   = 0.195 - 0.0002 * (T_K - 200.0)
    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            (mu * Cp / k),
    }

def get_Nitric_Acid_props(T_K: float, P_Pa: float) -> dict:
    """HNO3 — Pulled from https://webbook.nist.gov/cgi/cbook.cgi?ID=C7697372&Mask=1&Type=JANAFG&Table=on"""
    if (T_K <=1200) and (T_K >= 298):
        A=19.63229,     B=153.9599,     C=-115.8378,    D=32.87955
        E=-0.249114,    F=-146.8818,    G=247.7049,     H=-134.3060
        
    if T_K > 1200:
        A=97.45959,     B=5.429577,     C=-1.029688,    D=0.067950
        E=-12.29314,    F=-192.4912,    G=343.8051,     H=-134.3060
        

def get_Jet_A_props(T_K: float, P_Pa: float) -> dict:
    """Jet-A fuel approximated as C12H23."""
    ...