# fluid_equations.py

def get_CLF3_props(T_K: float, P_Pa: float) -> dict:
    """
    Chlorine trifluoride (ClF3). Gmelin Handbook of Inorganic Chemistry,
    8th ed. Valid range: 200-400 K (freeze 197 K, bp 285 K).
    Reference point: rho(298 K) = 1825 kg/m3 (NIST), mu near bp ~4.5e-4 Pa·s.
    """
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


def get_Gasoline_props(T_K: float, P_Pa: float) -> dict:
    """
    Automotive gasoline, approximated as a C7-C8 hydrocarbon blend
    (iso-octane proxy). Typical pump-gasoline reference values at 288 K
    (15 C): rho ~ 745 kg/m3, mu ~ 6.0e-4 Pa·s, Cp ~ 2020 J/kgK,
    k ~ 0.12 W/mK. Valid range: 211-450 K (approx freeze to safe upper
    bound well below auto-ignition).
    Source: SAE J1297 typical gasoline property tables; CRC Handbook
    hydrocarbon liquid correlations (iso-octane as surrogate).
    """
    T_ref = 288.15
    rho = 745.0 - 0.85 * (T_K - T_ref)
    mu  = 6.0e-4 * (T_ref / T_K) ** 2.0
    Cp  = 2020.0 + 3.0 * (T_K - T_ref)
    k   = 0.120 - 0.00010 * (T_K - T_ref)
    Pr  = (mu * Cp / k) if k else None
    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            Pr,
    }


def get_IPA_props(T_K: float, P_Pa: float) -> dict:
    """
    Pure isopropanol (2-propanol, C3H8O). NIST WebBook / CRC Handbook.
    Reference point: T_ref = 298.15 K -> rho = 785.1 kg/m3,
    mu = 2.043e-3 Pa·s, Cp = 2680 J/kgK, k = 0.135 W/mK.
    Valid range: 186-490 K (freeze 185 K, bp 356 K, well below crit 508 K).
    """
    T_ref = 298.15
    rho = 785.1 - 0.95 * (T_K - T_ref)
    mu  = 2.043e-3 * (T_ref / T_K) ** 2.3
    Cp  = 2680.0 + 4.2 * (T_K - T_ref)
    k   = 0.135 - 0.00018 * (T_K - T_ref)
    Pr  = (mu * Cp / k) if k else None
    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            Pr,
    }


def get_IPA90_props(T_K: float, P_Pa: float) -> dict:
    """
    90 wt% IPA / 10 wt% water blend. CRC Handbook aqueous-alcohol tables.
    Reference point: T_ref = 298.15 K -> rho = 813.0 kg/m3,
    mu = 2.55e-3 Pa·s (water raises viscosity above pure IPA near the
    azeotrope), Cp = 2640 J/kgK, k = 0.150 W/mK (water raises k).
    Valid range: 186-490 K.
    """
    T_ref = 298.15
    rho = 813.0 - 0.92 * (T_K - T_ref)
    mu  = 2.55e-3 * (T_ref / T_K) ** 2.2
    Cp  = 2640.0 + 4.0 * (T_K - T_ref)
    k   = 0.150 - 0.00017 * (T_K - T_ref)
    Pr  = (mu * Cp / k) if k else None
    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            Pr,
    }


def get_IPA80_props(T_K: float, P_Pa: float) -> dict:
    """
    80 wt% IPA / 20 wt% water blend. CRC Handbook aqueous-alcohol tables.
    Reference point: T_ref = 298.15 K -> rho = 850.0 kg/m3,
    mu = 3.05e-3 Pa·s, Cp = 2590 J/kgK, k = 0.165 W/mK.
    Valid range: 186-490 K.
    """
    T_ref = 298.15
    rho = 850.0 - 0.90 * (T_K - T_ref)
    mu  = 3.05e-3 * (T_ref / T_K) ** 2.1
    Cp  = 2590.0 + 3.8 * (T_K - T_ref)
    k   = 0.165 - 0.00016 * (T_K - T_ref)
    Pr  = (mu * Cp / k) if k else None
    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            Pr,
    }


def get_IPA70_props(T_K: float, P_Pa: float) -> dict:
    """
    70 wt% IPA / 30 wt% water blend. CRC Handbook aqueous-alcohol tables.
    Reference point: T_ref = 298.15 K -> rho = 876.0 kg/m3,
    mu = 3.45e-3 Pa·s, Cp = 2540 J/kgK, k = 0.180 W/mK.
    Valid range: 186-490 K.
    """
    T_ref = 298.15
    rho = 876.0 - 0.88 * (T_K - T_ref)
    mu  = 3.45e-3 * (T_ref / T_K) ** 2.0
    Cp  = 2540.0 + 3.6 * (T_K - T_ref)
    k   = 0.180 - 0.00015 * (T_K - T_ref)
    Pr  = (mu * Cp / k) if k else None
    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            Pr,
    }


def get_JetA_props(T_K: float, P_Pa: float) -> dict:
    """
    Jet-A aviation kerosene. ASTM D1655 typical property ranges.
    Reference point: T_ref = 288.15 K (15 C) -> rho = 808.0 kg/m3
    (mid-spec, range 775-840), mu = 1.8e-3 Pa·s, Cp = 2000 J/kgK,
    k = 0.110 W/mK. Valid range: 234-600 K (freeze ~233 K through
    superheated film-cooling conditions, well below the 430-570 K
    distillation range where composition begins to shift).
    Source: ASTM D1655 / CRC handbook kerosene-range hydrocarbon data.
    """
    T_ref = 288.15
    rho = 808.0 - 0.70 * (T_K - T_ref)
    mu  = 1.8e-3 * (T_ref / T_K) ** 2.4
    Cp  = 2000.0 + 3.3 * (T_K - T_ref)
    k   = 0.110 - 0.00009 * (T_K - T_ref)
    Pr  = (mu * Cp / k) if k else None
    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            Pr,
    }


def get_M20_props(T_K: float, P_Pa: float) -> dict:
    """
    M20: 20 wt% MMH / 80 wt% hydrazine (N2H4) hypergolic fuel blend.
    Properties computed by linear (ideal) weight-fraction mixing of the
    pure-component reference values at 298.15 K:
        Hydrazine(298 K): rho=1004 kg/m3, mu=0.913 cP, Cp=3084 J/kgK,
                          k=0.340 W/mK
        MMH(298 K):       rho=874  kg/m3, mu=0.850 cP, Cp=2980 J/kgK,
                          k=0.245 W/mK
    Mixing rule: prop_blend = 0.20*prop_MMH + 0.80*prop_N2H4
    Each pure-component property is then scaled with the same temperature
    dependence used for the pure hydrazine/MMH correlations (density linear,
    viscosity Andrade, Cp weak linear, k linear decreasing).
    Valid range: 169-430 K (bounded by the lower freeze point of the blend
    and well below either component's decomposition onset).
    Source: Rocket Propulsion Elements (Sutton) component tables; ideal
    mixing approximation — no published M20 blend data exists publicly.
    """
    T_ref = 298.15

    # Pure hydrazine (80 wt%)
    rho_N2H4 = 1004.0 - 0.85 * (T_K - T_ref)
    mu_N2H4  = 0.913e-3 * (T_ref / T_K) ** 1.9
    Cp_N2H4  = 3084.0 + 2.0 * (T_K - T_ref)
    k_N2H4   = 0.340 - 0.00035 * (T_K - T_ref)

    # Pure MMH (20 wt%)
    rho_MMH = 874.0 - 0.80 * (T_K - T_ref)
    mu_MMH  = 0.850e-3 * (T_ref / T_K) ** 1.9
    Cp_MMH  = 2980.0 + 2.4 * (T_K - T_ref)
    k_MMH   = 0.245 - 0.00028 * (T_K - T_ref)

    # Ideal weight-fraction mixing
    rho = 0.20 * rho_MMH + 0.80 * rho_N2H4
    mu  = 0.20 * mu_MMH  + 0.80 * mu_N2H4
    Cp  = 0.20 * Cp_MMH  + 0.80 * Cp_N2H4
    k   = 0.20 * k_MMH   + 0.80 * k_N2H4
    Pr  = (mu * Cp / k) if k else None

    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            Pr,
    }


def get_MON3_props(T_K: float, P_Pa: float) -> dict:
    """
    MON-3: nitrogen tetroxide (NTO) with 3 wt% dissolved nitric oxide (NO).
    Base properties from pure NTO at 298.15 K: rho=1443 kg/m3,
    mu=0.42 cP, Cp=1480 J/kgK, k=0.140 W/mK (CPIA/Rocket Propulsion
    Elements). NO addition lowers density slightly and depresses the
    freeze/boiling points; density correction approximated as
    -3.0 kg/m3 per wt% NO (linear, small-perturbation approximation).
    Valid range: 259-420 K (3 wt% NO depresses bp to ~289 K vs 294 K
    for pure NTO; upper bound stays below decomposition onset).
    Source: CPIA Publication 194 (Mixed Oxides of Nitrogen);
    Rocket Propulsion Elements (Sutton), NTO property tables.
    """
    T_ref = 298.15
    NO_wt_pct = 3.0

    rho_NTO_ref = 1443.0 - 3.0 * NO_wt_pct      # 1434.0 kg/m3 at T_ref
    rho = rho_NTO_ref - 1.95 * (T_K - T_ref)
    mu  = 0.42e-3 * (T_ref / T_K) ** 1.6
    Cp  = 1480.0 + 1.5 * (T_K - T_ref)
    k   = 0.140 - 0.00022 * (T_K - T_ref)
    Pr  = (mu * Cp / k) if k else None
    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            Pr,
    }


def get_MON15_props(T_K: float, P_Pa: float) -> dict:
    """
    MON-15: NTO with 15 wt% dissolved NO. Same base NTO properties and
    mixing approach as MON-3, with the larger NO fraction depressing
    density further and lowering bp to ~278 K.
    Valid range: 249-410 K.
    Source: CPIA Publication 194 (Mixed Oxides of Nitrogen).
    """
    T_ref = 298.15
    NO_wt_pct = 15.0

    rho_NTO_ref = 1443.0 - 3.0 * NO_wt_pct      # 1398.0 kg/m3 at T_ref
    rho = rho_NTO_ref - 1.90 * (T_K - T_ref)
    mu  = 0.40e-3 * (T_ref / T_K) ** 1.6
    Cp  = 1500.0 + 1.5 * (T_K - T_ref)
    k   = 0.137 - 0.00022 * (T_K - T_ref)
    Pr  = (mu * Cp / k) if k else None
    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            Pr,
    }


def get_MON25_props(T_K: float, P_Pa: float) -> dict:
    """
    MON-25: NTO with 25 wt% dissolved NO. Same base NTO properties and
    mixing approach as MON-3/MON-15, with the largest NO fraction in this
    series depressing density and bp the most (bp ~268 K).
    Valid range: 239-400 K.
    Source: CPIA Publication 194 (Mixed Oxides of Nitrogen).
    """
    T_ref = 298.15
    NO_wt_pct = 25.0

    rho_NTO_ref = 1443.0 - 3.0 * NO_wt_pct      # 1368.0 kg/m3 at T_ref
    rho = rho_NTO_ref - 1.85 * (T_K - T_ref)
    mu  = 0.38e-3 * (T_ref / T_K) ** 1.6
    Cp  = 1520.0 + 1.5 * (T_K - T_ref)
    k   = 0.133 - 0.00022 * (T_K - T_ref)
    Pr  = (mu * Cp / k) if k else None
    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            Pr,
    }


def get_N2F4_props(T_K: float, P_Pa: float) -> dict:
    """
    Dinitrogen tetrafluoride (N2F4). Exotic high-energy oxidiser, same
    chemical family as N2O4 but lighter and colder (bp 200 K vs 294 K).
    Reference point: T_ref = 190 K (near bp) -> rho ~ 1500 kg/m3
    (estimated from liquid molar volume, no direct experimental density
    series published), mu ~ 3.0e-4 Pa·s, Cp ~ 1000 J/kgK, k ~ 0.130 W/mK.
    Valid range: 110-305 K (freeze ~109 K, bp 200 K, supercrit 309 K).
    Source: estimated from N2O4 family correlations scaled by molar mass
    ratio (N2F4 = 104.0 g/mol vs N2O4 = 92.0 g/mol) — no direct literature
    transport-property series exists publicly; treat as order-of-magnitude.
    """
    T_ref = 190.0
    rho = 1500.0 - 2.10 * (T_K - T_ref)
    mu  = 3.0e-4 * (T_ref / T_K) ** 1.7
    Cp  = 1000.0 + 1.8 * (T_K - T_ref)
    k   = 0.130 - 0.00030 * (T_K - T_ref)
    Pr  = (mu * Cp / k) if k else None
    return {
        "density_kg_m3":             rho,
        "dynamic_viscosity_Pa_s":    mu,
        "specific_heat_J_kgK":       Cp,
        "thermal_conductivity_W_mK": k,
        "prandtl_number":            Pr,
    }