import numpy as np
# --------------------------
# Gross Primary Productivity (GPP) Model
# --------------------------
def gpp_model(PAR, RVI, Ts, GPP_max, k, alpha):
    """
    Calculate GPP using PAR (light), RVI (vegetation), and soil temperature Ts.
    """
    # Temperature function (FT)
    FT = np.where(Ts < -2, 0, np.where(Ts > 10, 1, (Ts + 2) / 12))  # linearly scaled from -2°C to 10°C
    return (GPP_max * PAR / (k + PAR)) * (RVI / (RVI + alpha)) * FT


# --------------------------
# Reco Model 2 (Karki et al., 2019)
# --------------------------
def reco_model_2(RVI, Ts, t1, a, b, T10=10, T0=-46):
    """
    Ecosystem respiration as a function of RVI and temperature.
    """
    exponent = b * (1 / (T10 - T0) - 1 / (Ts - T0))
    return t1 + (a * RVI) * np.exp(exponent)


# --------------------------
# Reco Model 3 (Rigney et al., 2018)
# --------------------------
def reco_model_3(WTD, Ts, t1, b, c, T10=10, T0=-46):
    """
    Ecosystem respiration using water table depth and temperature.
    """
    exponent = b * (1 / (T10 - T0) - 1 / (Ts - T0))
    return t1 * np.exp(exponent) + (WTD + c) ** 2


# --------------------------
# Reco Model 4 (This Study)
# --------------------------
def reco_model_4(RVI, WTD, Ts, WTD_max, t1, a, b, c, T10=10, T0=-46):
    """
    Enhanced Reco model including vegetation and water stress interaction.
    """
    exponent = b * (1 / (T10 - T0) - 1 / (Ts - T0))
    return t1 + (a * RVI) + ((WTD - WTD_max) * c) ** 2 * np.exp(exponent)


# CH4 emission model (mg CH₄/m²/h)
def ch4_emission(WTD, a=1.0, b=0.08):
    return a * np.exp(b * WTD)  # More CH4 when WTD is high (near 0 cm)