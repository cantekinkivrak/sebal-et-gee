"""Daily ET from SEBAL via the Evaporative Fraction (EF) method.

Key SEBAL assumption: EF = lambda_ET / (Rn - G) is approximately constant
through daytime hours. Daily ET is computed by scaling EF with the 24-hour
net radiation.
"""

import ee
from sebal.config import SEBAL_CONSTANTS


def compute_latent_heat_of_vaporization(t_air_k: ee.Image) -> ee.Image:
    """Temperature-dependent latent heat of vaporization (J/kg)."""
    # lambda = (2.501 - 0.00236 * (T - 273.15)) * 1e6
    return (t_air_k.subtract(273.15).multiply(-0.00236).add(2.501)
            .multiply(1e6))


def compute_latent_heat_flux(
    rn_minus_g: ee.Image,
    H: ee.Image,
) -> ee.Image:
    """Latent heat flux lambda_ET = Rn - G - H, floored at 0."""
    return rn_minus_g.subtract(H).max(ee.Image(0)).rename("LE")


def compute_evaporative_fraction(
    LE: ee.Image,
    rn_minus_g: ee.Image,
) -> ee.Image:
    """EF = LE / (Rn - G), clamped to [0, 1]."""
    return LE.divide(rn_minus_g).clamp(0, 1).rename("EF")


def compute_daily_et(
    LE: ee.Image,
    rn_minus_g: ee.Image,
    rn_24h: ee.Image,
    t_air_k: ee.Image,
) -> tuple[ee.Image, ee.Image]:
    """
    Compute daily ET (mm/day) via the Evaporative Fraction method.

    ET_daily = EF * 86400 * Rn_24h / (lambda * rho_water)

    Parameters
    ----------
    LE : ee.Image
        Instantaneous latent heat flux at overpass time (W/m^2).
    rn_minus_g : ee.Image
        Instantaneous available energy at overpass time (W/m^2).
    rn_24h : ee.Image
        24-hour mean net radiation (W/m^2) from ERA5-Land.
    t_air_k : ee.Image
        Near-surface air temperature (K) for lambda.

    Returns
    -------
    ET_daily : ee.Image
        Daily ET (mm/day), floored at 0.
    EF : ee.Image
        Evaporative fraction for diagnostics.
    """
    EF = compute_evaporative_fraction(LE, rn_minus_g)
    lam = compute_latent_heat_of_vaporization(t_air_k)

    # ET_daily (mm/day) = EF * (86400 s) * Rn_24h (W/m^2) / lam (J/kg)
    # Note: 1 kg water = 1 mm over 1 m^2, so dividing by rho_water (1000 kg/m^3)
    # is implicit in the mm conversion.
    ET_daily = (EF.multiply(rn_24h).multiply(86400)
                .divide(lam)
                .max(ee.Image(0))
                .rename("ET_daily_mm"))
    return ET_daily, EF
