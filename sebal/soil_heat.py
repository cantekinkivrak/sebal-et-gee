"""Soil heat flux G via the Bastiaanssen (1998) empirical formulation.

G / Rn = (T_s - 273.15) / alpha * (0.0038*alpha + 0.0074*alpha^2)
        * (1 - 0.98 * NDVI^4)

For water bodies (NDVI < 0), G/Rn is fixed at 0.5.

Reference:
    Bastiaanssen, W. G. M. (1998). Remote sensing in water resources
    management: The state of the art. J. Hydrol., 212-213, 198-212.
"""

import ee
from sebal.config import DEFAULT_PARAMS


def compute_soil_heat_flux(
    rn: ee.Image,
    lst_k: ee.Image,
    albedo: ee.Image,
    ndvi: ee.Image,
) -> tuple[ee.Image, ee.Image]:
    """
    Compute soil heat flux G and the G/Rn ratio.

    Returns
    -------
    G : ee.Image
        Soil heat flux (W/m^2).
    g_ratio : ee.Image
        G/Rn ratio, clamped to physical bounds.
    """
    water_ratio = DEFAULT_PARAMS["WATER_G_RATIO"]
    ratio_min = DEFAULT_PARAMS["G_RATIO_MIN"]
    ratio_max = DEFAULT_PARAMS["G_RATIO_MAX"]

    t_c = lst_k.subtract(273.15)
    term1 = t_c.divide(albedo)
    term2 = albedo.multiply(0.0038).add(albedo.pow(2).multiply(0.0074))
    term3 = ndvi.pow(4).multiply(0.98).multiply(-1).add(1)

    g_ratio = term1.multiply(term2).multiply(term3).rename("g_ratio")

    # Water override
    water = ndvi.lt(0)
    g_ratio = g_ratio.where(water, water_ratio)
    g_ratio = g_ratio.clamp(ratio_min, ratio_max)

    G = rn.multiply(g_ratio).rename("G")
    return G, g_ratio
