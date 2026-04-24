"""Net radiation computation for SEBAL.

Rn = (1 - alpha) * Rs_in + Rl_in - Rl_out

where Rl_out is derived from Stefan-Boltzmann with Landsat LST and
NDVI-based surface emissivity.
"""

import ee
from sebal.config import SEBAL_CONSTANTS


def compute_outgoing_longwave(
    lst_k: ee.Image,
    emissivity: ee.Image,
) -> ee.Image:
    """Stefan-Boltzmann outgoing longwave: R_l_out = eps * sigma * T_s^4."""
    sigma = SEBAL_CONSTANTS["SIGMA"]
    return (emissivity.multiply(sigma).multiply(lst_k.pow(4))
            .rename("Rl_out"))


def compute_net_radiation(
    albedo: ee.Image,
    rs_in: ee.Image,
    rl_in: ee.Image,
    lst_k: ee.Image,
    emissivity: ee.Image,
) -> ee.Image:
    """
    Compute Net Radiation Rn at the surface.

    Parameters
    ----------
    albedo : ee.Image
        Per-pixel broadband albedo.
    rs_in : ee.Image
        Incoming shortwave radiation (W/m^2) from ERA5-Land.
    rl_in : ee.Image
        Incoming longwave radiation (W/m^2) from ERA5-Land.
    lst_k : ee.Image
        Land surface temperature (K) from Landsat thermal.
    emissivity : ee.Image
        Surface emissivity (dimensionless).

    Returns
    -------
    ee.Image
        Single-band 'Rn' image in W/m^2.
    """
    rl_out = compute_outgoing_longwave(lst_k, emissivity)
    rn = (rs_in.multiply(ee.Image(1).subtract(albedo))
          .add(rl_in)
          .subtract(rl_out)
          .rename("Rn"))
    return rn


def compute_air_density(
    p_surf: ee.Image,
    t_air_k: ee.Image,
) -> ee.Image:
    """Ideal-gas air density: rho = P / (R_spec * T)."""
    r_spec = SEBAL_CONSTANTS["R_SPEC"]
    return p_surf.divide(t_air_k.multiply(r_spec)).rename("rho_air")
