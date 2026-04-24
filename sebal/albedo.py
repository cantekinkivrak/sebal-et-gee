"""Broadband surface albedo from Landsat 8/9 OLI using Tasumi (2008) coefficients.

Reference:
    Tasumi, M., Allen, R. G., & Trezza, R. (2008). At-surface reflectance
    and albedo from satellite for operational calculation of land surface
    energy balance. J. Hydrol. Eng., 13(2), 51-63.
"""

import ee
from sebal.config import DEFAULT_PARAMS


def compute_albedo_tasumi(mosaic: ee.Image) -> ee.Image:
    """
    Compute broadband albedo from atmospherically-corrected reflectance.

    alpha = 0.300*B2 + 0.277*B3 + 0.233*B4 + 0.143*B5 + 0.036*B6 + 0.012*B7

    Parameters
    ----------
    mosaic : ee.Image
        Landsat mosaic with scale factors applied (surface reflectance).

    Returns
    -------
    ee.Image
        Single-band 'albedo' image clamped to [0, 1].
    """
    coeffs = DEFAULT_PARAMS["ALBEDO_COEFFS"]
    terms = [mosaic.select(b).multiply(c) for b, c in coeffs.items()]
    albedo = terms[0]
    for t in terms[1:]:
        albedo = albedo.add(t)
    return albedo.rename("albedo").clamp(0.0, 1.0)


def compute_ndvi(mosaic: ee.Image) -> ee.Image:
    """Compute NDVI from Landsat NIR (B5) and Red (B4) bands."""
    return mosaic.normalizedDifference(["SR_B5", "SR_B4"]).rename("NDVI")


def compute_emissivity(ndvi: ee.Image) -> ee.Image:
    """Simple NDVI-based surface emissivity parameterization, clamped to [0.96, 0.99]."""
    return (ndvi.multiply(0.01).add(0.985)
            .clamp(0.96, 0.99).rename("emissivity"))
