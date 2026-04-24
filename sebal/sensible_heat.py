"""Sensible heat flux via the CIMEC iterative SEBAL scheme.

Steps:
1. Compute roughness length z0m from NDVI.
2. Compute friction velocity u* and neutral aerodynamic resistance r_ah.
3. Select cold and hot anchor pixels (extreme ET endpoints).
4. Iteratively solve H with Monin-Obukhov stability correction,
   enforcing H_cold = 0 and H_hot = Rn_hot - G_hot at each step.
"""

import math
import ee
from sebal.config import SEBAL_CONSTANTS, DEFAULT_PARAMS


# --- Roughness and aerodynamic resistance ---

def compute_roughness(ndvi: ee.Image) -> ee.Image:
    """Momentum roughness length z0m as a simple NDVI-based proxy."""
    coef = DEFAULT_PARAMS["Z0M_NDVI_COEF"]
    floor = DEFAULT_PARAMS["Z0M_FLOOR"]
    return ndvi.multiply(coef).max(ee.Image(floor)).rename("z0m")


def compute_neutral_aerodynamics(
    wind_10m: ee.Image,
    z0m: ee.Image,
) -> tuple[ee.Image, ee.Image, ee.Image]:
    """Neutral-stability friction velocity, 200m wind, and r_ah between z1 and z2."""
    k = SEBAL_CONSTANTS["K_VK"]
    z_blend = DEFAULT_PARAMS["Z_BLEND"]
    z_wind = DEFAULT_PARAMS["Z_WIND"]
    z1 = DEFAULT_PARAMS["Z1"]
    z2 = DEFAULT_PARAMS["Z2"]

    u_star_10 = wind_10m.multiply(k).divide(
        z0m.multiply(z_wind).log().subtract(z0m.log())
    )
    u_200 = u_star_10.divide(k).multiply(
        z0m.multiply(z_blend).log().subtract(z0m.log())
    )
    u_star = (u_200.multiply(k)
              .divide(ee.Image(z_blend).divide(z0m).log())
              .rename("u_star"))
    r_ah = (ee.Image(z2 / z1).log()
            .divide(u_star.multiply(k))
            .rename("r_ah"))
    return u_star, u_200, r_ah


# --- Anchor pixel selection ---

def select_anchor_pixels(
    ndvi: ee.Image,
    lst_k: ee.Image,
    albedo: ee.Image,
    rn: ee.Image,
    G: ee.Image,
    aoi_geom: ee.Geometry,
) -> tuple[dict, dict, ee.Image, ee.Image]:
    """
    Select hot and cold anchor pixels using NDVI/LST percentile thresholds
    over an agricultural mask.

    Returns
    -------
    cold_stats, hot_stats : dict
        Median values of NDVI, LST_K, albedo, Rn, G at the anchor pools.
    cold_mask, hot_mask : ee.Image
        Binary masks for the anchor candidate pixels.
    """
    ag_mask = (ndvi.gt(DEFAULT_PARAMS["AG_NDVI_MIN"])
               .And(albedo.lt(DEFAULT_PARAMS["AG_ALBEDO_MAX"]))
               .And(albedo.gt(DEFAULT_PARAMS["AG_ALBEDO_MIN"])))

    ag_stack = (ee.Image.cat([ndvi, lst_k.rename("LST_K"), albedo, rn, G])
                .updateMask(ag_mask))

    pcts = ag_stack.reduceRegion(
        reducer=ee.Reducer.percentile(DEFAULT_PARAMS["ANCHOR_PERCENTILES"]),
        geometry=aoi_geom,
        scale=90,
        maxPixels=int(1e9),
        bestEffort=True,
    ).getInfo()

    cold_mask = (ag_stack.select("NDVI").gt(pcts["NDVI_p95"])
                 .And(ag_stack.select("LST_K").lt(pcts["LST_K_p20"])))
    hot_mask = (ag_stack.select("NDVI").gt(pcts["NDVI_p5"])
                .And(ag_stack.select("NDVI").lt(pcts["NDVI_p20"]))
                .And(ag_stack.select("LST_K").gt(pcts["LST_K_p95"])))

    def _median_stats(img: ee.Image, mask: ee.Image) -> dict:
        return img.updateMask(mask).reduceRegion(
            reducer=ee.Reducer.median(),
            geometry=aoi_geom,
            scale=90,
            maxPixels=int(1e9),
            bestEffort=True,
        ).getInfo()

    cold_stats = _median_stats(ag_stack, cold_mask)
    hot_stats = _median_stats(ag_stack, hot_mask)
    return cold_stats, hot_stats, cold_mask, hot_mask


# --- Monin-Obukhov stability functions (unstable conditions) ---

def _psi_m_unstable(z: float, L: ee.Image) -> ee.Image:
    x = ee.Image(1).subtract(L.pow(-1).multiply(z).multiply(16)).pow(0.25)
    return (x.add(1).divide(2).log().multiply(2)
            .add(x.pow(2).add(1).divide(2).log())
            .subtract(x.atan().multiply(2))
            .add(ee.Image(math.pi / 2)))


def _psi_h_unstable(z: float, L: ee.Image) -> ee.Image:
    x = ee.Image(1).subtract(L.pow(-1).multiply(z).multiply(16)).pow(0.25)
    return x.pow(2).add(1).divide(2).log().multiply(2)


# --- CIMEC iterative H solution ---

def solve_sensible_heat_cimec(
    lst_k: ee.Image,
    rho_air: ee.Image,
    u_star_init: ee.Image,
    u_200: ee.Image,
    z0m: ee.Image,
    r_ah_init: ee.Image,
    rn_minus_g: ee.Image,
    cold_mask: ee.Image,
    hot_mask: ee.Image,
    cold_stats: dict,
    hot_stats: dict,
    aoi_geom: ee.Geometry,
    n_iter: int = None,
    verbose: bool = True,
) -> tuple[ee.Image, ee.Image, dict]:
    """
    Iterative CIMEC solution for sensible heat flux H.

    At each iteration:
    - H_hot is fixed at Rn_hot - G_hot (lambda_ET ~= 0 assumption)
    - H_cold is fixed at 0 (lambda_ET = Rn - G assumption)
    - dT linear fit re-derived from current anchor r_ah
    - Full-image H computed, Monin-Obukhov length updated
    - r_ah and u_star corrected for stability

    Returns
    -------
    H : ee.Image
        Final sensible heat flux (W/m^2), bounded by [0, Rn-G].
    r_ah : ee.Image
        Final stability-corrected aerodynamic resistance.
    fit : dict
        Final linear fit parameters {a, b, dT_hot, r_ah_hot}.
    """
    cp = SEBAL_CONSTANTS["CP"]
    k = SEBAL_CONSTANTS["K_VK"]
    g_grav = SEBAL_CONSTANTS["G_GRAV"]
    z_blend = DEFAULT_PARAMS["Z_BLEND"]
    z1 = DEFAULT_PARAMS["Z1"]
    z2 = DEFAULT_PARAMS["Z2"]
    h_floor = DEFAULT_PARAMS["H_FLOOR"]
    r_ah_floor = DEFAULT_PARAMS["R_AH_FLOOR"]

    if n_iter is None:
        n_iter = DEFAULT_PARAMS["N_ITER_H"]

    Rn_hot_val = ee.Number(hot_stats["Rn"])
    G_hot_val = ee.Number(hot_stats["G"])
    H_hot_fixed = Rn_hot_val.subtract(G_hot_val)

    Ts_cold = ee.Number(cold_stats["LST_K"])
    Ts_hot = ee.Number(hot_stats["LST_K"])

    def _anchor_median(img: ee.Image, mask: ee.Image) -> ee.Number:
        return ee.Number(img.updateMask(mask).reduceRegion(
            reducer=ee.Reducer.median(),
            geometry=aoi_geom,
            scale=90,
            maxPixels=int(1e9),
            bestEffort=True,
        ).values().get(0))

    r_ah = r_ah_init
    u_star = u_star_init
    dT = None
    b_slope = None
    a_intercept = None

    if verbose:
        print(f"H_hot_fixed = {H_hot_fixed.getInfo():.1f} W/m^2")

    for i in range(n_iter):
        r_ah_hot = _anchor_median(r_ah, hot_mask)
        rho_hot = _anchor_median(rho_air, hot_mask)

        dT_hot = H_hot_fixed.multiply(r_ah_hot).divide(rho_hot.multiply(cp))
        dT_cold = ee.Number(0.0)

        b_slope = dT_hot.subtract(dT_cold).divide(Ts_hot.subtract(Ts_cold))
        a_intercept = dT_cold.subtract(b_slope.multiply(Ts_cold))

        dT = lst_k.multiply(b_slope).add(a_intercept).rename("dT")

        H = (rho_air.multiply(cp).multiply(dT).divide(r_ah)
             .max(ee.Image(h_floor)).rename("H"))

        L = (rho_air.multiply(cp).multiply(lst_k).multiply(u_star.pow(3))
             .divide(H.multiply(k).multiply(g_grav)).multiply(-1)
             .rename("L"))

        unstable = L.lt(0)
        psi_m_200 = _psi_m_unstable(z_blend, L).updateMask(unstable).unmask(0)
        psi_h_2 = _psi_h_unstable(z2, L).updateMask(unstable).unmask(0)
        psi_h_01 = _psi_h_unstable(z1, L).updateMask(unstable).unmask(0)

        u_star = (u_200.multiply(k)
                  .divide(ee.Image(z_blend).divide(z0m).log().subtract(psi_m_200)))

        r_ah = (ee.Image(z2 / z1).log()
                .subtract(psi_h_2).add(psi_h_01)
                .divide(u_star.multiply(k))
                .max(ee.Image(r_ah_floor))
                .rename("r_ah"))

        if verbose:
            print(f"  Iter {i+1}/{n_iter}: "
                  f"b_slope={b_slope.getInfo():.4f}, "
                  f"dT_hot={dT_hot.getInfo():.2f} K, "
                  f"r_ah_hot={r_ah_hot.getInfo():.1f} s/m")

    H = (rho_air.multiply(cp).multiply(dT).divide(r_ah)
         .max(ee.Image(0))
         .min(rn_minus_g)
         .rename("H"))

    fit = {
        "a_intercept": a_intercept.getInfo() if a_intercept is not None else None,
        "b_slope": b_slope.getInfo() if b_slope is not None else None,
    }
    return H, r_ah, fit
