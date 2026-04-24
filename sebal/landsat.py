"""Landsat 8/9 Collection 2 Level 2 ingestion, cloud masking, and mosaic building."""

import ee
from sebal.config import DEFAULT_PARAMS


def apply_scale_factors(image: ee.Image) -> ee.Image:
    """Apply Landsat C2 L2 scale factors and add LST_K, LST_C bands."""
    optical = image.select("SR_B.").multiply(0.0000275).add(-0.2)
    thermal_k = image.select("ST_B10").multiply(0.00341802).add(149.0).rename("LST_K")
    thermal_c = thermal_k.subtract(273.15).rename("LST_C")
    return (image.addBands(optical, None, True)
                 .addBands(thermal_k)
                 .addBands(thermal_c))


def mask_clouds(image: ee.Image) -> ee.Image:
    """Mask clouds and shadows using QA_PIXEL bitmask."""
    qa = image.select("QA_PIXEL")
    cloud = qa.bitwiseAnd(1 << 3).neq(0)
    shadow = qa.bitwiseAnd(1 << 4).neq(0)
    return image.updateMask(cloud.Or(shadow).Not())


def build_landsat_mosaic(
    aoi_geom: ee.Geometry,
    anchor_date: str,
    window_days: int = 14,
    max_cloud: float = None,
) -> tuple[ee.Image, ee.ImageCollection, float]:
    """
    Build a Landsat 8+9 mosaic over the AOI in a temporal window.

    Parameters
    ----------
    aoi_geom : ee.Geometry
        Region of interest.
    anchor_date : str
        Center date of the window in 'YYYY-MM-DD' format.
    window_days : int
        Total window length in days (centered on anchor).
    max_cloud : float
        Maximum scene-level CLOUD_COVER filter (default from config).

    Returns
    -------
    mosaic : ee.Image
        Cloud-masked, scale-factor-applied Landsat mosaic clipped to AOI.
    merged : ee.ImageCollection
        Underlying L8+L9 collection used for the mosaic.
    coverage_pct : float
        Fraction of AOI covered by the mosaic (0-100).
    """
    if max_cloud is None:
        max_cloud = DEFAULT_PARAMS["LANDSAT_CLOUD_COVER_MAX"]

    bands = DEFAULT_PARAMS["LANDSAT_BANDS"]
    anchor = ee.Date(anchor_date)
    start = anchor.advance(-window_days // 2, "day")
    end = anchor.advance(window_days // 2, "day")

    collections = []
    for coll_id in DEFAULT_PARAMS["LANDSAT_COLLECTIONS"]:
        c = (ee.ImageCollection(coll_id)
             .filterBounds(aoi_geom)
             .filterDate(start, end)
             .filter(ee.Filter.lt("CLOUD_COVER", max_cloud))
             .select(bands))
        collections.append(c)

    merged = collections[0]
    for c in collections[1:]:
        merged = merged.merge(c)

    merged = merged.map(mask_clouds).map(apply_scale_factors).sort("CLOUD_COVER")
    mosaic = merged.mosaic().clip(aoi_geom)

    # Compute AOI coverage as the fraction of non-null LST_K pixels
    coverage = mosaic.select("LST_K").mask().reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=aoi_geom,
        scale=300,
        maxPixels=int(1e9),
    ).getInfo().get("LST_K", 0)
    coverage_pct = (coverage or 0) * 100.0

    return mosaic, merged, coverage_pct


def find_best_anchor(
    aoi_geom: ee.Geometry,
    candidate_dates: list[str],
    window_days: int = 14,
    max_cloud: float = None,
) -> tuple[str, float]:
    """Search multiple anchor dates, return the one with maximum AOI coverage."""
    best = (None, -1.0)
    for d in candidate_dates:
        _, _, cov = build_landsat_mosaic(aoi_geom, d, window_days, max_cloud)
        if cov > best[1]:
            best = (d, cov)
    return best
