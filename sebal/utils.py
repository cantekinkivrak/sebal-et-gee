"""Utility helpers for GEE initialization and geometry/data handling."""

from pathlib import Path
import ee
import geopandas as gpd
import pandas as pd


def init_ee(project: str) -> None:
    """Initialize Earth Engine with a given cloud project; authenticate if needed."""
    try:
        ee.Initialize(project=project)
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=project)


def load_aoi(geojson_path: str | Path) -> tuple[gpd.GeoDataFrame, ee.FeatureCollection, ee.Geometry]:
    """Load a GeoJSON AOI and convert to Earth Engine geometry."""
    gdf = gpd.read_file(str(geojson_path))
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    elif gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs("EPSG:4326")

    # Minimal geopandas -> ee conversion without geemap dependency
    features = []
    for _, row in gdf.iterrows():
        geom = row.geometry
        if geom is None:
            continue
        features.append(ee.Feature(ee.Geometry(geom.__geo_interface__),
                                   {k: v for k, v in row.items() if k != "geometry"}))
    aoi_ee = ee.FeatureCollection(features)
    return gdf, aoi_ee, aoi_ee.geometry()


def stations_to_feature_collection(
    stations_df: pd.DataFrame,
    lat_col: str = "lat",
    lon_col: str = "lon",
    name_col: str = "name",
) -> ee.FeatureCollection:
    """Convert a pandas DataFrame of station coordinates to an ee.FeatureCollection."""
    features = []
    for _, row in stations_df.iterrows():
        props = {k: v for k, v in row.items()
                 if k not in (lat_col, lon_col) and pd.notna(v)}
        pt = ee.Geometry.Point([float(row[lon_col]), float(row[lat_col])])
        features.append(ee.Feature(pt, props))
    return ee.FeatureCollection(features)


def extract_at_stations(
    image: ee.Image,
    stations_fc: ee.FeatureCollection,
    scale: int = 30,
    reducer: ee.Reducer = None,
) -> pd.DataFrame:
    """Extract image values at station locations; return a pandas DataFrame."""
    if reducer is None:
        reducer = ee.Reducer.mean()
    reduced = image.reduceRegions(
        collection=stations_fc,
        reducer=reducer,
        scale=scale,
    ).getInfo()

    rows = []
    for feat in reduced["features"]:
        row = feat["properties"].copy()
        rows.append(row)
    return pd.DataFrame(rows)


def export_to_drive(
    image: ee.Image,
    description: str,
    aoi_geom: ee.Geometry,
    folder: str = "sebal_et_gee_exports",
    scale: int = 30,
    max_pixels: int = int(1e10),
) -> ee.batch.Task:
    """Start a Drive export task and return the task handle."""
    task = ee.batch.Export.image.toDrive(
        image=image.toFloat(),
        description=description,
        folder=folder,
        region=aoi_geom,
        scale=scale,
        maxPixels=max_pixels,
        fileFormat="GeoTIFF",
    )
    task.start()
    return task
