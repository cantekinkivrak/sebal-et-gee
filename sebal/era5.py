"""ERA5-Land hourly reanalysis ingestion for SEBAL meteorology inputs."""

import ee
from sebal.config import DEFAULT_PARAMS


def load_era5_meteorology(
    aoi_geom: ee.Geometry,
    start_date: ee.Date,
    end_date: ee.Date,
    reference_hour_utc: int = None,
) -> ee.Image:
    """
    Build a window-mean ERA5-Land image at the given reference hour.

    This aggregates one hourly slice per day in the window and averages them,
    which is a standard SEBAL simplification for multi-scene Landsat mosaics.

    Parameters
    ----------
    aoi_geom : ee.Geometry
    start_date, end_date : ee.Date
    reference_hour_utc : int
        Hour of day (UTC) to sample; default from config (10 UTC for Thrace).

    Returns
    -------
    ee.Image
        Multi-band image containing T_air_K, T_air_C, wind_10m, Rs_in, Rl_in, P_surf.
    """
    if reference_hour_utc is None:
        reference_hour_utc = DEFAULT_PARAMS["ERA5_REFERENCE_HOUR_UTC"]

    era5_vars = DEFAULT_PARAMS["ERA5_VARS"]
    coll_id = DEFAULT_PARAMS["ERA5_COLLECTION"]

    n_days = int(end_date.difference(start_date, "day").getInfo())
    daily_snaps = []
    for i in range(n_days):
        d = start_date.advance(i, "day")
        h0 = d.update(hour=reference_hour_utc, minute=0, second=0)
        snap = (ee.ImageCollection(coll_id)
                .filterDate(h0, h0.advance(1, "hour"))
                .select(era5_vars)
                .first())
        daily_snaps.append(snap)

    era5_mean = ee.ImageCollection(daily_snaps).mean().clip(aoi_geom)
    return _derive_sebal_fields(era5_mean)


def _derive_sebal_fields(era5_mean: ee.Image) -> ee.Image:
    """Derive SEBAL-friendly meteorology from raw ERA5-Land bands."""
    T_air_K = era5_mean.select("temperature_2m").rename("T_air_K")
    T_air_C = T_air_K.subtract(273.15).rename("T_air_C")

    u10 = era5_mean.select("u_component_of_wind_10m")
    v10 = era5_mean.select("v_component_of_wind_10m")
    wind_10m = u10.pow(2).add(v10.pow(2)).sqrt().rename("wind_10m")

    # Accumulated J/m^2 per hour -> mean hourly W/m^2
    Rs_in = (era5_mean.select("surface_solar_radiation_downwards_hourly")
             .divide(3600).rename("Rs_in"))
    Rl_in = (era5_mean.select("surface_thermal_radiation_downwards_hourly")
             .divide(3600).rename("Rl_in"))
    P_surf = era5_mean.select("surface_pressure").rename("P_surf")

    return ee.Image.cat([T_air_K, T_air_C, wind_10m, Rs_in, Rl_in, P_surf])


def load_era5_daily_rn(
    aoi_geom: ee.Geometry,
    start_date: ee.Date,
    end_date: ee.Date,
) -> ee.Image:
    """
    Compute the 24-hour mean net radiation from ERA5-Land over the window.

    Used as the energy driver for daily ET via the Evaporative Fraction method.

    Returns
    -------
    ee.Image
        Single-band 'Rn_24h' image in W/m^2 (window mean, not instantaneous).
    """
    vars_ = DEFAULT_PARAMS["ERA5_NET_RAD_VARS"]
    coll_id = DEFAULT_PARAMS["ERA5_COLLECTION"]

    era5_full = (ee.ImageCollection(coll_id)
                 .filterDate(start_date, end_date)
                 .select(vars_))

    daily_mean = era5_full.mean().divide(3600).clip(aoi_geom)
    Rn_24h = (daily_mean.select("surface_net_solar_radiation_hourly")
              .add(daily_mean.select("surface_net_thermal_radiation_hourly"))
              .rename("Rn_24h"))
    return Rn_24h
