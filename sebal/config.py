"""Physical constants and default parameters used across SEBAL modules."""

SEBAL_CONSTANTS = {
    "SIGMA": 5.670374419e-8,   # Stefan-Boltzmann constant (W m^-2 K^-4)
    "CP": 1004.0,              # specific heat of air (J kg^-1 K^-1)
    "R_SPEC": 287.05,          # specific gas constant for dry air (J kg^-1 K^-1)
    "K_VK": 0.41,              # von Karman constant (dimensionless)
    "G_GRAV": 9.81,            # gravitational acceleration (m s^-2)
    "WATER_DENSITY": 1000.0,   # (kg m^-3)
}

DEFAULT_PARAMS = {
    # Aerodynamic
    "Z1": 0.1,                 # lower height for dT slab (m)
    "Z2": 2.0,                 # upper height for dT slab (m)
    "Z_BLEND": 200.0,          # blending height (m)
    "Z_WIND": 10.0,            # ERA5 wind reference height (m)

    # Landsat
    "LANDSAT_CLOUD_COVER_MAX": 80,
    "LANDSAT_BANDS": ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7",
                      "ST_B10", "QA_PIXEL"],
    "LANDSAT_COLLECTIONS": ["LANDSAT/LC08/C02/T1_L2", "LANDSAT/LC09/C02/T1_L2"],

    # ERA5-Land
    "ERA5_COLLECTION": "ECMWF/ERA5_LAND/HOURLY",
    "ERA5_REFERENCE_HOUR_UTC": 10,   # approximate Landsat overpass hour
    "ERA5_VARS": [
        "temperature_2m",
        "dewpoint_temperature_2m",
        "u_component_of_wind_10m",
        "v_component_of_wind_10m",
        "surface_solar_radiation_downwards_hourly",
        "surface_thermal_radiation_downwards_hourly",
        "surface_pressure",
    ],
    "ERA5_NET_RAD_VARS": [
        "surface_net_solar_radiation_hourly",
        "surface_net_thermal_radiation_hourly",
    ],

    # Albedo (Tasumi 2008 coefficients for L8/L9 OLI)
    "ALBEDO_COEFFS": {
        "SR_B2": 0.300, "SR_B3": 0.277, "SR_B4": 0.233,
        "SR_B5": 0.143, "SR_B6": 0.036, "SR_B7": 0.012,
    },

    # Roughness length parameterization (fallback linear with NDVI)
    "Z0M_NDVI_COEF": 0.06,
    "Z0M_FLOOR": 0.005,

    # SEBAL iteration
    "N_ITER_H": 8,
    "H_FLOOR": 1.0,
    "R_AH_FLOOR": 1.0,

    # Anchor pixel selection
    "ANCHOR_PERCENTILES": [5, 20, 80, 95],
    "AG_NDVI_MIN": 0.05,
    "AG_ALBEDO_MIN": 0.10,
    "AG_ALBEDO_MAX": 0.30,

    # Soil heat flux
    "G_RATIO_MIN": 0.05,
    "G_RATIO_MAX": 0.50,
    "WATER_G_RATIO": 0.50,
}
