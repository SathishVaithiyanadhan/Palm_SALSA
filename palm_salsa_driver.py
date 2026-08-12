##parallel procesing working
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PALM SALSA Driver Generator 
=================================================
Generate aerosol emission input files for PALM's SALSA module.
- 4 emission categories: traffic exhaust, road dust, wood combustion, other
- Date/time range filtering for specific simulation periods
- Category-specific GNFR sector mapping to preserve spatial structure
- 12 aerosol species support with dynamic bin calculation
- Uses LRU caching, GDAL Warp, precompiled regex, vectorized ops

References:
- Kurppa et al. (2019) GMD - PALM-SALSA implementation
- Kokkola et al. (2008) ACP - Original SALSA description
"""

import os
import glob
import datetime
import re
import time
from datetime import timedelta, datetime
import numpy as np
from netCDF4 import Dataset
from osgeo import gdal
from pyproj import Transformer
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# USER CONFIGURATION — Loaded from YAML file passed as argument
# =============================================================================
# Usage:
#   python palm_salsa_driver.py my_scenario.yaml

import yaml
import sys

if len(sys.argv) < 2:
    print("Usage: python palm_salsa_driver.py <config.yaml>")
    sys.exit(1)

CONFIG_FILE = sys.argv[1]

with open(CONFIG_FILE, "r") as f:
    _cfg = yaml.safe_load(f)

# ---- Data Paths ----
STATIC_FILE  = _cfg["paths"]["static_file"]
TIFF_DIR     = _cfg["paths"]["tiff_dir"]
OUTPUT_FILE  = _cfg["paths"]["output_file"]

# ---- Date/Time Range for Emissions ----
START_DATE = _cfg["time"]["start_date"]
END_DATE   = _cfg["time"]["end_date"]

# ---- Emission Categories ----
ACTIVE_OUTPUT_CATEGORIES = _cfg["active_categories"]

# ---- Species Selection ----
SPECIES_OUTPUT_MODE = _cfg["species_output_mode"]
if SPECIES_OUTPUT_MODE == "custom":
    CUSTOM_SPECIES_LIST = _cfg["custom_species_list"]

# ---- Bin Parameters (must match PALM namelist) ----
NBIN   = _cfg["bins"]["nbin"]
REGLIM = _cfg["bins"]["reglim"]
NF2A   = _cfg["bins"]["nf2a"]

# =============================================================================
# CATEGORY-SPECIFIC GNFR SECTOR MAPPING
# =============================================================================

GNFR_SECTOR_TO_CATEGORY = {
    'F_RoadTransport': 0,
    'A_PublicPower': 2,
    'B_Industry': 2,
    'C_OtherStationaryComb': 2,
    'D_Fugitives': 3,
    'E_Solvents': 3,
    'G_Shipping': 3,
    'H_Aviation': 3,
    'I_OffRoad': 3,
    'J_Waste': 3,
    'K_AgriLivestock': 3,
    'L_AgriOther': 3,
}

def get_category_from_band(band_name):
    """category lookup"""
    if band_name is None:
        return [3]
    categories = []
    for sector_pattern, category in GNFR_SECTOR_TO_CATEGORY.items():
        if sector_pattern in band_name:
            categories.append(category)
    return categories if categories else [3]

# =============================================================================
# SPECIES PROPERTIES
# =============================================================================

# Species properties for mass -> number conversion
# DENSITIES must match PALM SALSA internal values (salsa_mod.f90 lines 212-219):
#   arhooc=2000, arhobc=2000, arhodu=2650, arhoh2so4=1830,
#   arhoss=2165, arhohno3=1479, arhonh3=1530
species_properties = {
    "H2SO4": {"rho": 1830, "name": "Sulfuric acid", "molar_mass": 0.098},
    "OC":    {"rho": 2000, "name": "Organic carbon", "molar_mass": 0.012},
    "BC":    {"rho": 2000, "name": "Black carbon", "molar_mass": 0.012},
    "DU":    {"rho": 2650, "name": "Mineral dust", "molar_mass": 0.060},
    "SS":    {"rho": 2165, "name": "Sea salt", "molar_mass": 0.058},
    "HNO3":  {"rho": 1479, "name": "Nitric acid", "molar_mass": 0.063},
    "NH3":   {"rho": 1530, "name": "Ammonia", "molar_mass": 0.017},
    "PB":    {"rho": 11340, "name": "Lead", "molar_mass": 0.2072},
    "HG":    {"rho": 13534, "name": "Mercury", "molar_mass": 0.20059},
    "NI":    {"rho": 8908, "name": "Nickel", "molar_mass": 0.058693},
    "CD":    {"rho": 8650, "name": "Cadmium", "molar_mass": 0.112414},
    "AS":    {"rho": 5727, "name": "Arsenic", "molar_mass": 0.074922},
}

SPECIES_2A_FRACTION = {
    "H2SO4": 0.90, "OC": 0.70, "BC": 0.10, "DU": 0.10,
    "SS": 0.90, "NH3": 0.90, "HNO3": 0.90,
    "PB": 0.50, "HG": 0.50, "NI": 0.50, "CD": 0.50, "AS": 0.50,
}

SIZE_DISTRIBUTIONS = {
    0: {"name": "traffic exhaust", "modes": [
            {"Dg": 13.5e-9, "sigma": 1.6, "weight": 0.016},
            {"Dg": 60.0e-9, "sigma": 1.8, "weight": 0.800},
            {"Dg": 2.0e-7, "sigma": 1.6, "weight": 0.100},
            {"Dg": 2.5e-6, "sigma": 1.6, "weight": 2.239e-4}],  # brake/tire wear (coarse)
    },
    1: {"name": "road dust", "modes": [
            {"Dg": 1.4e-7, "sigma": 1.4, "weight": 0.5},
            {"Dg": 4.0e-6, "sigma": 1.6, "weight": 0.5}],
    },
    2: {"name": "wood combustion", "modes": [
            {"Dg": 5.4e-8, "sigma": 1.7, "weight": 0.6},
            {"Dg": 2.5e-7, "sigma": 1.6, "weight": 0.35},
            {"Dg": 2.0e-6, "sigma": 1.5, "weight": 1.413e-3}],  # small fly-ash/coarse fraction
    },
    3: {"name": "other", "modes": [
            {"Dg": 60.0e-9, "sigma": 1.7, "weight": 0.45},
            {"Dg": 3.0e-7, "sigma": 1.6, "weight": 0.40},
            {"Dg": 1.5e-6, "sigma": 1.6, "weight": 7.499e-2}],
    }
}

# =============================================================================
# OPTIONAL SIZE DISTRIBUTION OVERRIDE (from YAML config)
# =============================================================================
# If the config file contains a "size_distributions:" block it overrides the
# per-category lognormal modes above without touching code. Only the listed
# categories are overridden; the rest keep the defaults above.
#
#   size_distributions:
#     traffic: [{Dg: 1.35e-8, sigma: 1.6, weight: 0.016}, ...]
#     wood:    [{Dg: 5.4e-8, sigma: 1.7, weight: 0.6},
#               {Dg: 4.0e-7, sigma: 1.5, weight: 0.4}]
#     other:   [{Dg: 6.0e-8, sigma: 1.7, weight: 0.5}, ...]
#
_CATEGORY_NAME_TO_IDX = {'traffic': 0, 'dust': 1, 'wood': 2, 'other': 3}

if "size_distributions" in _cfg:
    for cat_name, modes in _cfg["size_distributions"].items():
        cat_idx = _CATEGORY_NAME_TO_IDX.get(cat_name)
        if cat_idx is None or cat_idx not in SIZE_DISTRIBUTIONS:
            print(f"WARNING: ignoring unknown size_distributions category '{cat_name}'")
            continue
        new_modes = []
        for m in modes:
            try:
                new_modes.append({
                    "Dg": float(m["Dg"]),
                    "sigma": float(m["sigma"]),
                    "weight": float(m["weight"]),
                })
            except (KeyError, TypeError, ValueError) as e:
                print(f"WARNING: bad size_distributions entry for '{cat_name}': "
                      f"{m} ({e}); keeping code default")
                new_modes = None
                break
        if new_modes:
            SIZE_DISTRIBUTIONS[cat_idx]["modes"] = new_modes
            print(f"Size distribution override applied for '{cat_name}' "
                  f"(category {cat_idx})")

SPECIES_CATEGORY_MAPPING = {
    "oc": {"target": "OC"}, "ec": {"target": "BC"}, # "bc": {"target": "BC"},
    "so4": {"target": "H2SO4"}, "no": {"target": "HNO3"}, "no2": {"target": "HNO3"},
    "nh3": {"target": "NH3"}, "othmin": {"target": "DU"},
    "pb": {"target": "PB"}, "hg": {"target": "HG"}, "ni": {"target": "NI"},
    "cd": {"target": "CD"}, "as": {"target": "AS"}, "na": {"target": "SS"},
}

BULK_SPECIES = ['ho2', 'h2o', 'o3', 'ro2', 'oh', 'rcho', 'n2o', 
                'nmvoc', 'voc', 'co', 'co2', 'ch4', 'pm25', 'pmcoarse',
                'pm2_5', 'pm10', 'nox',
                'no', 'no2', 'nh3', 'so4']

CONFIG_PROJ = "EPSG:25832"
DEFAULT_PROJ = "EPSG:4326"
transformer_to_utm = Transformer.from_crs(DEFAULT_PROJ, CONFIG_PROJ, always_xy=True)
transformer_to_wgs = Transformer.from_crs(CONFIG_PROJ, DEFAULT_PROJ, always_xy=True)

# =============================================================================
# CACHE CLEANUP
# =============================================================================

def clear_geotiff_cache():
    """Clear memory after processing"""
    import gc
    gc.collect()

# Pre-compiled regex for band parsing
BAND_PATTERN = re.compile(r'^([A-Za-z_]+)_h(\d{2})_(\d{8})$')

# =============================================================================
# DATE/TIME PARSING FUNCTIONS
# =============================================================================

def parse_date_range(start_str, end_str):
    """Parse date range and calculate hours"""
    start_dt = datetime.strptime(start_str, "%Y-%m-%d %H:%M:%S")
    end_dt = datetime.strptime(end_str, "%Y-%m-%d %H:%M:%S")
    time_diff = end_dt - start_dt
    n_hours = int(time_diff.total_seconds() / 3600) + 1
    print(f"Date range: {start_dt} to {end_dt}")
    print(f"Number of hours: {n_hours}")
    return start_dt, end_dt, n_hours


def create_all_time_steps(start_dt, end_dt):
    """Pre-compute ALL expected time steps"""
    current_dt = start_dt
    time_steps = []
    while current_dt <= end_dt:
        time_steps.append({
            'date': current_dt.strftime("%Y%m%d"),
            'hour': current_dt.hour,
            'hour_str': f"h{current_dt.hour:02d}",
            'datetime': current_dt
        })
        current_dt += timedelta(hours=1)
    return time_steps


def parse_band_description(desc):
    """band parsing with pre-compiled regex"""
    match = BAND_PATTERN.match(desc)
    if match:
        sector, hour_str, date_str = match.groups()
        hour = int(hour_str)
        band_dt = datetime(int(date_str[0:4]), int(date_str[4:6]), 
                          int(date_str[6:8]), hour, 0, 0)
        return sector, hour, date_str, band_dt
    return None


# =============================================================================
# BIN CALCULATION FUNCTIONS 
# =============================================================================

def calculate_palm_bin_diameters(reglim, nbin, nf2a):
    """Calculate bin diameters exactly as PALM does"""
    d_min, d_split, d_max = reglim
    nbin1, nbin2 = nbin
    
    has_2a = (nf2a > 0.0)
    has_2b = (nf2a < 1.0)
    
    if has_2a and has_2b:
        nbins_total = nbin1 + nbin2 + nbin2
    elif has_2a:
        nbins_total = nbin1 + nbin2
    elif has_2b:
        nbins_total = nbin1 + nbin2
    else:
        nbins_total = nbin1
    
    dmid = np.zeros(nbins_total)
    dlow = np.zeros(nbins_total)
    dhigh = np.zeros(nbins_total)
    
    labels = ['1'] * nbin1
    if has_2a:
        labels.extend(['2a'] * nbin2)
    if has_2b:
        labels.extend(['2b'] * nbin2)
    subrange_labels = np.array(labels, dtype=str)
    
    ratio_1 = d_split / d_min
    for i in range(nbin1):
        d_low = d_min * ratio_1**(i / nbin1)
        d_high = d_min * ratio_1**((i+1) / nbin1)
        dlow[i] = d_low; dhigh[i] = d_high
        dmid[i] = np.sqrt(d_low * d_high)
    
    ratio_2 = d_max / d_split
    bin_offset = nbin1
    
    if has_2a:
        for i in range(nbin2):
            idx = bin_offset + i
            d_low = d_split * ratio_2**(i / nbin2)
            d_high = d_split * ratio_2**((i+1) / nbin2)
            dlow[idx] = d_low; dhigh[idx] = d_high
            dmid[idx] = np.sqrt(d_low * d_high)
        bin_offset += nbin2
    
    if has_2b:
        for i in range(nbin2):
            idx = bin_offset + i
            d_low = d_split * ratio_2**(i / nbin2)
            d_high = d_split * ratio_2**((i+1) / nbin2)
            dlow[idx] = d_low; dhigh[idx] = d_high
            dmid[idx] = np.sqrt(d_low * d_high)
    
    return dmid, dlow, dhigh, subrange_labels, has_2a, has_2b, nbins_total


def lognormal_pdf(d, dg, sigma_g):
    """Log-normal PDF"""
    d_safe = np.maximum(d, 1e-30)
    return (1.0 / (d_safe * np.log(sigma_g) * np.sqrt(2 * np.pi))) * \
           np.exp(-(np.log(d_safe / dg)**2) / (2 * np.log(sigma_g)**2))


def get_size_distribution_fractions(category, bin_diameters, subrange_labels, nf2a, has_2a, has_2b):
    """Calculate number fractions for a given emission category.
    
    Uses the category-level size distribution (applies to ALL species
    within that category).  The 2a/2b soluble/insoluble split is applied
    later in the worker via SPECIES_2A_FRACTION.
    """
    cat_info = SIZE_DISTRIBUTIONS.get(category)
    if cat_info is None:
        raise ValueError(f"Unknown category: {category}")
    
    modes = cat_info["modes"]
    pdf_values = np.zeros(len(bin_diameters))
    for mode in modes:
        pdf_values += mode["weight"] * lognormal_pdf(bin_diameters, mode["Dg"], mode["sigma"])
    
    total = np.sum(pdf_values)
    if total > 0:
        return pdf_values / total
    # Fallback: uniform distribution
    return np.ones(len(bin_diameters)) / len(bin_diameters)


def extract_static_domain(static_nc):
    """Extract domain parameters from static file"""
    params = {}
    params['origin_time'] = static_nc.getncattr('origin_time')
    params['origin_lat'] = static_nc.getncattr('origin_lat')
    params['origin_lon'] = static_nc.getncattr('origin_lon')
    
    params['nx'] = len(static_nc.dimensions['x'])
    params['ny'] = len(static_nc.dimensions['y'])
    
    x_coords = static_nc.variables['x'][:]
    y_coords = static_nc.variables['y'][:]
    
    params['dx'] = x_coords[1] - x_coords[0] if len(x_coords) > 1 else 1.0
    params['dy'] = abs(y_coords[1] - y_coords[0]) if len(y_coords) > 1 else 1.0
    
    # Absolute UTM coordinates that the (relative) x/y arrays are referenced to.
    # Prefer the static file's origin_x / origin_y attributes; fall back to
    # transforming origin_lon/origin_lat if they are absent.
    try:
        origin_x_abs = static_nc.getncattr('origin_x')
    except AttributeError:
        origin_x_abs, _ = transformer_to_utm.transform(
            params['origin_lon'], params['origin_lat'])
    try:
        origin_y_abs = static_nc.getncattr('origin_y')
    except AttributeError:
        _, origin_y_abs = transformer_to_utm.transform(
            params['origin_lon'], params['origin_lat'])
    
    # Derive the domain bounds from the actual grid coordinates (x/y are the
    # cell-CENTRE offsets from the origin).  This keeps the GDAL Warp window
    # aligned with the PALM grid even when x/y start at 0 or at dx/2.  (The
    # previous code assumed x/y were centred on the origin point, which is off
    # by ~half the domain for the 128x128 statics -> emissions landed in the
    # wrong quadrant.)
    params['west']  = origin_x_abs + x_coords[0]  - params['dx'] / 2.0
    params['east']  = origin_x_abs + x_coords[-1] + params['dx'] / 2.0
    if y_coords[-1] >= y_coords[0]:      # y increasing northwards
        params['south'] = origin_y_abs + y_coords[0]  - params['dy'] / 2.0
        params['north'] = origin_y_abs + y_coords[-1] + params['dy'] / 2.0
    else:                                 # y decreasing (south first)
        params['south'] = origin_y_abs + y_coords[-1] - params['dy'] / 2.0
        params['north'] = origin_y_abs + y_coords[0]  + params['dy'] / 2.0
    params['origin_x'] = params['west']
    params['origin_y'] = params['north']
    
    params['lon_w'], params['lat_s'] = transformer_to_wgs.transform(
        params['west'], params['south'])
    params['lon_e'], params['lat_n'] = transformer_to_wgs.transform(
        params['east'], params['north'])
    
    return params


# =============================================================================
# MULTIPROCESSING WORKER: Module-level functions (pickle-safe)
# =============================================================================

_WORKER_PARAMS = None

def _init_worker(static_params, nx, ny, ntime, bin_diameters, subrange_labels,
                  size_distributions, nf2a, has_2a, has_2b, nbins_total,
                  start_dt, end_dt, time_steps, time_index_lookup,
                  selected_cat_indices, composition_name_list):
    """Initialize per-process worker data (runs once per worker).
    
    Stores all shared read-only data as module-level globals so that
    the worker function _process_tiff_file_wrapper can access them
    without pickling a TiffProcessor instance.
    """
    global _WORKER_PARAMS
    import numpy as np
    from osgeo import gdal
    
    # Build a TiffProcessor-like state dict
    _WORKER_PARAMS = {
        'static_params': static_params,
        'nx': nx, 'ny': ny, 'ntime': ntime,
        'bin_diameters': bin_diameters,
        'subrange_labels': subrange_labels,
        'size_distributions': size_distributions,
        'nf2a': nf2a, 'has_2a': has_2a, 'has_2b': has_2b,
        'nbins_total': nbins_total,
        'start_dt': start_dt, 'end_dt': end_dt,
        'time_steps': time_steps,
        'time_index_lookup': time_index_lookup,
        'selected_cat_indices': selected_cat_indices,
        'composition_name_list': composition_name_list,
    }
    
    # Pre-compute bin volumes and conversion factor for this worker
    _WORKER_PARAMS['bin_volumes'] = np.array([
        (4.0/3.0) * np.pi * (d/2.0)**3 for d in bin_diameters
    ], dtype=np.float32)
    _WORKER_PARAMS['conv_factor'] = np.float32(1.0 / 3600.0)
    
    # Per-worker local cache (module-level, private to this process)
    _WORKER_PARAMS['_geotiff_cache'] = {}
    _WORKER_PARAMS['_cache_order'] = []
    _WORKER_PARAMS['_MAX_CACHE_SIZE'] = 3


def _process_tiff_file_wrapper(tiff_file):
    """Module-level wrapper for processing a single TIFF file.
    
    This is a standalone function so multiprocessing.Pool can pickle it
    without serializing a TiffProcessor instance.  Uses _WORKER_PARAMS
    that was set by _init_worker in each worker process.
    """
    import os
    import numpy as np
    from osgeo import gdal
    
    p = _WORKER_PARAMS
    cache = p['_geotiff_cache']
    cache_order = p['_cache_order']
    max_cache = p['_MAX_CACHE_SIZE']
    
    filename = os.path.basename(tiff_file).lower()
    
    try:
        # --- identical logic to TiffProcessor.process_single_file ---
        
        # Skip bulk species
        if any(bulk in filename for bulk in BULK_SPECIES):
            return None
        
        # Find matching species
        species = None
        for key in SPECIES_CATEGORY_MAPPING.keys():
            if key in filename:
                species = key
                break
        if species is None:
            return None
        
        mapping = SPECIES_CATEGORY_MAPPING.get(species)
        if not mapping:
            return None
        
        target_species = mapping["target"]
        
        # Skip TIFF files whose target species is not in the output composition list
        # (e.g., skip PB, HG, NI, CD, AS if using 7-species config)
        if target_species not in p['composition_name_list']:
            print(f"  Skipping {filename}: species {target_species} not in output list")
            return None
        
        # Get band info (open file, parse bands within date range)
        ds = gdal.Open(tiff_file, gdal.GA_ReadOnly)
        if not ds:
            return None
        
        band_info = {}
        total_bands = ds.RasterCount
        for band_num in range(1, total_bands + 1):
            band = ds.GetRasterBand(band_num)
            desc = band.GetDescription()
            if desc:
                parsed = parse_band_description(desc)
                if parsed is not None:
                    sector, hour, date_str, band_dt = parsed
                    if p['start_dt'] <= band_dt <= p['end_dt']:
                        time_key = (date_str, hour)
                        time_idx = p['time_index_lookup'].get(time_key)
                        if time_idx is not None:
                            categories = get_category_from_band(desc)
                            if time_idx not in band_info:
                                band_info[time_idx] = []
                            band_info[time_idx].append({
                                'band_num': band_num,
                                'categories': categories
                            })
        ds = None
        
        if not band_info:
            print(f"\n  Skipping {filename}: No bands in date range")
            return None
        
        total_bands_to_process = sum(len(bands) for bands in band_info.values())
        print(f"\n{'='*60}")
        print(f"File: {filename}")
        print(f"  Species: {species} -> {target_species}")
        print(f"  Bands to process: {total_bands_to_process} (across {len(band_info)} time steps)")
        
        # Only allocate categories that actually appear in this file's bands
        # AND are in the user's selected output categories
        # (saves memory by avoiding unused category arrays)
        needed_categories = set()
        for bands_list in band_info.values():
            for entry in bands_list:
                for cat in entry['categories']:
                    needed_categories.add(cat)
        # Filter to only keep selected categories
        needed_categories = needed_categories & set(p['selected_cat_indices'])
        if not needed_categories:
            print(f"  Skipping {filename}: No bands match selected categories")
            return None
        category_emissions = {}
        for cat in needed_categories:
            category_emissions[cat] = np.zeros(
                (p['ntime'], p['ny'], p['nx'], p['nbins_total']), dtype=np.float32
            )
        
        f_2a = SPECIES_2A_FRACTION.get(target_species, 0.5)
        
        # --- Resample GeoTIFF (per-worker local cache) ---
        if tiff_file in cache:
            cache_order.remove(tiff_file)
            cache_order.insert(0, tiff_file)
            all_bands_data = cache[tiff_file]
        else:
            if not os.path.exists(tiff_file):
                raise FileNotFoundError(f"GeoTIFF file not found: {tiff_file}")
            print(f"  Resampling: {os.path.basename(tiff_file)}")
            warp_options = gdal.WarpOptions(
                format='MEM',
                outputBounds=[p['static_params']['west'], p['static_params']['south'],
                              p['static_params']['east'], p['static_params']['north']],
                width=p['static_params']['nx'],
                height=p['static_params']['ny'],
                dstSRS=CONFIG_PROJ,
                resampleAlg=gdal.GRA_Average   # mass-conserving resampling
            )
            resampled_ds = gdal.Warp('', tiff_file, options=warp_options)
            if resampled_ds is None:
                raise RuntimeError(f"Failed to resample GeoTIFF: {tiff_file}")
            num_bands = resampled_ds.RasterCount
            all_bands_data = []
            for bnum in range(1, num_bands + 1):
                arr = resampled_ds.GetRasterBand(bnum).ReadAsArray().astype(np.float32)
                arr = np.flipud(arr)
                all_bands_data.append(arr)
            resampled_ds = None
            
            # LRU cache per worker
            if len(cache_order) >= max_cache:
                oldest = cache_order.pop()
                del cache[oldest]
            cache[tiff_file] = all_bands_data
            cache_order.insert(0, tiff_file)
        
        # --- Process bands (logic identical to original) ---
        bands_processed = 0
        
        # Initialize mass tracking for this file (across all bands)
        mass_sums = {}
        
        for time_idx, bands in band_info.items():
            if time_idx >= p['ntime']:
                continue
            for band_entry in bands:
                band_num = band_entry['band_num']
                categories = band_entry['categories']
                
                if band_num - 1 < len(all_bands_data):
                    mass_data = all_bands_data[band_num - 1]
                else:
                    continue
                if np.all(mass_data == 0):
                    continue
                
                # Convert units: kg/m2/hr -> kg/m2/s
                mass_data = mass_data * p['conv_factor']
                bands_processed += 1
                
                for cat in categories:
                    # Track total mass for dynamic mass fraction computation
                    # (accumulated per category + species across all bands)
                    if cat not in mass_sums:
                        mass_sums[cat] = {}
                    if target_species not in mass_sums[cat]:
                        mass_sums[cat][target_species] = 0.0
                    mass_sums[cat][target_species] += np.sum(mass_data)
                    
                    # ── Category-level size distribution ──
                    # Every species in this category shares the same size distribution
                    size_fracs = p['size_distributions'][cat]  # pre-computed array
                    species_mass = mass_data
                    
                    adjusted_fracs = size_fracs.copy()
                    if p['has_2a'] and p['has_2b']:
                        mask_2a = p['subrange_labels'] == '2a'
                        mask_2b = p['subrange_labels'] == '2b'
                        adjusted_fracs[mask_2a] *= f_2a
                        adjusted_fracs[mask_2b] *= (1.0 - f_2a)
                        # Re-normalize to preserve total mass after 2a/2b split
                        total_frac = np.sum(adjusted_fracs)
                        if total_frac > 0:
                            adjusted_fracs /= total_frac
                    
                    # Compute number-weighted average particle volume
                    # This is the correct V_avg = Σ(nf_i × V_i) for the size distribution
                    V_avg = np.sum(adjusted_fracs * p['bin_volumes'])
                    rho = species_properties[target_species]["rho"]
                    
                    # Total number flux = total_mass / (density × average_volume)
                    # Distribute by number fractions to each bin
                    with np.errstate(divide='ignore', invalid='ignore'):
                        total_number_flux = np.where(
                            (V_avg > 0) & (species_mass > 0),
                            species_mass / (rho * V_avg), 0.0
                        )
                    
                    for bin_idx in range(p['nbins_total']):
                        if adjusted_fracs[bin_idx] <= 0:
                            continue
                        number_flux = total_number_flux * adjusted_fracs[bin_idx]
                        category_emissions[cat][time_idx, :, :, bin_idx] += number_flux
        
        print(f"  Processed: {bands_processed} bands")
        
        if bands_processed == 0:
            return None
        
        # Build result dict (only categories that were allocated)
        results = {}
        for cat, cat_data in category_emissions.items():
            if np.any(cat_data > 0):
                results[f"cat_{cat}"] = cat_data
        results['species'] = species
        results['target'] = target_species
        results['mass_sums'] = mass_sums  # for dynamic mass fraction computation
        return results
    
    except Exception as e:
        print(f"\n  ERROR processing {filename}: {e}")
        return None


# =============================================================================
# MAIN SALSA DRIVER CLASS
# =============================================================================

class SalsaDriver:
    """Generate PALM SALSA aerosol emission driver"""
    
    def __init__(self, static_file, tiff_dir, output_file, active_categories=None):
        print("=" * 70)
        print("PALM-SALSA Driver Generator")
        print("=" * 70)
        
        self.active_output_categories = ACTIVE_OUTPUT_CATEGORIES
        self.species_output_mode = SPECIES_OUTPUT_MODE
        self.nbin = NBIN
        self.reglim = REGLIM
        self.nf2a = NF2A
        
        # Parse date range
        self.start_dt, self.end_dt, self.ntime = parse_date_range(START_DATE, END_DATE)
        
        # Pre-compute time steps
        self.time_steps = create_all_time_steps(self.start_dt, self.end_dt)
        print(f"Pre-computed {len(self.time_steps)} time steps")
        
        # 4 CATEGORIES
        self.category_name_to_idx = {'traffic': 0, 'dust': 1, 'wood': 2, 'other': 3}
        
        # Determine species list
        if SPECIES_OUTPUT_MODE == "custom" and 'CUSTOM_SPECIES_LIST' in globals():
            self.composition_name_list = CUSTOM_SPECIES_LIST
        elif SPECIES_OUTPUT_MODE == "basic7":
            self.composition_name_list = ["H2SO4", "OC", "BC", "DU", "SS", "HNO3", "NH3"]
        else:
            self.composition_name_list = ["H2SO4", "OC", "BC", "DU", "SS", "HNO3", "NH3",
                                           "PB", "HG", "NI", "CD", "AS"]
        
        print(f"\nConfiguration:")
        print(f"  Date range: {self.start_dt} to {self.end_dt}")
        print(f"  Time steps: {self.ntime}")
        print(f"  nbin = {self.nbin}")
        print(f"  nf2a = {self.nf2a}")
        
        # Calculate bin structure
        (self.bin_diameters, self.bin_low, self.bin_high, 
         self.subrange_labels, self.has_2a, self.has_2b, self.nbins_total) = \
            calculate_palm_bin_diameters(self.reglim, self.nbin, self.nf2a)
        
        print(f"\nBin structure: {self.nbins_total} bins "
              f"(1: {self.nbin[0]}, 2a: {self.nbin[1]}, 2b: {self.nbin[1] if self.has_2b else 0})")
        
        self.selected_cat_indices = [self.category_name_to_idx[name] 
                                      for name in self.active_output_categories]
        self.selected_cat_names = [SIZE_DISTRIBUTIONS[idx]['name'] 
                                    for idx in self.selected_cat_indices]
        
        print(f"\nSelected categories: {self.selected_cat_names}")
        
        print("\nOpening static driver...")
        self.static_nc = Dataset(static_file, "r")
        self.output_file = output_file
        self.tiff_dir = tiff_dir
        
        self.static_params = extract_static_domain(self.static_nc)
        self.nx = self.static_params['nx']
        self.ny = self.static_params['ny']
        self.static_x = self.static_nc.variables["x"][:]
        self.static_y = self.static_nc.variables["y"][:]
        
        print(f"Domain: {self.nx} x {self.ny}")

        # Print all bins
        print(f"\nBin details ({self.nbins_total} bins):")
        for i, (d, label) in enumerate(zip(self.bin_diameters, self.subrange_labels)):
            print(f"  Bin {i+1:2d} [{label}]: {d*1e9:.1f} nm "
                  f"[{self.bin_low[i]*1e9:.1f}-{self.bin_high[i]*1e9:.1f}] nm")
        
        # Calculate category-level size distributions (one per category, shared by all species)
        print("\nCalculating size distributions (category-level)...")
        self.size_distributions = {}
        for cat in [0, 1, 2, 3]:
            dist = get_size_distribution_fractions(
                cat, self.bin_diameters, self.subrange_labels,
                self.nf2a, self.has_2a, self.has_2b
            )
            self.size_distributions[cat] = dist
        
        print(f"\nCreating output file: {output_file}")
        self.nc_file = Dataset(output_file, "w", format="NETCDF4")
        
        self.define_dimensions()
        self.write_global_attributes()
        self.add_variables()
        self.validate_emissions()
        self.finalize()
    
    def write_global_attributes(self):
        """Write global attributes"""
        print("\nWriting global attributes...")
        for attr in self.static_nc.ncattrs():
            try:
                setattr(self.nc_file, attr, self.static_nc.getncattr(attr))
            except:
                pass
        
        self.nc_file.creation_date = str(datetime.now())
        self.nc_file.description = (
            f"Aerosol input (SALSA driver) for PALM - 4 categories. "
            f"Period: {START_DATE} to {END_DATE}. "
            f"nbin={self.nbin}, nf2a={self.nf2a}, bins={self.nbins_total}"
        )
        self.nc_file.title = "PALM SALSA aerosol input"
        self.nc_file.institution = "University of Augsburg"
        self.nc_file.author = "Sathish Kumar Vaithiyanadhan"
        self.nc_file.emission_start = START_DATE
        self.nc_file.emission_end = END_DATE
        self.nc_file.composition_name = ', '.join(self.composition_name_list)
        self.nc_file.nbin = f"{self.nbin[0]}, {self.nbin[1]}"
        self.nc_file.nf2a = f"{self.nf2a}"
        self.nc_file.total_bins = f"{self.nbins_total}"

    def define_dimensions(self):
        """Define NetCDF dimensions"""
        print("Defining dimensions...")
        self.nncat = len(self.selected_cat_indices)
        self.ncomposition_index = len(self.composition_name_list)
        self.nmax_string_length = 25
        
        self.nc_file.createDimension("x", self.nx)
        self.nc_file.createDimension("y", self.ny)
        self.nc_file.createDimension("time", self.ntime)
        self.nc_file.createDimension("ncat", self.nncat)
        self.nc_file.createDimension("composition_index", self.ncomposition_index)
        self.nc_file.createDimension("max_string_length", self.nmax_string_length)
        self.nc_file.createDimension("Dmid", self.nbins_total)

        x = self.nc_file.createVariable("x", "f4", ("x",))
        x[:] = self.static_x; x.units = "m"

        y = self.nc_file.createVariable("y", "f4", ("y",))
        y[:] = self.static_y; y.units = "m"

        t = self.nc_file.createVariable("time", "f4", ("time",))
        t[:] = np.arange(0, self.ntime * 3600, 3600); t.units = "s"
        #t = self.nc_file.createVariable("time", "f4", ("time",))
        #t[:] = np.arange(1, self.ntime + 1) * 3600  # Start from 3600 (1 hour), not 0
        #t.units = "s"

        Dmid = self.nc_file.createVariable("Dmid", "f4", ("Dmid",))
        Dmid[:] = self.bin_diameters; Dmid.units = "m"

        ncat_coord = self.nc_file.createVariable("ncat", "i4", ("ncat",))
        ncat_coord[:] = np.arange(1, self.nncat + 1)
        
        comp_index = self.nc_file.createVariable("composition_index", "i4", 
                                                  ("composition_index",))
        comp_index[:] = np.arange(1, self.ncomposition_index + 1)
        
        max_str = self.nc_file.createVariable("max_string_length", "i4", 
                                               ("max_string_length",))
        max_str[:] = np.arange(1, self.nmax_string_length + 1)

    def add_variables(self):
        """Add emission variables"""
        print("\nAdding emission variables...")

        all_cat_names = ["traffic exhaust", "road dust", "wood combustion", "other"]
        filtered_names = [all_cat_names[i] for i in self.selected_cat_indices]
        
        nc_cat_name = self.nc_file.createVariable("emission_category_name", "S1", 
                                                   ("ncat", "max_string_length"))
        for i, name in enumerate(filtered_names):
            chars = list(name.ljust(self.nmax_string_length))
            nc_cat_name[i, :] = np.array(list(chars), dtype="S1")

        nc_cat_idx = self.nc_file.createVariable("emission_category_index", "i1", ("ncat",))
        nc_cat_idx[:] = np.arange(1, self.nncat + 1)

        nc_comp_name = self.nc_file.createVariable("composition_name", "S1",
                                                    ("composition_index", "max_string_length"))
        for i, name in enumerate(self.composition_name_list):
            chars = list(name.ljust(self.nmax_string_length))
            nc_comp_name[i, :] = np.array(list(chars), dtype="S1")

        # Create emission_mass_fracs variable (filled dynamically after TIFF processing)
        nc_mass_fracs = self.nc_file.createVariable("emission_mass_fracs", "f4", 
                                                    ("ncat", "composition_index"), 
                                                    fill_value=-9999.0)
        nc_mass_fracs.units = "1"

        # Create emission_number_fracs — filled AFTER emission generation
        # with the NUMBER-WEIGHTED average (not equal-weight) so PALM
        # reconstructs correct mass from the combined multi-species flux
        nc_num_fracs = self.nc_file.createVariable("emission_number_fracs", "f4", 
                                                    ("ncat", "Dmid"), 
                                                    fill_value=-9999.0)
        nc_num_fracs.units = "1"

        # Emission values
        nc_aerosol = self.nc_file.createVariable("aerosol_emission_values", "f4",
                                                  ("time", "y", "x", "ncat"), 
                                                  fill_value=-9999.0)
        nc_aerosol.units = "#/m2/s"; nc_aerosol.lod = 2

        self.generate_emission_data(nc_aerosol, nc_num_fracs)

    def generate_emission_data(self, nc_var, nc_num_fracs):
        """Generate emission data"""
        print("\n" + "=" * 70)
        print("PROCESSING EMISSION TIFF FILES")
        print(f"Date range: {self.start_dt} to {self.end_dt}")
        print("=" * 70)
        
        tiff_patterns = [os.path.join(self.tiff_dir, "emission_*_temporal.tif")]
        tiff_files = []
        for pattern in tiff_patterns:
            files = glob.glob(pattern)
            if files:
                tiff_files.extend(files)
        
        tiff_files = list(set(tiff_files))
        print(f"\nFound {len(tiff_files)} unique TIFF files")
        
        if len(tiff_files) == 0:
            print("WARNING: No TIFF files found!")
            nc_var[:] = np.zeros((self.ntime, self.ny, self.nx, self.nncat))
            return
        
# Cap workers to prevent OOM: each worker allocates ~1 GB for a 512x512x26x10 grid
        # (4 category arrays x 260 MB each + resampled GeoTIFF data).
        # With too many parallel workers, the system runs out of memory and kills the process.
        num_cpus = cpu_count()
        # Safe limit: each file result returned can be ~260 MB (one non-zero category);
        # use at most 4 workers to keep peak memory under ~5 GB
        max_safe_workers = min(4, num_cpus)
        num_processes = max(1, min(max_safe_workers, len(tiff_files)))
        print(f"Using {num_processes} process{'es' if num_processes > 1 else ''}" +
              (f" (capped from {num_cpus} to prevent OOM)" if num_processes < num_cpus else ""))
        
        # Pre-compute time index lookup for workers
        time_index_lookup = {}
        for idx, ts in enumerate(self.time_steps):
            key = (ts['date'], ts['hour'])
            time_index_lookup[key] = idx
        
        # Initialize accumulators (only for selected categories)
        category_accumulators = {}
        for cat in self.selected_cat_indices:
            category_accumulators[cat] = np.zeros((self.ntime, self.ny, self.nx, self.nbins_total), 
                                                  dtype=np.float32)
        
        processed_files = 0
        failed_files = 0
        
        # Aggregate mass sums across all TIFF files for dynamic mass fraction computation
        # Structure: mass_sums_agg[cat][species_name] = total_mass (kg/m2/s summed over domain+time)
        mass_sums_agg = {}
        
        with Pool(
            processes=num_processes,
            initializer=_init_worker,
            initargs=(
                self.static_params, self.nx, self.ny, self.ntime,
                self.bin_diameters, self.subrange_labels,
                self.size_distributions, self.nf2a, self.has_2a, self.has_2b,
                self.nbins_total, self.start_dt, self.end_dt,
                self.time_steps, time_index_lookup,
                self.selected_cat_indices,
                self.composition_name_list,
            ),
            # Restart worker after each file to free ~1 GB of accumulators + cached
            # resampled data.  Without this, memory accumulates across files.
            maxtasksperchild=1,
        ) as pool:
            # Use imap_unordered for streaming progress; each file is independent
            for result in pool.imap_unordered(_process_tiff_file_wrapper, tiff_files):
                if result is not None:
                    processed_files += 1
                    # Aggregate bin accumulators
                    for cat_key, cat_data in result.items():
                        if cat_key.startswith('cat_'):
                            cat_idx = int(cat_key.split('_')[1])
                            category_accumulators[cat_idx] += cat_data
                    # Aggregate mass sums for dynamic fraction computation
                    if 'mass_sums' in result:
                        for cat, species_dict in result['mass_sums'].items():
                            if cat not in mass_sums_agg:
                                mass_sums_agg[cat] = {}
                            for species_name, mass_val in species_dict.items():
                                mass_sums_agg[cat][species_name] = \
                                    mass_sums_agg[cat].get(species_name, 0.0) + mass_val
                else:
                    failed_files += 1
        
        print(f"\nProcessed {processed_files} files successfully" +
              (f", {failed_files} failed" if failed_files else ""))
        
        # ============================================================
        # DYNAMIC MASS FRACTION COMPUTATION
        # ============================================================
        # Compute mass fractions from actual TIFF data:
        #   mass_frac[cat][species] = total_mass[cat][species] / total_mass[cat]
        # ============================================================
        computed_mass_fracs = np.zeros(
            (self.nncat, self.ncomposition_index), dtype=np.float64
        )
        
        print(f"\n{'='*60}")
        print("COMPUTED FRACTIONS (from TIFF data)")
        print("="*60)
        
        # --- Mass fractions (per species, per category) ---
        print(f"\n  --- MASS FRACTIONS ---")
        print(f"  (fraction of total category mass attributed to each species)")
        for new_cat_idx, old_cat_idx in enumerate(self.selected_cat_indices):
            cat_name = self.selected_cat_names[new_cat_idx]
            cat_mass_sums = mass_sums_agg.get(old_cat_idx, {})
            # Only sum mass for species that are in the output composition list
            # (trace metals like PB, HG, etc. are not in the output, so exclude them)
            total_cat_mass = sum(cat_mass_sums.get(name, 0.0) 
                                 for name in self.composition_name_list)
            
            if total_cat_mass > 0:
                print(f"\n  Category {old_cat_idx} ({cat_name}):")
                print(f"    Total mass: {total_cat_mass:.6e} kg/m2/s "
                      f"(output species only)")
                for comp_idx, comp_name in enumerate(self.composition_name_list):
                    species_mass = cat_mass_sums.get(comp_name, 0.0)
                    frac = species_mass / total_cat_mass
                    computed_mass_fracs[new_cat_idx, comp_idx] = frac
                    if frac > 0:
                        print(f"    {comp_name:<8}: {frac:.6f}  ({species_mass:.6e} kg/m2/s)")
                mass_sum_check = np.sum(computed_mass_fracs[new_cat_idx, :])
                print(f"    {'Sum':<8}: {mass_sum_check:.10f}  (normalized: {abs(mass_sum_check - 1.0) < 1e-6})")
            else:
                print(f"\n  Category {old_cat_idx} ({cat_name}): NO EMISSION DATA")
        
        # Write computed mass fractions to NetCDF variable
        self.nc_file.variables['emission_mass_fracs'][:] = \
            computed_mass_fracs.astype(np.float32)
        
        # ============================================================
        
        # Create output
        aerosol_emission_values = np.zeros((self.ntime, self.ny, self.nx, self.nncat), 
                                           dtype=np.float32)
        
        print(f"\nEmission totals:")
        n_total_cells = self.ntime * self.nx * self.ny
        for new_cat_idx, old_cat_idx in enumerate(self.selected_cat_indices):
            aerosol_emission_values[:, :, :, new_cat_idx] = np.sum(
                category_accumulators[old_cat_idx], axis=3)
            
            total = np.sum(aerosol_emission_values[:, :, :, new_cat_idx])
            avg_per_cell = total / n_total_cells
            max_val = np.max(aerosol_emission_values[:, :, :, new_cat_idx])
            print(f"  {self.selected_cat_names[new_cat_idx]:<20}:")
            print(f"    Total:   {total:.2e} #/s (sum over all times × cells)")
            print(f"    Average: {avg_per_cell:.2e} #/m²/s (per cell per timestep)")
            print(f"    Max:     {max_val:.2e} #/m²/s")
        
        nc_var[:] = aerosol_emission_values
        
        # ============================================================
        # NUMBER-WEIGHTED EMISSION NUMBER FRACTIONS
        # ============================================================
        # PALM's LOD=2 format has ONE nf per category, but we process
        # each species with its own nf.  For PALM to reconstruct mass
        # correctly it needs the NUMBER-WEIGHTED average:
        #   nf_file[bin] = Σ_species(N_species × nf_species[bin]) / N_total
        # where N_species is the total number emitted per species.
        # This is derived from the per-bin accumulators directly.
        emission_num_fracs = np.zeros((self.nncat, self.nbins_total), dtype=np.float64)
        for new_cat_idx, old_cat_idx in enumerate(self.selected_cat_indices):
            # Sum over time, y, x to get total number per bin
            total_per_bin = np.sum(category_accumulators[old_cat_idx], axis=(0, 1, 2))
            total_per_bin = total_per_bin.astype(np.float64)
            bin_sum = np.sum(total_per_bin)
            if bin_sum > 0:
                nf_weighted = total_per_bin / bin_sum
            else:
                nf_weighted = np.ones(self.nbins_total) / self.nbins_total
            emission_num_fracs[new_cat_idx, :] = nf_weighted
            
            # Print per-bin number fractions for validation
            cat_name = self.selected_cat_names[new_cat_idx]
            print(f"\n  --- NUMBER FRACTIONS (number-weighted, cat {old_cat_idx}: {cat_name}) ---")
            for bin_idx in range(self.nbins_total):
                frac = nf_weighted[bin_idx]
                if frac > 0:
                    d_nm = self.bin_diameters[bin_idx] * 1e9
                    label = self.subrange_labels[bin_idx]
                    print(f"    Bin {bin_idx+1:2d} [{label}]: {frac:.6f}  ({d_nm:.1f} nm)")
            print(f"    {'Sum':<8}: {np.sum(nf_weighted):.10f}")
        
        nc_num_fracs[:] = emission_num_fracs.astype(np.float32)
        
        # Clear cache
        clear_geotiff_cache()

    def validate_emissions(self):
        """Validate emission values"""
        print("\n" + "=" * 70)
        print("VALIDATION")
        print("=" * 70)
        
        with Dataset(self.output_file, "r") as nc:
            emissions = nc.variables["aerosol_emission_values"][:]
            mass_fracs = nc.variables["emission_mass_fracs"][:]
            num_fracs = nc.variables["emission_number_fracs"][:]
            
            print(f"  Shape: {emissions.shape}")
            print(f"\n  Mass fraction sums:")
            for cat in range(self.nncat):
                print(f"    {self.selected_cat_names[cat]:<20}: {np.sum(mass_fracs[cat, :]):.6f}")
            print(f"\n  Number fraction sums:")
            for cat in range(self.nncat):
                print(f"    {self.selected_cat_names[cat]:<20}: {np.sum(num_fracs[cat, :]):.6f}")

    def finalize(self):
        """Close files"""
        print("\nFinalizing...")
        self.static_nc.close()
        self.nc_file.close()
        print(f"Created: {self.output_file}")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    total_start = time.time()
    
    print("\n" + "=" * 70)
    print("STARTING PALM-SALSA DRIVER GENERATION")
    print("=" * 70)
    print(f"Config: {CONFIG_FILE}")
    print(f"Date: {START_DATE} to {END_DATE}")
    print(f"Categories: {ACTIVE_OUTPUT_CATEGORIES}")
    print(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    driver = SalsaDriver(STATIC_FILE, TIFF_DIR, OUTPUT_FILE)
    
    total_time = time.time() - total_start
    print("\n" + "=" * 70)
    print(f"COMPLETED in {timedelta(seconds=int(total_time))}")
    print("=" * 70)