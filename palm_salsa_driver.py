#12 species
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PALM SALSA Driver Generator
============================
Generate aerosol emission input files for PALM's SALSA module using the methodology
described in Kurppa et al. (2019) "Implementation of the sectional aerosol module 
SALSA2.0 into the PALM model system 6.0"

MODIFIED: Extended to support 12 aerosol species (original 7 + Pb, Hg, Ni, Cd, As)
to match the modified salsa_mod.f90 with maxspec = 12.
DYNAMIC BIN SUPPORT

Handles all nf2a cases dynamically:
  - nf2a = 1.0: nbin1 + nbin2 bins
  - nf2a = 0.0: nbin1 + nbin2 bins
  - 0.0 < nf2a < 1.0: nbin1 + 2*nbin2 bins

Key scientific references:
- Kurppa et al. (2019) GMD - Main reference for PALM-SALSA implementation
- Kokkola et al. (2008) ACP - Original SALSA description
- Kumar et al. (2009) AE - Traffic emission factors and size distributions
- Zhang et al. (2001) AE - Dry deposition and size distributions
- Dallmann et al. (2014) ACP - Chemical composition of traffic emissions
"""

import os
import glob
import datetime
import re
import time
from datetime import timedelta
import numpy as np
import rasterio
from netCDF4 import Dataset
from rasterio.windows import Window
from scipy.ndimage import zoom
from pyproj import Transformer
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# USER CONFIGURATION OPTIONS
# =============================================================================

# Select which emission categories to include in the output
# Options: 'traffic', 'dust', 'wood' (any combination)
ACTIVE_OUTPUT_CATEGORIES = ['traffic', 'dust', 'wood']  # All three
# ACTIVE_OUTPUT_CATEGORIES = ['traffic', 'wood']        # Exclude road dust
# ACTIVE_OUTPUT_CATEGORIES = ['traffic']                # Only traffic

# Select species for output file
# Option 1: Original 7 species
# SPECIES_OUTPUT_MODE = "basic7"

# Option 2: Extended 12 species (includes Pb, Hg, Ni, Cd, As)
# SPECIES_OUTPUT_MODE = "extended12"  # <-- Change this to "basic7" if needed
# When using CUSTOM_SPECIES_LIST, set this to None or leave commented
SPECIES_OUTPUT_MODE = "custom"  # Add this line - set to None when using custom list

# Or manually specify custom list (uncomment to use):
# CUSTOM_SPECIES_LIST = ["H2SO4", "OC", "BC","DU", "SS", "NH3", "HNO3"]
CUSTOM_SPECIES_LIST = ["H2SO4", "OC", "BC","DU", "SS", "NH3", "HNO3", "PB", "HG", "NI", "CD", "AS"]  # Extended list

# =============================================================================
# CONSTANTS - DYNAMIC BIN CONFIGURATION
# =============================================================================

# PALM/SALSA configuration parameters - CHANGE THESE AS NEEDED!
NBIN = [3, 7]                              # [nbin1, nbin2] - ANY combination works!
REGLIM = [3.9e-8, 1.56e-7, 1.0e-5]        # [d_min, d_split, d_max]
NF2A = 0.75                                # 0.0 to 1.0

# These will be calculated dynamically based on NBIN and NF2A
NBIN1 = NBIN[0]
NBIN2 = NBIN[1]
# Total bins calculated in calculate_palm_bin_diameters()

# Species properties for mass → number conversion
species_properties = {
    "H2SO4": {"rho": 1840, "name": "Sulfuric acid", "molar_mass": 0.098},
    "OC":    {"rho": 1500, "name": "Organic carbon", "molar_mass": 0.012},
    "BC":    {"rho": 1800, "name": "Black carbon", "molar_mass": 0.012},
    "DU":    {"rho": 2650, "name": "Mineral dust", "molar_mass": 0.060},
    "SS":    {"rho": 2200, "name": "Sea salt", "molar_mass": 0.058},
    "HNO3":  {"rho": 1500, "name": "Nitric acid", "molar_mass": 0.063},
    "NH3":   {"rho": 1700, "name": "Ammonia", "molar_mass": 0.017},
    "PB":    {"rho": 11340, "name": "Lead", "molar_mass": 0.2072},
    "HG":    {"rho": 13534, "name": "Mercury", "molar_mass": 0.20059},
    "NI":    {"rho": 8908, "name": "Nickel", "molar_mass": 0.058693},
    "CD":    {"rho": 8650, "name": "Cadmium", "molar_mass": 0.112414},
    "AS":    {"rho": 5727, "name": "Arsenic", "molar_mass": 0.074922},
}

# Species-specific 2a partitioning (only used when 0 < nf2a < 1)
SPECIES_2A_FRACTION = {
    "H2SO4": 0.90, "OC": 0.70, "BC": 0.10, "DU": 0.10,
    "SS": 0.90, "NH3": 0.90, "HNO3": 0.90,
    "PB": 0.50, "HG": 0.50, "NI": 0.50, "CD": 0.50, "AS": 0.50,
}

# Log-normal size distribution parameters
SIZE_DISTRIBUTIONS = {
    0: {"name": "Traffic exhaust", "modes": [
            {"Dg": 20e-9, "sigma": 1.5, "weight": 0.4},
            {"Dg": 60e-9, "sigma": 1.8, "weight": 0.6}]},
    1: {"name": "Road dust", "modes": [
            {"Dg": 2.5e-6, "sigma": 2.0, "weight": 1.0}]},
    2: {"name": "Wood combustion", "modes": [
            {"Dg": 80e-9, "sigma": 1.6, "weight": 1.0}]}
}

# Chemical composition
CATEGORY_COMPOSITION = {
    0: {"H2SO4": 0.04, "OC": 0.48, "BC": 0.48, "DU": 0.00,
        "SS": 0.00, "HNO3": 0.00, "NH3": 0.00,
        "PB": 0.0001, "HG": 0.00001, "NI": 0.0002, "CD": 0.00005, "AS": 0.00004},
    1: {"H2SO4": 0.01, "OC": 0.05, "BC": 0.02, "DU": 0.90,
        "SS": 0.02, "HNO3": 0.00, "NH3": 0.00,
        "PB": 0.0005, "HG": 0.00002, "NI": 0.0003, "CD": 0.0001, "AS": 0.00008},
    2: {"H2SO4": 0.02, "OC": 0.70, "BC": 0.28, "DU": 0.00,
        "SS": 0.00, "HNO3": 0.00, "NH3": 0.00,
        "PB": 0.00005, "HG": 0.00001, "NI": 0.0001, "CD": 0.00003, "AS": 0.0001}
}

SPECIES_CATEGORY_MAPPING = {
    "oc": {"target": "OC", "categories": [0, 2]},
    "bc": {"target": "BC", "categories": [0, 2]},
    "ec": {"target": "BC", "categories": [0, 2]},
    "so2": {"target": "H2SO4", "categories": [0, 2]},
    "no": {"target": "HNO3", "categories": [0, 2]},
    "no2": {"target": "HNO3", "categories": [0, 2]},
    "nh3": {"target": "NH3", "categories": [0, 2]},
    #"pm10": {"target": "DU", "categories": [0, 1, 2]},
    "othmin": {"target": "DU", "categories": [0, 1, 2]},
    "pb": {"target": "PB", "categories": [0, 1]},
    "hg": {"target": "HG", "categories": [0, 1]},
    "ni": {"target": "NI", "categories": [0, 1]},
    "cd": {"target": "CD", "categories": [0, 1]},
    "as": {"target": "AS", "categories": [0, 1]},
    "na": {"target": "SS", "categories": [1]},
}

BULK_SPECIES = ['ho2', 'h2o', 'o3', 'ro2', 'oh', 'rcho', 'n2o', 
                'nmvoc', 'voc', 'co', 'co2', 'ch4', 'pm25', 'pmcoarse']

# Projection
CONFIG_PROJ = "EPSG:25832"
DEFAULT_PROJ = "EPSG:4326"
transformer_to_utm = Transformer.from_crs(DEFAULT_PROJ, CONFIG_PROJ, always_xy=True)
transformer_to_wgs = Transformer.from_crs(CONFIG_PROJ, DEFAULT_PROJ, always_xy=True)


# =============================================================================
# FULLY DYNAMIC BIN CALCULATION - WORKS FOR ANY NBIN!
# =============================================================================

def calculate_palm_bin_diameters(reglim, nbin, nf2a):
    """
    Calculate bin diameters exactly as PALM does.
    FULLY DYNAMIC - works for any nbin = [n1, n2]!
    
    Parameters:
    -----------
    reglim : list
        [d_min, d_split, d_max] - 3 values defining 2 subranges
    nbin : list
        [nbin1, nbin2] - ANY number of bins in subrange 1 and 2
    nf2a : float
        Fraction of subrange 2 allocated to 2a (0.0 to 1.0)
    
    Returns:
    --------
    dmid : array - Bin mid diameters in meters
    dlow : array - Bin lower boundaries
    dhigh : array - Bin upper boundaries
    subrange_labels : array - Labels for each bin
    has_2a : bool - Whether subrange 2a exists
    has_2b : bool - Whether subrange 2b exists
    nbins_total : int - Total number of bins
    """
    d_min, d_split, d_max = reglim
    nbin1, nbin2 = nbin
    
    # Determine which subranges exist based on nf2a
    has_2a = (nf2a > 0.0)
    has_2b = (nf2a < 1.0)
    
    # Calculate total bins dynamically
    if has_2a and has_2b:
        nbins_total = nbin1 + nbin2 + nbin2  # Both 2a and 2b
    elif has_2a:
        nbins_total = nbin1 + nbin2  # Only 2a
    elif has_2b:
        nbins_total = nbin1 + nbin2  # Only 2b
    else:
        nbins_total = nbin1  # Neither (should not happen)
    
    dmid = np.zeros(nbins_total)
    dlow = np.zeros(nbins_total)
    dhigh = np.zeros(nbins_total)
    
    # Build label list dynamically
    labels = ['1'] * nbin1
    if has_2a:
        labels.extend(['2a'] * nbin2)
    if has_2b:
        labels.extend(['2b'] * nbin2)
    subrange_labels = np.array(labels, dtype=str)
    
    # Subrange 1: d_min to d_split (nbin1 bins)
    ratio_1 = d_split / d_min
    for i in range(nbin1):
        d_low = d_min * ratio_1**(i / nbin1)
        d_high = d_min * ratio_1**((i+1) / nbin1)
        dlow[i] = d_low
        dhigh[i] = d_high
        dmid[i] = np.sqrt(d_low * d_high)
    
    # Subrange 2: d_split to d_max (nbin2 bins each for 2a/2b)
    ratio_2 = d_max / d_split
    bin_offset = nbin1
    
    # Subrange 2a (if exists)
    if has_2a:
        for i in range(nbin2):
            idx = bin_offset + i
            d_low = d_split * ratio_2**(i / nbin2)
            d_high = d_split * ratio_2**((i+1) / nbin2)
            dlow[idx] = d_low
            dhigh[idx] = d_high
            dmid[idx] = np.sqrt(d_low * d_high)
        bin_offset += nbin2
    
    # Subrange 2b (if exists)
    if has_2b:
        for i in range(nbin2):
            idx = bin_offset + i
            d_low = d_split * ratio_2**(i / nbin2)
            d_high = d_split * ratio_2**((i+1) / nbin2)
            dlow[idx] = d_low
            dhigh[idx] = d_high
            dmid[idx] = np.sqrt(d_low * d_high)
    
    return dmid, dlow, dhigh, subrange_labels, has_2a, has_2b, nbins_total


def lognormal_pdf(d, dg, sigma_g):
    """Log-normal probability density function"""
    # Add small epsilon to avoid division by zero
    d_safe = np.maximum(d, 1e-30)
    return (1.0 / (d_safe * np.log(sigma_g) * np.sqrt(2 * np.pi))) * \
           np.exp(-(np.log(d_safe / dg)**2) / (2 * np.log(sigma_g)**2))


def get_size_distribution_fractions(category, bin_diameters, subrange_labels, nf2a, has_2a, has_2b):
    """
    Calculate number fractions in each PALM bin.
    FULLY DYNAMIC - works for any bin configuration!
    """
    if category not in SIZE_DISTRIBUTIONS:
        raise ValueError(f"Unknown category: {category}")
    
    dist_info = SIZE_DISTRIBUTIONS[category]
    modes = dist_info["modes"]
    
    # Calculate base PDF at each bin diameter
    pdf_values = np.zeros(len(bin_diameters))
    
    for mode in modes:
        dg = mode["Dg"]
        sigma = mode["sigma"]
        weight = mode["weight"]
        mode_pdf = lognormal_pdf(bin_diameters, dg, sigma)
        pdf_values += weight * mode_pdf
    
    fractions = pdf_values.copy()
    
    # Apply nf2a splitting only if both modes exist
    if has_2a and has_2b:
        mask_2a = subrange_labels == '2a'
        mask_2b = subrange_labels == '2b'
        fractions[mask_2a] = fractions[mask_2a] * nf2a
        fractions[mask_2b] = fractions[mask_2b] * (1.0 - nf2a)
    # If only one mode exists, no splitting needed
    
    # Normalize
    total = np.sum(fractions)
    if total > 0:
        fractions = fractions / total
    
    return fractions


def mass_to_number_conversion(mass_flux, species, bin_diameter, size_fraction):
    """Convert mass flux to number flux for a specific size bin"""
    if species not in species_properties:
        raise ValueError(f"Unknown species: {species}")
    
    rho = species_properties[species]["rho"]
    bin_mass_flux = mass_flux * size_fraction
    volume = (4.0/3.0) * np.pi * (bin_diameter/2.0)**3
    particle_mass = rho * volume
    
    with np.errstate(divide='ignore', invalid='ignore'):
        number_flux = np.where((particle_mass > 0) & (bin_mass_flux > 0), 
                               bin_mass_flux / particle_mass, 
                               0.0)
    
    return number_flux


def extract_hour_from_band_name(band_name):
    if not band_name:
        return None
    hour_match = re.search(r'h(\d{1,2})', band_name)
    if hour_match:
        hour = int(hour_match.group(1))
        if 0 <= hour <= 23:
            return hour
        else:
            return 0 if hour == 24 else None
    return None


def pattern_match(band_name, patterns):
    for pattern in patterns:
        regex_pattern = pattern.replace('*', '.*')
        if re.match(regex_pattern, band_name, re.IGNORECASE):
            return True
    return False


def get_crs_from_netcdf(nc_file):
    try:
        if hasattr(nc_file, 'crs'):
            crs_str = nc_file.crs
            if 'EPSG' in crs_str:
                epsg_code = crs_str.split('EPSG:')[-1].split()[0]
                return f"EPSG:{epsg_code}"
        return CONFIG_PROJ
    except:
        return CONFIG_PROJ


def extract_static_domain(static_nc):
    params = {}
    params['origin_time'] = static_nc.getncattr('origin_time')
    params['origin_lat'] = static_nc.getncattr('origin_lat')
    params['origin_lon'] = static_nc.getncattr('origin_lon')
    
    center_x, center_y = transformer_to_utm.transform(
        params['origin_lon'], params['origin_lat'])
    
    params['nx'] = len(static_nc.dimensions['x'])
    params['ny'] = len(static_nc.dimensions['y'])
    
    x_coords = static_nc.variables['x'][:]
    y_coords = static_nc.variables['y'][:]
    
    params['dx'] = x_coords[1] - x_coords[0] if len(x_coords) > 1 else 1.0
    params['dy'] = abs(y_coords[1] - y_coords[0]) if len(y_coords) > 1 else 1.0
    
    half_nx = (params['nx'] - 1) * params['dx'] / 2
    half_ny = (params['ny'] - 1) * params['dy'] / 2
    
    params['west'] = center_x - half_nx
    params['east'] = center_x + half_nx
    params['south'] = center_y - half_ny
    params['north'] = center_y + half_ny
    
    params['origin_x'] = params['west']
    params['origin_y'] = params['north']
    
    params['lon_w'], params['lat_s'] = transformer_to_wgs.transform(
        params['west'], params['south'])
    params['lon_e'], params['lat_n'] = transformer_to_wgs.transform(
        params['east'], params['north'])
    
    return params


# =============================================================================
# TIFF PROCESSOR CLASS - DYNAMIC BIN SUPPORT
# =============================================================================

class TiffProcessor:
    def __init__(self, static_params, static_crs, active_categories, 
                 nx, ny, ntime, bin_diameters, bin_low, bin_high, subrange_labels,
                 size_distributions, nf2a, has_2a, has_2b, nbins_total):
        self.static_params = static_params
        self.static_crs = static_crs
        self.active_categories = active_categories
        self.nx = nx
        self.ny = ny
        self.ntime = ntime
        self.bin_diameters = bin_diameters
        self.bin_low = bin_low
        self.bin_high = bin_high
        self.subrange_labels = subrange_labels
        self.nbins = nbins_total  # Dynamic!
        self.size_distributions = size_distributions
        self.nf2a = nf2a
        self.has_2a = has_2a
        self.has_2b = has_2b
        
        self.static_xmin = static_params['west']
        self.static_xmax = static_params['east']
        self.static_ymin = static_params['south']
        self.static_ymax = static_params['north']
    
    def get_tiff_extent(self, tiff_file):
        with rasterio.open(tiff_file) as src:
            transform = src.transform
            width = src.width
            height = src.height
            tiff_crs = src.crs.to_string() if src.crs else CONFIG_PROJ
            left = transform[2]
            top = transform[5]
            right = left + transform[0] * width
            bottom = top + transform[4] * height
            return left, bottom, right, top, transform, tiff_crs, width, height

    def convert_extent_to_target_crs(self, left, bottom, right, top, source_crs, target_crs):
        if source_crs == target_crs:
            return left, bottom, right, top
        try:
            transformer = Transformer.from_crs(source_crs, target_crs, always_xy=True)
            left_top = transformer.transform(left, top)
            right_top = transformer.transform(right, top)
            right_bottom = transformer.transform(right, bottom)
            left_bottom = transformer.transform(left, bottom)
            x_coords = [left_top[0], right_top[0], right_bottom[0], left_bottom[0]]
            y_coords = [left_top[1], right_top[1], right_bottom[1], left_bottom[1]]
            new_left = min(x_coords)
            new_right = max(x_coords)
            new_bottom = min(y_coords)
            new_top = max(y_coords)
            return new_left, new_bottom, new_right, new_top
        except Exception as e:
            print(f"Error converting CRS: {e}")
            return left, bottom, right, top

    def crop_tiff_to_static_domain(self, tiff_file, band_idx=1):
        (tiff_left, tiff_bottom, tiff_right, tiff_top, 
         tiff_transform, tiff_crs, tiff_width, tiff_height) = self.get_tiff_extent(tiff_file)
        
        if tiff_crs != self.static_crs:
            tiff_left, tiff_bottom, tiff_right, tiff_top = self.convert_extent_to_target_crs(
                tiff_left, tiff_bottom, tiff_right, tiff_top, tiff_crs, self.static_crs
            )
        
        if (tiff_right < self.static_xmin or tiff_left > self.static_xmax or
            tiff_top < self.static_ymin or tiff_bottom > self.static_ymax):
            return np.zeros((self.ny, self.nx))
        
        overlap_left = max(tiff_left, self.static_xmin)
        overlap_right = min(tiff_right, self.static_xmax)
        overlap_bottom = max(tiff_bottom, self.static_ymin)
        overlap_top = min(tiff_top, self.static_ymax)
        
        pixel_width = tiff_transform[0]
        pixel_height = abs(tiff_transform[4])
        
        col_start = int((overlap_left - tiff_left) / pixel_width)
        row_start = int((tiff_top - overlap_top) / pixel_height)
        cols_to_read = int((overlap_right - overlap_left) / pixel_width)
        rows_to_read = int((overlap_top - overlap_bottom) / pixel_height)
        
        col_start = max(0, min(col_start, tiff_width - 1))
        row_start = max(0, min(row_start, tiff_height - 1))
        cols_to_read = max(1, min(cols_to_read, tiff_width - col_start))
        rows_to_read = max(1, min(rows_to_read, tiff_height - row_start))
        
        with rasterio.open(tiff_file) as src:
            window = Window(col_start, row_start, cols_to_read, rows_to_read)
            arr = src.read(band_idx, window=window)
            
            if arr.shape != (self.ny, self.nx):
                if arr.shape[0] > 0 and arr.shape[1] > 0:
                    zoom_factors = (self.ny / arr.shape[0], self.nx / arr.shape[1])
                    arr = zoom(arr, zoom_factors, order=1)
                else:
                    arr = np.zeros((self.ny, self.nx))
            
            arr = arr / (365.25 * 24 * 3600)
            arr = np.flipud(arr)
            return arr

    def get_band_names(self, tiff_file):
        with rasterio.open(tiff_file) as src:
            return src.descriptions

    def process_single_file(self, tiff_file):
        filename = os.path.basename(tiff_file).lower()
        
        if any(bulk in filename for bulk in BULK_SPECIES):
            return None
        
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
        categories = mapping["categories"]
        
        if not categories:
            return None
        
        print(f"\nProcessing {filename}")
        print(f"  Species: {species} → {target_species}")
        
        band_names = self.get_band_names(tiff_file)
        if not band_names or all(name is None for name in band_names):
            return None
        
        # Dynamic bin count!
        category_emissions = {
            cat: np.zeros((self.ntime, self.ny, self.nx, self.nbins))
            for cat in [0, 1, 2]
        }
        
        f_2a = SPECIES_2A_FRACTION.get(target_species, 0.5)
        
        bands_processed = 0
        for band_idx, band_name in enumerate(band_names, 1):
            if not band_name:
                continue
            if not pattern_match(band_name, self.active_categories):
                continue
            
            hour = extract_hour_from_band_name(band_name)
            if hour is None or hour >= self.ntime:
                hour = (band_idx - 1) % self.ntime
            
            mass_data = self.crop_tiff_to_static_domain(tiff_file, band_idx)
            if np.all(mass_data == 0):
                continue
            
            bands_processed += 1
            
            for cat in categories:
                mass_frac = CATEGORY_COMPOSITION[cat].get(target_species, 0.0)
                if mass_frac <= 0:
                    continue
                
                size_fracs_base = self.size_distributions[cat]
                species_mass = mass_data * mass_frac
                
                for bin_idx in range(self.nbins):
                    base_frac = size_fracs_base[bin_idx]
                    label = self.subrange_labels[bin_idx]
                    
                    if self.has_2a and self.has_2b:
                        if label == '2a':
                            adjusted_frac = base_frac * f_2a
                        elif label == '2b':
                            adjusted_frac = base_frac * (1.0 - f_2a)
                        else:
                            adjusted_frac = base_frac
                    else:
                        adjusted_frac = base_frac
                    
                    if adjusted_frac <= 0:
                        continue
                    
                    number_flux = mass_to_number_conversion(
                        species_mass, target_species,
                        self.bin_diameters[bin_idx], adjusted_frac
                    )
                    category_emissions[cat][hour, :, :, bin_idx] += number_flux
        
        if bands_processed == 0:
            return None
        
        results = {}
        for cat in categories:
            results[f"cat_{cat}"] = category_emissions[cat]
        results['species'] = species
        results['target'] = target_species
        
        return results


# =============================================================================
# MAIN SALSA DRIVER CLASS - FULLY DYNAMIC
# =============================================================================

class SalsaDriver:
    def __init__(self, static_file, tiff_dir, output_file, active_categories=None):
        print("=" * 70)
        print("PALM-SALSA Driver Generator - DYNAMIC BIN VERSION")
        print("=" * 70)
        
        self.active_output_categories = ACTIVE_OUTPUT_CATEGORIES
        self.species_output_mode = SPECIES_OUTPUT_MODE
        self.nbin = NBIN
        self.reglim = REGLIM
        self.nf2a = NF2A
        
        self.category_name_to_idx = {'traffic': 0, 'dust': 1, 'wood': 2}
        
        # Determine species list
        if SPECIES_OUTPUT_MODE == "custom" and 'CUSTOM_SPECIES_LIST' in globals():
            self.composition_name_list = CUSTOM_SPECIES_LIST
        elif SPECIES_OUTPUT_MODE == "basic7":
            self.composition_name_list = ["H2SO4", "OC", "BC", "DU", "SS", "HNO3", "NH3"]
        else:
            self.composition_name_list = ["H2SO4", "OC", "BC", "DU", "SS", "HNO3", "NH3",
                                           "PB", "HG", "NI", "CD", "AS"]
        
        print(f"\nConfiguration (DYNAMIC - works for ANY nbin!):")
        print(f"  nbin = {self.nbin}")
        print(f"  reglim = {self.reglim}")
        print(f"  nf2a = {self.nf2a}")
        
        # Calculate bin structure dynamically
        (self.bin_diameters, self.bin_low, self.bin_high, 
         self.subrange_labels, self.has_2a, self.has_2b, self.nbins_total) = \
            calculate_palm_bin_diameters(self.reglim, self.nbin, self.nf2a)
        
        self.nbin1 = self.nbin[0]
        self.nbin2 = self.nbin[1]
        
        # Summary
        print(f"\nBin structure for nf2a = {self.nf2a}:")
        print(f"  Subrange 1: {self.nbin1} bins (always present)")
        if self.has_2a:
            print(f"  Subrange 2a: {self.nbin2} bins (nf2a = {self.nf2a})")
        if self.has_2b:
            print(f"  Subrange 2b: {self.nbin2} bins (1-nf2a = {1-self.nf2a:.2f})")
        print(f"  TOTAL BINS: {self.nbins_total}")
        
        if self.nf2a == 1.0:
            print(f"  NOTE: nf2a = 1.0 → NO subrange 2b")
        elif self.nf2a == 0.0:
            print(f"  NOTE: nf2a = 0.0 → NO subrange 2a")
        
        print(f"\nSpecies: {len(self.composition_name_list)} species")
        
        self.selected_cat_indices = [self.category_name_to_idx[name] 
                                      for name in self.active_output_categories]
        self.selected_cat_names = [SIZE_DISTRIBUTIONS[idx]['name'] 
                                    for idx in self.selected_cat_indices]
        
        print(f"\nSelected categories: {self.selected_cat_names}")
        
        print("\nOpening static driver...")
        self.static_nc = Dataset(static_file, "r")
        self.output_file = output_file
        self.tiff_dir = tiff_dir
        self.active_categories = active_categories if active_categories else ['*']
        
        self.static_params = extract_static_domain(self.static_nc)
        self.nx = self.static_params['nx']
        self.ny = self.static_params['ny']
        self.ntime = 24
        self.static_x = self.static_nc.variables["x"][:]
        self.static_y = self.static_nc.variables["y"][:]
        self.static_crs = get_crs_from_netcdf(self.static_nc)
        
        print(f"Domain: {self.nx} x {self.ny} grid points")
        
        # Print all bins
        print(f"\nBin details ({self.nbins_total} bins):")
        for i, (d, label) in enumerate(zip(self.bin_diameters, self.subrange_labels)):
            print(f"  Bin {i+1:2d} [{label}]: {d*1e9:.1f} nm "
                  f"[{self.bin_low[i]*1e9:.1f}-{self.bin_high[i]*1e9:.1f}] nm")
        
        # Calculate size distributions
        print("\nCalculating size distributions...")
        self.size_distributions = {}
        for cat in [0, 1, 2]:
            dist = get_size_distribution_fractions(
                cat, self.bin_diameters, self.subrange_labels, 
                self.nf2a, self.has_2a, self.has_2b
            )
            self.size_distributions[cat] = dist
            print(f"  {SIZE_DISTRIBUTIONS[cat]['name']}: sum={np.sum(dist):.6f}")
        
        if self.has_2a and self.has_2b:
            print(f"\nSpecies-specific 2a fractions:")
            for spec in self.composition_name_list:
                if spec in SPECIES_2A_FRACTION:
                    f = SPECIES_2A_FRACTION[spec]
                    print(f"  {spec:<6}: {f*100:.0f}% to 2a, {(1-f)*100:.0f}% to 2b")
        
        print(f"\nCreating output file: {output_file}")
        self.nc_file = Dataset(output_file, "w", format="NETCDF4")
        
        self.define_dimensions()
        self.write_global_attributes()
        self.add_variables()
        self.validate_emissions()
        self.finalize()
    
    def write_global_attributes(self):
        print("\nWriting global attributes...")
        for attr in self.static_nc.ncattrs():
            try:
                setattr(self.nc_file, attr, self.static_nc.getncattr(attr))
            except:
                pass
        
        self.nc_file.creation_date = str(datetime.datetime.now())
        self.nc_file.description = (
            f"Aerosol input (SALSA driver) for PALM. "
            f"nbin={self.nbin}, reglim={self.reglim}, nf2a={self.nf2a}. "
            f"Total bins={self.nbins_total}."
        )
        self.nc_file.title = "PALM input file for SALSA aerosol module"
        self.nc_file.institution = "University of Augsburg"
        self.nc_file.author = "Sathish Kumar Vaithiyanadhan"
        self.nc_file.palm_version = "6.0"
        self.nc_file.composition_name = ', '.join(self.composition_name_list)
        self.nc_file.nbin = f"{self.nbin[0]}, {self.nbin[1]}"
        self.nc_file.reglim = f"{self.reglim[0]:.2e}, {self.reglim[1]:.2e}, {self.reglim[2]:.2e}"
        self.nc_file.nf2a = f"{self.nf2a}"
        self.nc_file.total_bins = f"{self.nbins_total}"
        self.nc_file.nbin1 = f"{self.nbin1}"
        self.nc_file.nbin2 = f"{self.nbin2}"
        self.nc_file.has_2a = f"{self.has_2a}"
        self.nc_file.has_2b = f"{self.has_2b}"

    def define_dimensions(self):
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
        self.nc_file.createDimension("Dmid", self.nbins_total)  # DYNAMIC!

        x = self.nc_file.createVariable("x", "f4", ("x",))
        x[:] = self.static_x
        x.units = "m"

        y = self.nc_file.createVariable("y", "f4", ("y",))
        y[:] = self.static_y
        y.units = "m"

        t = self.nc_file.createVariable("time", "f4", ("time",))
        t[:] = np.arange(0, 24 * 3600, 3600)
        t.units = "s"

        Dmid = self.nc_file.createVariable("Dmid", "f4", ("Dmid",))
        Dmid[:] = self.bin_diameters
        Dmid.units = "m"
        Dmid.long_name = "geometric mean diameter of aerosol size bin"
        
        Dlow = self.nc_file.createVariable("Dlow", "f4", ("Dmid",))
        Dlow[:] = self.bin_low
        Dlow.units = "m"
        
        Dhigh = self.nc_file.createVariable("Dhigh", "f4", ("Dmid",))
        Dhigh[:] = self.bin_high
        Dhigh.units = "m"

        ncat_coord = self.nc_file.createVariable("ncat", "i4", ("ncat",))
        ncat_coord[:] = np.arange(1, self.nncat + 1)
        
        comp_index = self.nc_file.createVariable("composition_index", "i4", ("composition_index",))
        comp_index[:] = np.arange(1, self.ncomposition_index + 1)
        
        max_str = self.nc_file.createVariable("max_string_length", "i4", ("max_string_length",))
        max_str[:] = np.arange(1, self.nmax_string_length + 1)

    def add_variables(self):
        print("\nAdding emission variables...")

        all_cat_names = ["traffic exhaust", "road dust", "wood combustion"]
        filtered_names = [all_cat_names[i] for i in self.selected_cat_indices]
        
        nc_cat_name = self.nc_file.createVariable("emission_category_name", "S1", 
                                                   ("ncat", "max_string_length"))
        for i, name in enumerate(filtered_names):
            chars = list(name.ljust(self.nmax_string_length))
            nc_cat_name[i, :] = np.array(list(chars), dtype="S1")

        nc_cat_idx = self.nc_file.createVariable("emission_category_index", "i4", ("ncat",))
        nc_cat_idx[:] = np.arange(1, self.nncat + 1)

        nc_comp_name = self.nc_file.createVariable("composition_name", "S1",
                                                    ("composition_index", "max_string_length"))
        for i, name in enumerate(self.composition_name_list):
            chars = list(name.ljust(self.nmax_string_length))
            nc_comp_name[i, :] = np.array(list(chars), dtype="S1")

        emission_mass_fracs = np.zeros((self.nncat, self.ncomposition_index))
        for new_cat_idx, old_cat_idx in enumerate(self.selected_cat_indices):
            for comp_idx, comp_name in enumerate(self.composition_name_list):
                emission_mass_fracs[new_cat_idx, comp_idx] = CATEGORY_COMPOSITION[old_cat_idx].get(comp_name, 0.0)
        
        for cat in range(self.nncat):
            total = np.sum(emission_mass_fracs[cat, :])
            if total > 0 and abs(total - 1.0) > 1e-6:
                emission_mass_fracs[cat, :] /= total
        
        nc_mass_fracs = self.nc_file.createVariable("emission_mass_fracs", "f4", 
                                                     ("ncat", "composition_index"), fill_value=-9999.0)
        nc_mass_fracs[:] = emission_mass_fracs
        nc_mass_fracs.units = "1"

        emission_num_fracs = np.zeros((self.nncat, self.nbins_total))
        for new_cat_idx, old_cat_idx in enumerate(self.selected_cat_indices):
            emission_num_fracs[new_cat_idx, :] = self.size_distributions[old_cat_idx]
        
        nc_num_fracs = self.nc_file.createVariable("emission_number_fracs", "f4", 
                                                    ("ncat", "Dmid"), fill_value=-9999.0)
        nc_num_fracs[:] = emission_num_fracs
        nc_num_fracs.units = "1"

        nc_aerosol = self.nc_file.createVariable("aerosol_emission_values", "f4",
                                                  ("time", "y", "x", "ncat"), fill_value=-9999.0)
        nc_aerosol.units = "#/m2/s"
        nc_aerosol.lod = 2

        self.generate_emission_data(nc_aerosol)

    def generate_emission_data(self, nc_var):
        print("\n" + "=" * 70)
        print("PROCESSING EMISSION TIFF FILES")
        print("=" * 70)
        
        tiff_patterns = [os.path.join(self.tiff_dir, "emission_*_temporal.tif")]
        tiff_files = []
        for pattern in tiff_patterns:
            files = glob.glob(pattern)
            if files:
                tiff_files.extend(files)
                print(f"Found {len(files)} files")
        
        tiff_files = list(set(tiff_files))
        print(f"\nTotal TIFF files: {len(tiff_files)}")
        
        if len(tiff_files) == 0:
            print("WARNING: No TIFF files found!")
            nc_var[:] = np.zeros((self.ntime, self.ny, self.nx, self.nncat))
            return
        
        processor = TiffProcessor(
            self.static_params, self.static_crs, self.active_categories,
            self.nx, self.ny, self.ntime,
            self.bin_diameters, self.bin_low, self.bin_high, self.subrange_labels,
            self.size_distributions, self.nf2a, self.has_2a, self.has_2b, self.nbins_total
        )
        
        num_processes = max(1, min(cpu_count(), len(tiff_files)))
        print(f"\nUsing {num_processes} processes")
        
        category_accumulators = {
            0: np.zeros((self.ntime, self.ny, self.nx, self.nbins_total)),
            1: np.zeros((self.ntime, self.ny, self.nx, self.nbins_total)),
            2: np.zeros((self.ntime, self.ny, self.nx, self.nbins_total))
        }
        
        processed_files = 0
        with Pool(processes=num_processes) as pool:
            results = pool.map(processor.process_single_file, tiff_files)
            for result in results:
                if result is not None:
                    processed_files += 1
                    for cat_key, cat_data in result.items():
                        if cat_key.startswith('cat_'):
                            cat_idx = int(cat_key.split('_')[1])
                            category_accumulators[cat_idx] += cat_data
        
        print(f"\nProcessed {processed_files} files")
        
        aerosol_emission_values = np.zeros((self.ntime, self.ny, self.nx, self.nncat))
        for new_cat_idx, old_cat_idx in enumerate(self.selected_cat_indices):
            aerosol_emission_values[:, :, :, new_cat_idx] = np.sum(category_accumulators[old_cat_idx], axis=3)
            total = np.sum(aerosol_emission_values[:, :, :, new_cat_idx])
            print(f"  {self.selected_cat_names[new_cat_idx]}: {total:.2e} #/s")
        
        nc_var[:] = aerosol_emission_values

    def validate_emissions(self):
        print("\n" + "=" * 70)
        print("VALIDATION")
        print("=" * 70)
        with Dataset(self.output_file, "r") as nc:
            emissions = nc.variables["aerosol_emission_values"][:]
            print(f"  Shape: {emissions.shape}")
            print(f"  Total emissions: {np.sum(emissions):.2e} #/s")
            print(f"  Dmid dimension: {nc.dimensions['Dmid'].size} bins")

    def finalize(self):
        print("\nFinalizing...")
        self.static_nc.close()
        self.nc_file.close()
        print(f"Successfully created: {self.output_file}")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    total_start = time.time()
    
    # MODIFY THESE PATHS
    static_file = "/home/vaithisa/palm_model_system-v25.10/JOBS/smallegu/INPUT/smallegu_static"
    tiff_dir = "/home/vaithisa/Downscale_Emissions_simple/downscale/"
    output_file = "/home/vaithisa/palm_model_system-v25.10/JOBS/smallegu/smallegu_salsa"
    
    active_categories = ['A_PublicPower', 'B_Industry', 'C_OtherStationaryComb',
                         'D_Fugitives', 'E_Solvents', 'F_RoadTransport',
                         'G_Shipping', 'H_Aviation', 'I_OffRoad', 'J_Waste',
                         'K_AgriLivestock', 'L_AgriOther']
    
    print("\n" + "=" * 70)
    print("STARTING PALM-SALSA DRIVER - DYNAMIC BIN VERSION")
    print("=" * 70)
    print(f"nbin = {NBIN}  ← CHANGE THIS FOR DIFFERENT BIN COUNTS!")
    print(f"reglim = {REGLIM}")
    print(f"nf2a = {NF2A}")
    
    driver = SalsaDriver(static_file, tiff_dir, output_file, active_categories)
    
    total_time = time.time() - total_start
    print("\n" + "=" * 70)
    print(f"COMPLETED in {timedelta(seconds=int(total_time))}")
    print("=" * 70)