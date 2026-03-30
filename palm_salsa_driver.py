###
#12 species
#!/usr/bin/env python3
"""
PALM SALSA Driver Generator
============================
Generate aerosol emission input files for PALM's SALSA module using the methodology
described in Kurppa et al. (2019) "Implementation of the sectional aerosol module 
SALSA2.0 into the PALM model system 6.0"

MODIFIED: Extended to support 12 aerosol species (original 7 + Pb, Hg, Ni, Cd, As)
to match the modified salsa_mod.f90 with maxspec = 12.

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
# CONSTANTS AND SCIENTIFIC PARAMETERS
# =============================================================================

# Species properties for mass → number conversion (from Kokkola et al., 2008)
# EXTENDED: Added Pb, Hg, Ni, Cd, As
species_properties = {
    "H2SO4": {"rho": 1840, "name": "Sulfuric acid", "molar_mass": 0.098},   # kg/mol
    "OC":    {"rho": 1500, "name": "Organic carbon", "molar_mass": 0.012},  # kg/mol (approx)
    "BC":    {"rho": 1800, "name": "Black carbon", "molar_mass": 0.012},    # kg/mol
    "DU":    {"rho": 2650, "name": "Mineral dust", "molar_mass": 0.060},    # kg/mol (approx)
    "SS":    {"rho": 2200, "name": "Sea salt", "molar_mass": 0.058},        # kg/mol (NaCl)
    "HNO3":  {"rho": 1500, "name": "Nitric acid", "molar_mass": 0.063},     # kg/mol
    "NH3":   {"rho": 1700, "name": "Ammonia", "molar_mass": 0.017},         # kg/mol
    # NEW TRACE METAL SPECIES (added to match salsa_mod.f90)
    "PB":    {"rho": 11340, "name": "Lead", "molar_mass": 0.2072},          # kg/mol (Pb)
    "HG":    {"rho": 13534, "name": "Mercury", "molar_mass": 0.20059},      # kg/mol (Hg)
    "NI":    {"rho": 8908, "name": "Nickel", "molar_mass": 0.058693},       # kg/mol (Ni)
    "CD":    {"rho": 8650, "name": "Cadmium", "molar_mass": 0.112414},      # kg/mol (Cd)
    "AS":    {"rho": 5727, "name": "Arsenic", "molar_mass": 0.074922},      # kg/mol (As)
}

# PALM/SALSA configuration parameters (from your namelist)
REGLIM = [3.9e-8, 5.0e-8, 2.5e-6]  # [dmin1, dmin2, dmax2] in meters
NBIN = [1, 7]                       # [nbin1, nbin2] - 8 total bins
NBINS_TOTAL = NBIN[0] + NBIN[1]      # 8 bins

# Log-normal size distribution parameters for different emission categories
# Based on Kumar et al. (2009), Zhang et al. (2001), and literature
SIZE_DISTRIBUTIONS = {
    0: {  # Traffic exhaust (nucleation + accumulation mode)
        "name": "Traffic exhaust",
        "modes": [
            {"Dg": 20e-9, "sigma": 1.5, "weight": 0.4},   # Nucleation mode
            {"Dg": 60e-9, "sigma": 1.8, "weight": 0.6},   # Accumulation mode
        ],
        "reference": "Kumar et al. (2009)"
    },
    1: {  # Road dust (coarse mode)
        "name": "Road dust",
        "modes": [
            {"Dg": 2.5e-6, "sigma": 2.0, "weight": 1.0},  # Coarse mode
        ],
        "reference": "Zhang et al. (2001)"
    },
    2: {  # Wood combustion (accumulation mode)
        "name": "Wood combustion",
        "modes": [
            {"Dg": 80e-9, "sigma": 1.6, "weight": 1.0},   # Accumulation mode
        ],
        "reference": "Hays et al. (2002)"
    }
}

# Chemical composition for each emission category (mass fractions)
# EXTENDED: Added Pb, Hg, Ni, Cd, As for all categories
# Based on Kurppa et al. (2019), Dallmann et al. (2014), and literature
CATEGORY_COMPOSITION = {
    0: {  # Traffic exhaust
        "H2SO4": 0.04,   # 4% - sulfuric acid
        "OC":    0.48,   # 48% - organic carbon
        "BC":    0.48,   # 48% - black carbon
        "DU":    0.00,   # 0% - mineral dust
        "SS":    0.00,   # 0% - sea salt
        "HNO3":  0.00,   # 0% - nitric acid (gas phase, condenses later)
        "NH3":   0.00,   # 0% - ammonia (gas phase, condenses later)
        # NEW TRACE METALS (trace amounts from fuel and lubricants)
        "PB":    0.0001, # 0.01% - Lead (from legacy fuel, brake wear)
        "HG":    0.00001,# 0.001% - Mercury (trace from fuel)
        "NI":    0.0002, # 0.02% - Nickel (from fuel, lubricants)
        "CD":    0.00005,# 0.005% - Cadmium (trace from tires)
        "AS":    0.00004,# 0.004% - Arsenic (trace from fuel)
        "reference": "Dallmann et al. (2014), Kurppa et al. (2019), US EPA AP-42"
    },
    1: {  # Road dust (mineral dust)
        "H2SO4": 0.01,   # 1% - sulfuric acid
        "OC":    0.05,   # 5% - organic carbon (from road debris)
        "BC":    0.02,   # 2% - black carbon (from tire wear)
        "DU":    0.90,   # 90% - mineral dust
        "SS":    0.02,   # 2% - sea salt (road salt in winter)
        "HNO3":  0.00,   # 0% - nitric acid
        "NH3":   0.00,   # 0% - ammonia
        # NEW TRACE METALS (from road wear, brake wear, tire wear)
        "PB":    0.0005, # 0.05% - Lead (from brake pads, legacy road paint)
        "HG":    0.00002,# 0.002% - Mercury (trace)
        "NI":    0.0003, # 0.03% - Nickel (from tire wear, road surface)
        "CD":    0.0001, # 0.01% - Cadmium (from tire wear)
        "AS":    0.00008,# 0.008% - Arsenic (trace)
        "reference": "US EPA AP-42, Zhang et al. (2001), Thorpe & Harrison (2008)"
    },
    2: {  # Wood combustion
        "H2SO4": 0.02,   # 2% - sulfuric acid
        "OC":    0.70,   # 70% - organic carbon
        "BC":    0.28,   # 28% - black carbon
        "DU":    0.00,   # 0% - mineral dust
        "SS":    0.00,   # 0% - sea salt
        "HNO3":  0.00,   # 0% - nitric acid
        "NH3":   0.00,   # 0% - ammonia
        # NEW TRACE METALS (from wood itself and treated wood)
        "PB":    0.00005,# 0.005% - Lead (trace from soil/wood)
        "HG":    0.00001,# 0.001% - Mercury (trace)
        "NI":    0.0001, # 0.01% - Nickel (trace from wood)
        "CD":    0.00003,# 0.003% - Cadmium (trace)
        "AS":    0.0001, # 0.01% - Arsenic (from treated wood)
        "reference": "McDonald et al. (2000), Fine et al. (2001), Naeher et al. (2007)"
    }
}

# Species to category mapping - using ONLY speciated species (NO double counting)
# EXTENDED: Added mapping for Pb, Hg, Ni, Cd, As
SPECIES_CATEGORY_MAPPING = {
    # Traffic-related species
    "oc":      {"target": "OC",    "categories": [0, 2]},     # OC from traffic and wood
    "bc":      {"target": "BC",    "categories": [0, 2]},     # BC from traffic and wood
    "ec":      {"target": "BC",    "categories": [0, 2]},     # Elemental carbon = BC
    "so2":     {"target": "H2SO4", "categories": [0, 2]},     # SO2 → H2SO4 (traffic, wood)
    "no":      {"target": "HNO3",  "categories": [0, 2]},     # NO → HNO3
    "no2":     {"target": "HNO3",  "categories": [0, 2]},     # NO2 → HNO3
    "nh3":     {"target": "NH3",   "categories": [0, 2]},     # NH3 (traffic, wood)
    
    # Road dust related (heavy metals and minerals)
    "pm10":    {"target": "DU",    "categories": [0, 1, 2]},  # PM10 (treated as dust)
    "pm2_5":   {"target": "DU",    "categories": [0, 1, 2]},  # PM2.5 (treated as dust)
    
    # NEW TRACE METAL SPECIES
    "pb":      {"target": "PB",    "categories": [0, 1, 2]},  # Lead - all categories
    "hg":      {"target": "HG",    "categories": [0, 1, 2]},  # Mercury - all categories
    "ni":      {"target": "NI",    "categories": [0, 1, 2]},  # Nickel - all categories
    "cd":      {"target": "CD",    "categories": [0, 1, 2]},  # Cadmium - all categories
    "as":      {"target": "AS",    "categories": [0, 1, 2]},  # Arsenic - all categories
    
    # Sea salt (if applicable)
    "na":      {"target": "SS",    "categories": [1]},        # Not used in inland
}

# Bulk species to SKIP (avoid double counting)
BULK_SPECIES = ['ho2', 'h2o', 'o3', 'ro2', 'oh', 'rcho', 'n2o', 'nmvoc', 'voc', 'co', 'co2', 'ch4']

# Projection configurations
CONFIG_PROJ = "EPSG:25832"  # UTM Zone 32N (Augsburg region)
DEFAULT_PROJ = "EPSG:4326"  # WGS84

# Coordinate transformers
transformer_to_utm = Transformer.from_crs(DEFAULT_PROJ, CONFIG_PROJ, always_xy=True)
transformer_to_wgs = Transformer.from_crs(CONFIG_PROJ, DEFAULT_PROJ, always_xy=True)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def calculate_palm_bin_diameters(reglim, nbin):
    """
    Calculate bin diameters exactly as PALM does (Kokkola et al., 2008)
    
    Parameters:
    -----------
    reglim : list
        [dmin1, dmin2, dmax2] - size range limits in meters
    nbin : list
        [nbin1, nbin2] - number of bins in each subrange
    
    Returns:
    --------
    dmid : array
        Bin mid diameters in meters
    dlow : array
        Bin lower boundaries
    dhigh : array
        Bin upper boundaries
    """
    dmin1, dmin2, dmax2 = reglim
    nbin1, nbin2 = nbin
    
    nbins_total = nbin1 + nbin2
    dmid = np.zeros(nbins_total)
    dlow = np.zeros(nbins_total)
    dhigh = np.zeros(nbins_total)
    
    # Subrange 1: logarithmic spacing
    ratio_d1 = dmin2 / dmin1
    for i in range(nbin1):
        d_low = dmin1 * ratio_d1**(i / nbin1)
        d_high = dmin1 * ratio_d1**((i+1) / nbin1)
        dlow[i] = d_low
        dhigh[i] = d_high
        dmid[i] = np.sqrt(d_low * d_high)
    
    # Subrange 2: logarithmic spacing  
    ratio_d2 = dmax2 / dmin2
    for i in range(nbin2):
        idx = nbin1 + i
        d_low = dmin2 * ratio_d2**(i / nbin2)
        d_high = dmin2 * ratio_d2**((i+1) / nbin2)
        dlow[idx] = d_low
        dhigh[idx] = d_high
        dmid[idx] = np.sqrt(d_low * d_high)
    
    return dmid, dlow, dhigh


def lognormal_pdf(d, dg, sigma_g):
    """
    Log-normal probability density function
    
    Parameters:
    -----------
    d : array
        Diameters to evaluate PDF at [m]
    dg : float
        Geometric mean diameter [m]
    sigma_g : float
        Geometric standard deviation
    
    Returns:
    --------
    pdf : array
        Probability density values
    """
    return (1.0 / (d * np.log(sigma_g) * np.sqrt(2 * np.pi))) * \
           np.exp(-(np.log(d / dg)**2) / (2 * np.log(sigma_g)**2))


def get_size_distribution_fractions(category, bin_diameters):
    """
    Calculate number fractions in each size bin for a given emission category
    
    Uses multi-modal log-normal distributions based on literature values.
    
    Parameters:
    -----------
    category : int
        Emission category (0=traffic, 1=road dust, 2=wood)
    bin_diameters : array
        Bin mid diameters in meters
    
    Returns:
    --------
    fractions : array
        Number fractions for each bin (normalized to sum to 1)
    """
    if category not in SIZE_DISTRIBUTIONS:
        raise ValueError(f"Unknown category: {category}")
    
    dist_info = SIZE_DISTRIBUTIONS[category]
    modes = dist_info["modes"]
    
    # Calculate PDF at each bin diameter for each mode
    pdf_values = np.zeros(len(bin_diameters))
    
    for mode in modes:
        dg = mode["Dg"]
        sigma = mode["sigma"]
        weight = mode["weight"]
        
        mode_pdf = lognormal_pdf(bin_diameters, dg, sigma)
        pdf_values += weight * mode_pdf
    
    # Normalize to sum to 1
    fractions = pdf_values / np.sum(pdf_values)
    
    return fractions


def mass_to_number_conversion(mass_flux, species, bin_diameter, size_fraction):
    """
    Convert mass flux to number flux for a specific size bin
    
    Parameters:
    -----------
    mass_flux : float or array
        Total mass flux [kg/m²/s] for this species
    species : str
        Chemical species name (must be in species_properties)
    bin_diameter : float
        Diameter of the size bin [m]
    size_fraction : float
        Fraction of total mass in this bin (from size distribution)
    
    Returns:
    --------
    number_flux : float or array
        Number flux in this bin [#/m²/s]
    """
    if species not in species_properties:
        raise ValueError(f"Unknown species: {species}")
    
    rho = species_properties[species]["rho"]
    
    # Mass in this bin
    bin_mass_flux = mass_flux * size_fraction
    
    # Volume of a single spherical particle in this bin
    volume = (4.0/3.0) * np.pi * (bin_diameter/2.0)**3
    
    # Mass of a single particle
    particle_mass = rho * volume
    
    # Number of particles (avoid division by zero)
    with np.errstate(divide='ignore', invalid='ignore'):
        number_flux = np.where((particle_mass > 0) & (bin_mass_flux > 0), 
                               bin_mass_flux / particle_mass, 
                               0.0)
    
    return number_flux


def extract_hour_from_band_name(band_name):
    """
    Extract hour information from band name.
    Expected format: F_RoadTransport_h00_YYYYMMDD -> 0
                     F_RoadTransport_h23_YYYYMMDD -> 23
    
    Returns:
    --------
    hour : int (0-23) or None if not found
    """
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
    """
    Check if a band name matches any of the provided patterns.
    Patterns can use wildcard '*' to match any characters.
    """
    for pattern in patterns:
        regex_pattern = pattern.replace('*', '.*')
        if re.match(regex_pattern, band_name, re.IGNORECASE):
            return True
    return False


def get_crs_from_netcdf(nc_file):
    """
    Extract CRS information from a PALM NetCDF file.
    """
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
    """
    Extract domain parameters from static driver file
    """
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
# TIFF PROCESSOR CLASS (for parallel processing)
# =============================================================================

class TiffProcessor:
    """
    Helper class to process TIFF files in parallel with proper aerosol physics
    """
    
    def __init__(self, static_params, static_crs, active_categories, 
                 nx, ny, ntime, bin_diameters, size_distributions):
        """
        Initialize the processor
        
        Parameters:
        -----------
        static_params : dict
            Domain parameters from static file
        static_crs : str
            CRS of the static domain
        active_categories : list
            List of active emission category patterns
        nx, ny, ntime : int
            Grid dimensions
        bin_diameters : array
            Bin mid diameters in meters
        size_distributions : dict
            Pre-calculated size distributions for each category
        """
        self.static_params = static_params
        self.static_crs = static_crs
        self.active_categories = active_categories
        self.nx = nx
        self.ny = ny
        self.ntime = ntime
        self.bin_diameters = bin_diameters
        self.nbins = len(bin_diameters)
        
        self.size_distributions = size_distributions
        
        self.static_xmin = static_params['west']
        self.static_xmax = static_params['east']
        self.static_ymin = static_params['south']
        self.static_ymax = static_params['north']
    
    def get_tiff_extent(self, tiff_file):
        """Get the spatial extent of a TIFF file and its CRS"""
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
        """Convert extent coordinates from source CRS to target CRS"""
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
            print(f"Error converting CRS from {source_crs} to {target_crs}: {e}")
            return left, bottom, right, top

    def crop_tiff_to_static_domain(self, tiff_file, band_idx=1):
        """
        Crop TIFF data to match the static domain extent for a specific band
        """
        (tiff_left, tiff_bottom, tiff_right, tiff_top, 
         tiff_transform, tiff_crs, tiff_width, tiff_height) = self.get_tiff_extent(tiff_file)
        
        if tiff_crs != self.static_crs:
            tiff_left, tiff_bottom, tiff_right, tiff_top = self.convert_extent_to_target_crs(
                tiff_left, tiff_bottom, tiff_right, tiff_top, tiff_crs, self.static_crs
            )
        
        if (tiff_right < self.static_xmin or tiff_left > self.static_xmax or
            tiff_top < self.static_ymin or tiff_bottom > self.static_ymax):
            print(f"  Warning: No overlap between TIFF and static domain")
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
                    print(f"  Warning: Empty array after cropping")
                    arr = np.zeros((self.ny, self.nx))
            
            arr = arr / (365.25 * 24 * 3600)
            arr = np.flipud(arr) 
            
            return arr

    def get_band_names(self, tiff_file):
        """Get all band names from a TIFF file"""
        with rasterio.open(tiff_file) as src:
            return src.descriptions

    def process_single_file(self, tiff_file):
        """
        Process a single TIFF file with proper aerosol physics
        
        This method implements the scientific approach:
        1. Skip bulk species to avoid double counting
        2. Read mass emissions for each species
        3. Map species to target chemical component and emission categories
        4. Distribute mass across size bins according to log-normal distributions
        5. Convert mass to number concentration per bin
        6. Aggregate by category
        """
        filename = os.path.basename(tiff_file).lower()
        
        if any(bulk in filename for bulk in BULK_SPECIES):
            print(f"Skipping {filename}: bulk species")
            return None
        
        species = None
        for key in SPECIES_CATEGORY_MAPPING.keys():
            if key in filename:
                species = key
                break
        
        if species is None:
            print(f"Skipping {filename}: No matching species found")
            return None
        
        mapping = SPECIES_CATEGORY_MAPPING.get(species)
        if not mapping:
            return None
        
        target_species = mapping["target"]
        categories = mapping["categories"]
        
        if not categories:
            print(f"Species {species} not assigned to any emission category, skipping")
            return None
        
        print(f"\nProcessing {filename}")
        print(f"  Species: {species} → {target_species}")
        print(f"  Categories: {categories}")
        
        band_names = self.get_band_names(tiff_file)
        if not band_names or all(name is None for name in band_names):
            print(f"  Warning: No band names found, skipping")
            return None
        
        # Initialize arrays for each category and bin
        # Shape: (ntime, ny, nx, nbins)
        category_emissions = {
            cat: np.zeros((self.ntime, self.ny, self.nx, self.nbins))
            for cat in [0, 1, 2]
        }
        
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
            if bands_processed <= 3:
                print(f"  Processing hour {hour:2d} from band {band_idx}")
            
            for cat in categories:
                mass_frac = CATEGORY_COMPOSITION[cat].get(target_species, 0.0)
                
                if mass_frac <= 0:
                    continue
                
                size_fracs = self.size_distributions[cat]
                species_mass = mass_data * mass_frac
                
                for bin_idx, size_frac in enumerate(size_fracs):
                    if size_frac <= 0:
                        continue
                    
                    number_flux = mass_to_number_conversion(
                        species_mass,
                        target_species,
                        self.bin_diameters[bin_idx],
                        size_frac
                    )
                    
                    category_emissions[cat][hour, :, :, bin_idx] += number_flux
        
        if bands_processed == 0:
            print(f"  Warning: No valid bands processed for {filename}")
            return None
        
        results = {}
        for cat in categories:
            results[f"cat_{cat}"] = category_emissions[cat]
        
        results['species'] = species
        results['target'] = target_species
        results['file'] = filename
        results['bands_processed'] = bands_processed
        
        return results


# =============================================================================
# MAIN SALSA DRIVER CLASS
# =============================================================================

class SalsaDriver:
    """
    Generate a SALSA driver NetCDF file for PALM using GeoTIFF emissions.
    
    This class implements the methodology described in Kurppa et al. (2019)
    for creating aerosol emission input files for PALM's SALSA module.
    
    MODIFIED: Extended to support 12 species (original 7 + Pb, Hg, Ni, Cd, As)
    """
    
    def __init__(self, static_file, tiff_dir, output_file, active_categories=None):
        """
        Initialize the SALSA driver generator
        
        Parameters:
        -----------
        static_file : str
            Path to PALM static file
        tiff_dir : str
            Directory containing emission TIFF files
        output_file : str
            Path for output NetCDF file
        active_categories : list, optional
            List of active emission category patterns
        """
        print("=" * 70)
        print("PALM-SALSA Driver Generator (12 Species Version)")
        print("=" * 70)
        
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
        print(f"Static file CRS: {self.static_crs}")
        
        self.static_xmin = self.static_params['west']
        self.static_xmax = self.static_params['east']
        self.static_ymin = self.static_params['south']
        self.static_ymax = self.static_params['north']
        
        print(f"Domain dimensions: {self.nx} x {self.ny} grid points")
        print(f"Domain extent: x=[{self.static_xmin:.1f}, {self.static_xmax:.1f}] m")
        print(f"               y=[{self.static_ymin:.1f}, {self.static_ymax:.1f}] m")
        
        self.bin_diameters, self.bin_low, self.bin_high = calculate_palm_bin_diameters(REGLIM, NBIN)
        print(f"\nSALSA size bins ({NBINS_TOTAL} total):")
        for i, d in enumerate(self.bin_diameters):
            print(f"  Bin {i+1}: {d:.2e} m ({d*1e9:.1f} nm) - "
                  f"range [{self.bin_low[i]*1e9:.1f}, {self.bin_high[i]*1e9:.1f}] nm")
        
        print("\nCalculating size distributions for each category:")
        self.size_distributions = {}
        for cat in [0, 1, 2]:
            dist = get_size_distribution_fractions(cat, self.bin_diameters)
            self.size_distributions[cat] = dist
            print(f"  Category {cat} ({SIZE_DISTRIBUTIONS[cat]['name']}):")
            print(f"    Bin fractions: " + " ".join([f"{f:.3f}" for f in dist]))
        
        print(f"\nCreating output file: {output_file}")
        self.nc_file = Dataset(output_file, "w", format="NETCDF4")
        
        self.define_dimensions()
        self.write_global_attributes()
        
        self.add_variables()
        
        self.validate_emissions()
        
        self.finalize()
    
    def write_global_attributes(self):
        """Write global attributes to NetCDF file"""
        print("\nWriting global attributes...")
        
        for attr in self.static_nc.ncattrs():
            try:
                setattr(self.nc_file, attr, self.static_nc.getncattr(attr))
            except:
                pass
        
        self.nc_file.creation_date = str(datetime.datetime.now())
        self.nc_file.description = (
            "Aerosol input (SALSA driver) for PALM to simulate aerosol particle "
            "concentrations, size distributions and chemical compositions. "
            "EXTENDED: Includes 12 species (H2SO4, OC, BC, DU, SS, HNO3, NH3, "
            "PB, HG, NI, CD, AS) to match modified salsa_mod.f90 with maxspec=12."
        )
        self.nc_file.title = "PALM input file for SALSA aerosol module (12 species)"
        self.nc_file.institution = "Chair of Model-based Environmental Exposure Science, University of Augsburg"
        self.nc_file.lod = 2
        self.nc_file.author = "Sathish Kumar Vaithiyanadhan"
        self.nc_file.palm_version = "6.0"
        self.nc_file.active_categories = ', '.join(self.active_categories)
        
        # UPDATED: Include all 12 species
        self.nc_file.composition_name = "H2SO4, OC, BC, DU, SS, HNO3, NH3, PB, HG, NI, CD, AS"
        self.nc_file.reglim = f"{REGLIM[0]:.2e}, {REGLIM[1]:.2e}, {REGLIM[2]:.2e}"
        self.nc_file.nbin = f"{NBIN[0]}, {NBIN[1]}"
        
        for cat in [0, 1, 2]:
            comp = CATEGORY_COMPOSITION[cat]
            attr_name = f"category_{cat}_composition"
            comp_items = [(k, v) for k, v in comp.items() if k != 'reference']
            comp_str = ", ".join([f"{k}={v:.5f}" for k, v in comp_items])
            setattr(self.nc_file, attr_name, comp_str)

    def define_dimensions(self):
        """Define NetCDF dimensions"""
        print("Defining dimensions...")
        
        # UPDATED: ncomposition_index = 12 (for 12 species)
        self.nncat = 3                     # Number of emission categories
        self.ncomposition_index = 12      # Number of chemical components (increased from 7 to 12)
        self.nmax_string_length = 25
        
        self.nc_file.createDimension("x", self.nx)
        self.nc_file.createDimension("y", self.ny)
        self.nc_file.createDimension("time", self.ntime)
        self.nc_file.createDimension("ncat", self.nncat)
        self.nc_file.createDimension("composition_index", self.ncomposition_index)
        self.nc_file.createDimension("max_string_length", self.nmax_string_length)
        self.nc_file.createDimension("Dmid", NBINS_TOTAL)

        x = self.nc_file.createVariable("x", "f4", ("x",))
        x[:] = self.static_x
        x.units = "m"
        x.long_name = "distance to origin in x-direction"

        y = self.nc_file.createVariable("y", "f4", ("y",))
        y[:] = self.static_y
        y.units = "m"
        y.long_name = "distance to origin in y-direction"

        t = self.nc_file.createVariable("time", "f4", ("time",))
        t[:] = np.arange(0, 24 * 3600, 3600)
        t.units = "s"
        t.long_name = "time in seconds"

        Dmid = self.nc_file.createVariable("Dmid", "f4", ("Dmid",))
        Dmid[:] = self.bin_diameters
        Dmid.units = "m"
        Dmid.long_name = "geometric mean diameter of aerosol size bin"
        
        Dlow = self.nc_file.createVariable("Dlow", "f4", ("Dmid",))
        Dlow[:] = self.bin_low
        Dlow.units = "m"
        Dlow.long_name = "lower diameter of aerosol size bin"
        
        Dhigh = self.nc_file.createVariable("Dhigh", "f4", ("Dmid",))
        Dhigh[:] = self.bin_high
        Dhigh.units = "m"
        Dhigh.long_name = "upper diameter of aerosol size bin"

        ncat_coord = self.nc_file.createVariable("ncat", "i4", ("ncat",))
        ncat_coord[:] = np.arange(1, self.nncat + 1)
        ncat_coord.units = ""
        ncat_coord.long_name = "emission category index"
        
        comp_index_coord = self.nc_file.createVariable("composition_index", "i4", 
                                                        ("composition_index",))
        comp_index_coord[:] = np.arange(1, self.ncomposition_index + 1)
        comp_index_coord.units = ""
        comp_index_coord.long_name = "composition index"
        
        max_str_len_coord = self.nc_file.createVariable("max_string_length", "i4", 
                                                         ("max_string_length",))
        max_str_len_coord[:] = np.arange(1, self.nmax_string_length + 1)
        max_str_len_coord.units = ""
        max_str_len_coord.long_name = "maximum string length"

    def add_variables(self):
        """Add emission variables to NetCDF file"""
        print("\nAdding emission variables...")
        print(f"Active categories: {self.active_categories}")

        emission_category_name_list = [
            "traffic exhaust",
            "road dust", 
            "wood combustion"
        ]
        nc_emission_category_name = self.nc_file.createVariable(
            "emission_category_name", "S1", ("ncat", "max_string_length")
        )
        for i, name in enumerate(emission_category_name_list):
            chars = list(name.ljust(self.nmax_string_length))
            nc_emission_category_name[i, :] = np.array(list(chars), dtype="S1")
        nc_emission_category_name.long_name = "emission category name"

        nc_emission_category_index = self.nc_file.createVariable(
            "emission_category_index", "i1", ("ncat",)
        )
        nc_emission_category_index[:] = np.arange(1, self.nncat + 1)
        nc_emission_category_index.long_name = "emission category index"
        nc_emission_category_index.units = ""

        # UPDATED: Composition names list extended to 12 species
        composition_name_list = ["H2SO4", "OC", "BC", "DU", "SS", "HNO3", "NH3", 
                                 "PB", "HG", "NI", "CD", "AS"]
        nc_composition_name = self.nc_file.createVariable(
            "composition_name", "S1",
            ("composition_index", "max_string_length")
        )
        for i, name in enumerate(composition_name_list):
            chars = list(name.ljust(self.nmax_string_length))
            nc_composition_name[i, :] = np.array(list(chars), dtype="S1")
        nc_composition_name.long_name = "aerosol composition name"

        # UPDATED: emission_mass_fracs array now has 12 columns
        emission_mass_fracs = np.zeros((self.nncat, self.ncomposition_index))
        for cat in [0, 1, 2]:
            for comp_idx, comp_name in enumerate(composition_name_list):
                emission_mass_fracs[cat, comp_idx] = CATEGORY_COMPOSITION[cat].get(comp_name, 0.0)
        
        # Normalize to ensure sum = 1 (safety check)
        for cat in range(self.nncat):
            total = np.sum(emission_mass_fracs[cat, :])
            if total > 0 and abs(total - 1.0) > 1e-6:
                emission_mass_fracs[cat, :] /= total
        
        nc_emission_mass_fracs = self.nc_file.createVariable(
            "emission_mass_fracs", "f4", ("ncat", "composition_index"), 
            fill_value=-9999.0
        )
        nc_emission_mass_fracs[:] = emission_mass_fracs
        nc_emission_mass_fracs.long_name = "mass fractions of chemical components in aerosol emissions"
        nc_emission_mass_fracs.units = "1"
        nc_emission_mass_fracs.coordinates = "ncat composition_index"

        emission_number_fracs = np.zeros((self.nncat, NBINS_TOTAL))
        for cat in [0, 1, 2]:
            emission_number_fracs[cat, :] = self.size_distributions[cat]
        
        nc_emission_number_fracs = self.nc_file.createVariable(
            "emission_number_fracs", "f4", ("ncat", "Dmid"), 
            fill_value=-9999.0
        )
        nc_emission_number_fracs[:] = emission_number_fracs
        nc_emission_number_fracs.long_name = "number fractions of aerosol size bins in aerosol emissions"
        nc_emission_number_fracs.units = "1"
        nc_emission_number_fracs.coordinates = "ncat Dmid"

        nc_aerosol_emission_values = self.nc_file.createVariable(
            "aerosol_emission_values", "f4",
            ("time", "y", "x", "ncat"), fill_value=-9999.0
        )
        nc_aerosol_emission_values.units = "#/m2/s"
        nc_aerosol_emission_values.long_name = (
            "aerosol emission values (total number concentration across all size bins)"
        )
        nc_aerosol_emission_values.source = (
            "Based on Kumar et al. (2009) for traffic, Zhang et al. (2001) for road dust, "
            "and literature values for wood combustion. Extended to include trace metals "
            "Pb, Hg, Ni, Cd, As based on US EPA AP-42 and Thorpe & Harrison (2008)."
        )
        nc_aerosol_emission_values.lod = 2
        nc_aerosol_emission_values.coordinates = "time y x ncat"

        self.generate_emission_data(nc_aerosol_emission_values)

    def generate_emission_data(self, nc_var):
        """
        Generate emission data by processing all TIFF files
        """
        print("\n" + "=" * 70)
        print("PROCESSING EMISSION TIFF FILES")
        print("=" * 70)
        
        tiff_patterns = [
            os.path.join(self.tiff_dir, "emission_*_temporal.tif"),
        ]
        
        tiff_files = []
        for pattern in tiff_patterns:
            files = glob.glob(pattern)
            if files:
                tiff_files.extend(files)
                print(f"Found {len(files)} files with pattern: {pattern}")
        
        tiff_files = list(set(tiff_files))
        print(f"\nTotal unique TIFF files to process: {len(tiff_files)}")
        
        if len(tiff_files) == 0:
            print("WARNING: No TIFF files found! Creating empty emission data.")
            aerosol_emission_values = np.zeros((self.ntime, self.ny, self.nx, self.nncat))
            nc_var[:] = aerosol_emission_values
            return
        
        processor = TiffProcessor(
            self.static_params,
            self.static_crs,
            self.active_categories,
            self.nx,
            self.ny,
            self.ntime,
            self.bin_diameters,
            self.size_distributions
        )
        
        num_processes = max(1, min(cpu_count(), len(tiff_files)))
        print(f"\nUsing {num_processes} processes for parallel processing")
        
        category_accumulators = {
            0: np.zeros((self.ntime, self.ny, self.nx, NBINS_TOTAL)),
            1: np.zeros((self.ntime, self.ny, self.nx, NBINS_TOTAL)),
            2: np.zeros((self.ntime, self.ny, self.nx, NBINS_TOTAL))
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
        
        print(f"\nProcessed {processed_files} TIFF files successfully")
        print(f"Skipped {len(tiff_files) - processed_files} files")
        
        aerosol_emission_values = np.zeros((self.ntime, self.ny, self.nx, self.nncat))
        
        for cat in [0, 1, 2]:
            aerosol_emission_values[:, :, :, cat] = np.sum(category_accumulators[cat], axis=3)
            
            total_emissions = np.sum(aerosol_emission_values[:, :, :, cat])
            max_emissions = np.max(aerosol_emission_values[:, :, :, cat])
            nonzero_mask = aerosol_emission_values[:, :, :, cat] > 0
            mean_nonzero = np.mean(aerosol_emission_values[nonzero_mask, cat]) if np.any(nonzero_mask) else 0
            
            print(f"\nCategory {cat} ({SIZE_DISTRIBUTIONS[cat]['name']}):")
            print(f"  Total emissions: {total_emissions:.2e} #/s")
            print(f"  Max grid cell: {max_emissions:.2e} #/m²/s")
            print(f"  Mean (non-zero): {mean_nonzero:.2e} #/m²/s")
        
        nc_var[:] = aerosol_emission_values

    def validate_emissions(self):
        """
        Validate emission values against literature expectations
        """
        print("\n" + "=" * 70)
        print("VALIDATION OF EMISSION VALUES")
        print("=" * 70)
        
        with Dataset(self.output_file, "r") as nc:
            emissions = nc.variables["aerosol_emission_values"][:]
            
            total_by_category = np.sum(emissions, axis=(0, 1, 2))
            max_by_category = np.max(emissions, axis=(0, 1, 2))
            nonzero_mask = emissions > 0
            mean_by_category = np.array([
                np.mean(emissions[nonzero_mask[:, :, :, cat], cat]) 
                if np.any(nonzero_mask[:, :, :, cat]) else 0
                for cat in range(3)
            ])
            
            print("\nEmission Statistics:")
            print("-" * 50)
            
            for cat in [0, 1, 2]:
                cat_name = SIZE_DISTRIBUTIONS[cat]['name']
                print(f"\nCategory {cat}: {cat_name}")
                print(f"  Total emissions: {total_by_category[cat]:.2e} #/s over domain")
                print(f"  Max grid cell: {max_by_category[cat]:.2e} #/m²/s")
                if mean_by_category[cat] > 0:
                    print(f"  Mean (non-zero cells): {mean_by_category[cat]:.2e} #/m²/s")
            
            print("\nValidation against literature:")
            
            traffic_emissions = emissions[:, :, :, 0]
            traffic_nonzero = traffic_emissions[traffic_emissions > 0]
            
            if len(traffic_nonzero) > 0:
                avg_traffic_per_cell = np.mean(traffic_nonzero)
                expected_range = [1e8, 1e9]
                
                print(f"\nTraffic emissions:")
                print(f"  Modeled mean: {avg_traffic_per_cell:.2e} #/m²/s")
                
                if expected_range[0] <= avg_traffic_per_cell <= expected_range[1]:
                    print("  Traffic emissions within expected range")
                elif avg_traffic_per_cell < expected_range[0]:
                    print("  Traffic emissions below expected range (may be underestimation)")
                else:
                    print("  Traffic emissions above expected range (may need checking)")
            
            zero_cells = np.sum(emissions == 0, axis=(0, 1, 2))
            total_cells = emissions.shape[0] * emissions.shape[1] * emissions.shape[2]
            
            print(f"\nZero value check:")
            for cat in [0, 1, 2]:
                cat_name = SIZE_DISTRIBUTIONS[cat]['name']
                zero_pct = 100 * zero_cells[cat] / total_cells
                print(f"  {cat_name}: {zero_pct:.1f}% zero cells")
            
            hourly_totals = np.sum(emissions, axis=(1, 2, 3))
            print(f"\nDiurnal variation:")
            print(f"  Min hour: {np.argmin(hourly_totals)} (total: {np.min(hourly_totals):.2e})")
            print(f"  Max hour: {np.argmax(hourly_totals)} (total: {np.max(hourly_totals):.2e})")

    def finalize(self):
        """Close NetCDF files"""
        print("\nFinalizing and closing files...")
        self.static_nc.close()
        self.nc_file.close()
        print(f"Successfully created: {self.output_file}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    static_file = "/home/vaithisa/palm_model_system-v25.10/JOBS/Customized/INPUT/Customized_static"
    tiff_dir = "/home/vaithisa/Downscale_Emissions_simple/downscale/"
    output_file = "/home/vaithisa/palm_model_system-v25.10/JOBS/Customized/Customized_salsa"
    
    active_categories = [
        'A_PublicPower',
        'B_Industry',
        'C_OtherStationaryComb',
        'D_Fugitives',
        'E_Solvents',
        'F_RoadTransport',
        'G_Shipping',
        'H_Aviation',
        'I_OffRoad',
        'J_Waste',
        'K_AgriLivestock',
        'L_AgriOther',
    ]
    
    print("\n" + "=" * 70)
    print("STARTING PALM-SALSA DRIVER GENERATION (12 SPECIES)")
    print("=" * 70)
    print(f"Static file: {static_file}")
    print(f"TIFF directory: {tiff_dir}")
    print(f"Output file: {output_file}")
    print(f"Active categories: {len(active_categories)} patterns")
    
    driver = SalsaDriver(
        static_file=static_file,
        tiff_dir=tiff_dir,
        output_file=output_file,
        active_categories=active_categories
    )
    
    print("\n" + "=" * 70)
    print("PROCESS COMPLETED SUCCESSFULLY")
    print("=" * 70)