import os
import glob
import math
import datetime
import re
import numpy as np
import rasterio
from netCDF4 import Dataset
from rasterio.windows import Window
from scipy.ndimage import zoom
from pyproj import Transformer, CRS

# === Species properties for mass → number conversion ===
species_properties = {
    "OC":   {"rho": 1500, "d": 0.1e-6},   # 0.1 µm
    "BC":   {"rho": 1800, "d": 0.05e-6},  # 0.05 µm
    "DU":   {"rho": 2650, "d": 1.0e-6},   # 1 µm
    "SS":   {"rho": 2200, "d": 0.5e-6},   # 0.5 µm
    "NH3":  {"rho": 1700, "d": 0.05e-6},  # 0.05 µm
    "H2SO4":{"rho": 1840, "d": 0.02e-6},  # 20 nm (unused here, set to 0)
    "HNO3": {"rho": 1500, "d": 0.05e-6},  # (unused here, set to 0)
}

# Projection configurations
config_proj = "EPSG:25832"  # UTM Zone 32N
default_proj = "EPSG:4326"  # WGS84

# Coordinate transformers
transformer_to_utm = Transformer.from_crs(default_proj, config_proj, always_xy=True)
transformer_to_wgs = Transformer.from_crs(config_proj, default_proj, always_xy=True)

def mass_to_number(mass_flux, species):
    """
    Convert mass flux [kg/m²/hour] → number flux [#/m²/s]
    """
    props = species_properties[species]
    rho = props["rho"]
    d = props["d"]
    volume = (4.0/3.0) * np.pi * (d/2.0)**3
    mass_flux_si = mass_flux / 3600.0  # convert hour → second
    number_flux = mass_flux_si / (rho * volume)
    return number_flux

def pattern_match(band_name, patterns):
    """
    Check if a band name matches any of the provided patterns.
    Patterns can use wildcard '*' to match any characters.
    """
    for pattern in patterns:
        # Convert wildcard pattern to regex
        regex_pattern = pattern.replace('*', '.*')
        if re.match(regex_pattern, band_name):
            return True
    return False

def extract_hour_from_band_name(band_name):
    """
    Extract hour information from band name.
    Expected format: F_RoadTransport_h01_YYYYMMDD or similar
    Returns hour as integer (0-23) or None if not found.
    """
    if not band_name:
        return None
    
    # Try to find hour pattern: h01, h02, ..., h24
    hour_match = re.search(r'h(\d{1,2})', band_name)
    if hour_match:
        hour = int(hour_match.group(1))
        # Convert 1-24 to 0-23 format
        return hour - 1 if hour >= 1 else 0
    return None

def get_crs_from_netcdf(nc_file):
    """
    Extract CRS information from a PALM NetCDF file.
    PALM typically stores CRS as a global attribute.
    """
    try:
        # Try to get CRS from global attributes
        if hasattr(nc_file, 'crs'):
            crs_str = nc_file.crs
            # Check if it's an EPSG code
            if 'EPSG' in crs_str:
                epsg_code = crs_str.split('EPSG:')[-1].split()[0]
                return f"EPSG:{epsg_code}"
        # Default to UTM Zone 32N if not specified
        return config_proj
    except:
        return config_proj  # Default fallback

def extract_static_domain(static_nc):
    """
    Extract domain parameters from static driver file using the same method
    as in the chemistry driver code
    """
    params = {}
    
    # Extract basic attributes
    params['origin_time'] = static_nc.getncattr('origin_time')
    params['origin_lat'] = static_nc.getncattr('origin_lat')
    params['origin_lon'] = static_nc.getncattr('origin_lon')
    
    # Convert center point to UTM coordinates
    center_x, center_y = transformer_to_utm.transform(
        params['origin_lon'], params['origin_lat'])
    
    # Get grid dimensions and resolution
    params['nx'] = len(static_nc.dimensions['x'])
    params['ny'] = len(static_nc.dimensions['y'])
    
    x_coords = static_nc.variables['x'][:]
    y_coords = static_nc.variables['y'][:]
    
    params['dx'] = x_coords[1] - x_coords[0] if len(x_coords) > 1 else 1.0
    params['dy'] = abs(y_coords[1] - y_coords[0]) if len(y_coords) > 1 else 1.0
    
    # Calculate domain boundaries relative to center point
    half_nx = (params['nx'] - 1) * params['dx'] / 2
    half_ny = (params['ny'] - 1) * params['dy'] / 2
    
    params['west'] = center_x - half_nx
    params['east'] = center_x + half_nx
    params['south'] = center_y - half_ny
    params['north'] = center_y + half_ny
    
    # Update origin coordinates to match PALM convention
    params['origin_x'] = params['west']
    params['origin_y'] = params['north']
    
    # Convert boundaries back to WGS84 for reference
    params['lon_w'], params['lat_s'] = transformer_to_wgs.transform(
        params['west'], params['south'])
    params['lon_e'], params['lat_n'] = transformer_to_wgs.transform(
        params['east'], params['north'])
        
    return params

class SalsaDriver:
    """Generate a SALSA driver NetCDF file for PALM using GeoTIFF emissions."""

    def __init__(self, static_file, tiff_dir, output_file="urban_environment_salsa_salsa.nc", 
                 active_categories=None):
        print("Opening static driver...")
        self.static_nc = Dataset(static_file, "r")
        print("Creating salsa driver file...")
        self.nc_file = Dataset(output_file, "w", format="NETCDF4")
        self.tiff_dir = tiff_dir
        
        # Set active categories to process (default to all if not specified)
        self.active_categories = active_categories if active_categories else ['*']
        
        # Extract domain parameters using the chemistry driver method
        self.static_params = extract_static_domain(self.static_nc)
        
        # Get domain dimensions
        self.nx = self.static_params['nx']
        self.ny = self.static_params['ny']
        
        # Get coordinates from static file
        self.static_x = self.static_nc.variables["x"][:]
        self.static_y = self.static_nc.variables["y"][:]
        
        # Get CRS from static file
        self.static_crs = get_crs_from_netcdf(self.static_nc)
        print(f"Static file CRS: {self.static_crs}")
        
        # Get the extent of the static domain
        self.static_xmin = self.static_params['west']
        self.static_xmax = self.static_params['east']
        self.static_ymin = self.static_params['south']
        self.static_ymax = self.static_params['north']
        
        print(f"Static domain extent (native CRS): x={self.static_xmin}:{self.static_xmax}, y={self.static_ymin}:{self.static_ymax}")

    def write_global_attributes(self):
        print("Writing global attributes...")
        for attr in self.static_nc.ncattrs():
            setattr(self.nc_file, attr, self.static_nc.getncattr(attr))
        
        # Set specific global attributes as requested
        self.nc_file.creation_date = str(datetime.datetime.now())
        self.nc_file.description = "Aerosol input (SALSA driver) for PALM to simulate the aerosol particle concentrations, size distributions and chemical compositions."
        self.nc_file.title = "PALM input file for SALSA aerosol module"
        self.nc_file.institution = "Chair of Model-based Environmental Exposure Science, University of Augsburg"
        self.nc_file.author = "Sathish Kumar Vaithiyanadhan"
        self.nc_file.palm_version = "6.0"
        self.nc_file.active_categories = ', '.join(self.active_categories)

    def define_dimensions(self):
        print("Defining dimensions...")

        self.ntime = 24
        self.nncat = 3
        self.ncomposition_index = 7
        self.nmax_string_length = 25

        # === Create dimensions ===
        self.nc_file.createDimension("x", self.nx)
        self.nc_file.createDimension("y", self.ny)
        self.nc_file.createDimension("time", self.ntime)
        self.nc_file.createDimension("ncat", self.nncat)
        self.nc_file.createDimension("composition_index", self.ncomposition_index)
        self.nc_file.createDimension("max_string_length", self.nmax_string_length)
        self.nc_file.createDimension("Dmid", 8)

        # === Coordinates ===
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
        Dmid[:] = np.linspace(0.01e-6, 2.5e-6, 8)  # dummy bins
        Dmid.units = "m"

        # === Dimension coordinate variables ===
        # Create coordinate variables for categorical dimensions
        ncat_coord = self.nc_file.createVariable("ncat", "i4", ("ncat",))
        ncat_coord[:] = np.arange(1, self.nncat + 1)  # 1, 2, 3
        ncat_coord.units = ""
        ncat_coord.long_name = "emission category index"
        
        comp_index_coord = self.nc_file.createVariable("composition_index", "i4", ("composition_index",))
        comp_index_coord[:] = np.arange(1, self.ncomposition_index + 1)  # 1, 2, 3, 4, 5, 6, 7
        comp_index_coord.units = ""
        comp_index_coord.long_name = "composition index"
        
        max_str_len_coord = self.nc_file.createVariable("max_string_length", "i4", ("max_string_length",))
        max_str_len_coord[:] = np.arange(1, self.nmax_string_length + 1)  # 1, 2, 3, ..., 25
        max_str_len_coord.units = ""
        max_str_len_coord.long_name = "maximum string length"

    def add_variables(self):
        print("Adding emission variables...")
        print(f"Active categories: {self.active_categories}")

        # === Emission categories ===
        emission_category_name_list = ["traffic exhaust", "road dust", "wood combustion"]
        nc_emission_category_name = self.nc_file.createVariable(
            "emission_category_name", "S1", ("ncat", "max_string_length")
        )
        for i, name in enumerate(emission_category_name_list):
            chars = list(name.ljust(self.nmax_string_length))
            nc_emission_category_name[i, :] = np.array(list(chars), dtype="S1")
        nc_emission_category_name.long_name = "emission category name"

        # Emission category index (as a data variable, not coordinate)
        nc_emission_category_index = self.nc_file.createVariable(
            "emission_category_index", "i1", ("ncat",)
        )
        nc_emission_category_index[:] = np.arange(1, self.nncat + 1)  # 1-indexed
        nc_emission_category_index.long_name = "emission category index"
        nc_emission_category_index.units = ""

        # === Composition names ===
        composition_name_list = ["H2SO4", "OC", "BC", "DU", "SS", "HNO3", "NH3"]
        nc_composition_name = self.nc_file.createVariable(
            "composition_name", "S1",
            ("composition_index", "max_string_length")
        )
        for i, name in enumerate(composition_name_list):
            chars = list(name.ljust(self.nmax_string_length))
            nc_composition_name[i, :] = np.array(list(chars), dtype="S1")
        nc_composition_name.long_name = "aerosol composition name"

        # === Emission mass fractions ===
        # Define default mass fractions (these should be customized based on your data)
        emission_mass_fracs = np.zeros((self.nncat, self.ncomposition_index))
        # Traffic exhaust: mostly OC, BC, some NH3
        emission_mass_fracs[0, :] = [0.0, 0.5, 0.3, 0.0, 0.0, 0.0, 0.2]  # H2SO4, OC, BC, DU, SS, HNO3, NH3
        # Road dust: mostly DU
        emission_mass_fracs[1, :] = [0.0, 0.1, 0.0, 0.9, 0.0, 0.0, 0.0]
        # Wood combustion: mostly OC, BC
        emission_mass_fracs[2, :] = [0.0, 0.7, 0.3, 0.0, 0.0, 0.0, 0.0]
        
        nc_emission_mass_fracs = self.nc_file.createVariable(
            "emission_mass_fracs", "f4", ("ncat", "composition_index"), fill_value=-9999.0
        )
        nc_emission_mass_fracs[:] = emission_mass_fracs
        nc_emission_mass_fracs.long_name = "mass fractions of chemical components in aerosol emissions"
        nc_emission_mass_fracs.units = "1"
        nc_emission_mass_fracs.coordinates = "ncat composition_index"

        # === Emission number fractions ===
        # Define default size distribution (these should be customized based on your data)
        emission_number_fracs = np.zeros((self.nncat, 8))
        # All categories use a log-normal distribution centered around 0.1 µm
        for i in range(self.nncat):
            emission_number_fracs[i, :] = [0.1, 0.2, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02]
            emission_number_fracs[i, :] /= np.sum(emission_number_fracs[i, :])  # Normalize to sum=1
        
        nc_emission_number_fracs = self.nc_file.createVariable(
            "emission_number_fracs", "f4", ("ncat", "Dmid"), fill_value=-9999.0
        )
        nc_emission_number_fracs[:] = emission_number_fracs
        nc_emission_number_fracs.long_name = "number fractions of aerosol size bins in aerosol emissions"
        nc_emission_number_fracs.units = "1"
        nc_emission_number_fracs.coordinates = "ncat Dmid"

        # === Aerosol emissions ===
        nc_aerosol_emission_values = self.nc_file.createVariable(
            "aerosol_emission_values", "f4",
            ("time", "y", "x", "ncat"), fill_value=-9999.0
        )
        nc_aerosol_emission_values.units = "#/m2/s"
        nc_aerosol_emission_values.long_name = "aerosol emission values"
        nc_aerosol_emission_values.source = "Based on Kumar et al., 2009: Comparative study of measured and modelled number concentrations of nanoparticles in an urban street canyon"
        nc_aerosol_emission_values.lod = 2
        nc_aerosol_emission_values.coordinates = "time y x ncat"

        # --- Initialize arrays for each hour and category ---
        # Create separate arrays for each emission category and each hour
        traffic_by_hour = np.zeros((self.ntime, self.ny, self.nx))
        road_dust_by_hour = np.zeros((self.ntime, self.ny, self.nx))
        wood_by_hour = np.zeros((self.ntime, self.ny, self.nx))
        
        # Define which species to process and how to map them to categories
        species_mapping = {
            "oc": {"target": "OC", "category": [0, 2]},  # OC goes to traffic (0) and wood (2)
            "bc": {"target": "BC", "category": [0, 2]},  # BC goes to traffic (0) and wood (2)
            "na": {"target": "SS", "category": []},      # SS not assigned to any category
            "nh3": {"target": "NH3", "category": [0]},   # NH3 goes to traffic (0)
            "pb": {"target": "DU", "category": [1]},     # DU goes to road dust (1)
            "cd": {"target": "DU", "category": [1]},     # DU goes to road dust (1)
            "hg": {"target": "DU", "category": [1]},     # DU goes to road dust (1)
            "as": {"target": "DU", "category": [1]},     # DU goes to road dust (1)
            "ni": {"target": "DU", "category": [1]},     # DU goes to road dust (1)
            "othmin": {"target": "DU", "category": [1]}  # DU goes to road dust (1)
        }
        
        # Process each TIFF file
        processed_files = 0
        for tiff_file in glob.glob(os.path.join(self.tiff_dir, "emission_*_temporal.tif")):
            species = os.path.basename(tiff_file).split("_")[1].lower()
            
            # Skip species that are not in our mapping
            if species not in species_mapping:
                print(f"Skipping {species}, not mapped to any aerosol species.")
                continue
                
            print(f"Processing {species} from {tiff_file}")
            processed_files += 1
            
            # Get all band names from the TIFF file
            band_names = self.get_band_names(tiff_file)
            if not band_names or all(name is None for name in band_names):
                print(f"Warning: No band names found in {tiff_file}, skipping")
                continue
            
            # Process each band that matches active categories and has valid hour info
            target_species = species_mapping[species]["target"]
            target_categories = species_mapping[species]["category"]
            
            if not target_categories:
                print(f"Warning: Species {species} not assigned to any emission category, skipping")
                continue
                
            for band_idx, band_name in enumerate(band_names, 1):
                if not band_name:
                    continue
                    
                # Check if band matches active categories
                if not pattern_match(band_name, self.active_categories):
                    continue
                
                # Extract hour from band name
                hour = extract_hour_from_band_name(band_name)
                if hour is None or hour >= self.ntime:
                    print(f"Warning: Could not extract valid hour from band '{band_name}', skipping")
                    continue
                
                print(f"  Processing hour {hour} from band {band_idx}: {band_name}")
                
                # Crop the TIFF to match the static domain for this band
                arr = self.crop_tiff_to_static_domain(tiff_file, band_idx)
                converted_arr = mass_to_number(arr, target_species)
                
                # Fix orientation
                converted_arr = np.flipud(converted_arr)
                
                # Add to the appropriate emission categories
                for category in target_categories:
                    if category == 0:  # Traffic exhaust
                        traffic_by_hour[hour, :, :] += converted_arr
                    elif category == 1:  # Road dust
                        road_dust_by_hour[hour, :, :] += converted_arr
                    elif category == 2:  # Wood combustion
                        wood_by_hour[hour, :, :] += converted_arr

        print(f"Processed {processed_files} TIFF files with emission data")
        
        # Check if we have any data
        total_traffic = np.sum(traffic_by_hour)
        total_road_dust = np.sum(road_dust_by_hour)
        total_wood = np.sum(wood_by_hour)
        
        print(f"Total emissions by category:")
        print(f"  Traffic exhaust: {total_traffic:.2e} #/m²/s")
        print(f"  Road dust: {total_road_dust:.2e} #/m²/s")
        print(f"  Wood combustion: {total_wood:.2e} #/m²/s")
        
        if total_traffic == 0 and total_road_dust == 0 and total_wood == 0:
            print("WARNING: All emission values are zero! Check your input data and active categories.")

        # --- Assign to final output array ---
        aerosol_emission_values = np.zeros((self.ntime, self.ny, self.nx, self.nncat))

        for t in range(self.ntime):
            aerosol_emission_values[t, :, :, 0] = traffic_by_hour[t, :, :]  # Traffic exhaust
            aerosol_emission_values[t, :, :, 1] = road_dust_by_hour[t, :, :]  # Road dust
            aerosol_emission_values[t, :, :, 2] = wood_by_hour[t, :, :]  # Wood combustion

        nc_aerosol_emission_values[:] = aerosol_emission_values

    def get_tiff_extent(self, tiff_file):
        """Get the spatial extent of a TIFF file and its CRS"""
        with rasterio.open(tiff_file) as src:
            transform = src.transform
            width = src.width
            height = src.height
            tiff_crs = src.crs.to_string() if src.crs else config_proj
            
            # Calculate corner coordinates in the TIFF's native CRS
            left = transform[2]
            top = transform[5]
            right = left + transform[0] * width
            bottom = top + transform[4] * height
            
            return left, bottom, right, top, transform, tiff_crs

    def convert_extent_to_target_crs(self, left, bottom, right, top, source_crs, target_crs):
        """Convert extent coordinates from source CRS to target CRS"""
        if source_crs == target_crs:
            return left, bottom, right, top
            
        try:
            transformer = Transformer.from_crs(source_crs, target_crs, always_xy=True)
            
            # Convert all four corners
            left_top = transformer.transform(left, top)
            right_top = transformer.transform(right, top)
            right_bottom = transformer.transform(right, bottom)
            left_bottom = transformer.transform(left, bottom)
            
            # Find the min/max coordinates in the target CRS
            x_coords = [left_top[0], right_top[0], right_bottom[0], left_bottom[0]]
            y_coords = [left_top[1], right_top[1], right_bottom[1], left_bottom[1]]
            
            new_left = min(x_coords)
            new_right = max(x_coords)
            new_bottom = min(y_coords)
            new_top = max(y_coords)
            
            return new_left, new_bottom, new_right, new_top
            
        except Exception as e:
            print(f"Error converting CRS from {source_crs} to {target_crs}: {e}")
            return left, bottom, right, top  # Return original if conversion fails

    def crop_tiff_to_static_domain(self, tiff_file, band_idx=1):
        """Crop TIFF data to match the static domain extremet for a specific band"""
        # Get TIFF extent and transform
        tiff_left, tiff_bottom, tiff_right, tiff_top, tiff_transform, tiff_crs = self.get_tiff_extent(tiff_file)
        
        # Convert TIFF extent to static CRS if needed
        if tiff_crs != self.static_crs:
            tiff_left, tiff_bottom, tiff_right, tiff_top = self.convert_extent_to_target_crs(
                tiff_left, tiff_bottom, tiff_right, tiff_top, tiff_crs, self.static_crs
            )
        
        # Check if there's any overlap
        if (tiff_right < self.static_xmin or tiff_left > self.static_xmax or
            tiff_top < self.static_ymin or tiff_bottom > self.static_ymax):
            print(f"Warning: No overlap between TIFF {tiff_file} and static domain")
            return np.zeros((self.ny, self.nx))
        
        # Calculate the overlapping region
        overlap_left = max(tiff_left, self.static_xmin)
        overlap_right = min(tiff_right, self.static_xmax)
        overlap_bottom = max(tiff_bottom, self.static_ymin)
        overlap_top = min(tiff_top, self.static_ymax)
        
        # Calculate pixel size (assuming regular grid)
        pixel_width = tiff_transform[0]
        pixel_height = abs(tiff_transform[4])
        
        # Calculate row and column offsets for cropping
        col_start = int((overlap_left - tiff_left) / pixel_width)
        row_start = int((tiff_top - overlap_top) / pixel_height)
        
        # Calculate number of columns and rows to read
        cols_to_read = int((overlap_right - overlap_left) / pixel_width)
        rows_to_read = int((overlap_top - overlap_bottom) / pixel_height)
        
        # Ensure we don't exceed TIFF boundaries
        col_start = max(0, col_start)
        row_start = max(0, row_start)
        
        with rasterio.open(tiff_file) as src:
            # Read only the portion that overlaps with the static domain for the specified band
            window = Window(col_start, row_start, cols_to_read, rows_to_read)
            arr = src.read(band_idx, window=window)
            
            # Resize to match static domain dimensions if needed
            if arr.shape != (self.ny, self.nx):
                if arr.shape[0] > 0 and arr.shape[1] > 0:  # Ensure we have data to resize
                    zoom_factors = (self.ny / arr.shape[0], self.nx / arr.shape[1])
                    arr = zoom(arr, zoom_factors, order=1)  # linear interpolation
                else:
                    print(f"Warning: Empty array after cropping for {tiff_file}")
                    arr = np.zeros((self.ny, self.nx))
            
            return arr

    def get_band_names(self, tiff_file):
        """Get all band names from a TIFF file"""
        with rasterio.open(tiff_file) as src:
            return src.descriptions

    def finalize(self):
        print("Closing files...")
        self.static_nc.close()
        self.nc_file.close()


if __name__ == "__main__":
    static_file = "/home/vaithisa/GEO4PALM-main/JOBS/Augsburg_10/OUTPUT/Augsburg_3_static"
    tiff_dir = "/home/vaithisa/Downscale_Emissions/Downscale_Winter_10m"
    
    # Define which categories to process
    active_categories = [
        'A_PublicPower', 
        'B_Industry', 
        'C_OtherStationaryComb', 
        #'D_Fugitives',
        #'E_Solvents', 
        'F_RoadTransport', 
        #'G_Shipping', 
        #'H_Aviation',
        #'I_OffRoad', 
        #'J_Waste', 
        #'K_AgriLivestock', 
        #'L_AgriOther',
        #'SumAllSectors'
    ]
    
    driver = SalsaDriver(static_file, tiff_dir, "/home/vaithisa/Palm_SALSA/Augsburg_32_salsa", active_categories)
    driver.write_global_attributes()
    driver.define_dimensions()
    driver.add_variables()
    driver.finalize()