# SALSA Aerosol module (LOD2) from GRETA emission inventory for PALM Simulations

This repository provides a modular workflow to generate aerosol input for the [PALM modeling system](https://gitlab.palm-model.org/releases/palm_model_system/-/releases) using downscaled GRETA emission inventories. The scripts process spatiotemporal emission data into CF-compliant NetCDF files compatible with PALM's LOD2 aerosol module.

LOD 2: Gridded preprocessed hourly (other temporal intervals will be possible in later versions) emission information that is already temporally disaggregated must be supplied by the user. IMPORTANT: In this mode, the initial date of the simulation has to coincide with the first day for which emission values are available - source: https://palm.muk.uni-hannover.de/trac/wiki/doc/app/chemdesc

---

# Attributes and Dimensions

To setup the attributes and the dimensions for the salsa driver, based on the latest [SALSA Aerosol module](https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/aerosol) documentation. This salsa aerosol driver follows the PIDS for the PALM version 25.04. 

# Dimensions

- **x**: Spatial dimension in the x-direction (meters, from static driver).
- **y**: Spatial dimension in the y-direction (meters, from static driver).
- **time**: Temporal dimension (24 hours, in seconds).
- **ncat**: Number of emission categories (3: traffic exhaust, road dust, wood combustion).
- **composition_index**: Number of aerosol chemical components (7: H2SO4, OC, BC, DU, SS, HNO3, NH3).
- **max_string_length**: Length of string arrays for category and composition names (25 characters).
- **Dmid**: Number of aerosol size bins (8 bins, from 0.01 µm to 2.5 µm).

# Variables

- **x (f4, x)**: Distance to origin in x-direction (m).
- **y (f4, y)**: Distance to origin in y-direction (m).
- **time (f4, time)**: Time in seconds (hourly, 0 to 86,400 s).
- **ncat (i4, ncat)**: Emission category indices (1, 2, 3).
- **max_string_length (i4, max_string_length)**: Character positions for string arrays (1 to 25).
- **composition_index (i4, composition_index)**: Aerosol composition indices (1 to 7).
- **emission_category_name (S1, ncat, max_string_length)**: Names of emission categories (e.g., "traffic exhaust", "road dust", "wood combustion").
- **emission_category_index (i4, ncat)**: Indices of emission categories (1, 2, 3).
- **composition_name (S1, composition_index, max_string_length)**: Names of aerosol components (e.g., "H2SO4", "OC", "BC").
- **emission_mass_fracs (f4, ncat, composition_index)**: Mass fractions of chemical components for each emission category (dimensionless).
- **emission_number_fracs (f4, ncat, Dmid)**: Number fractions of aerosol size bins for each emission category (dimensionless).
- **aerosol_emission_values (f4, time, y, x, ncat)**: Aerosol emission values in number flux (#/m²/s).



## Key features

The SALSA aerosol driver integrates gridded emission data from the GRETA inventory with PALM’s urban microclimate simulations, providing the following features:

- **Area of Interest (AOI) Extraction**: Extracts grid dimensions and spatial extent from the PALM static driver file.
- **Multiple Sector Handling**: Processes emission data from multiple sectors (e.g., F_RoadTransport, C_OtherStationaryComb) in the GRETA inventory.
- **Multiple Species Handling**: Maps GRETA species (e.g., oc, bc, nh3, pb) to SALSA aerosol components (e.g., H₂SO₄, OC, BC, DU (dust & metals), SS (sea salt / Na), HNO₃, NH₃ ).
- **Source-Category Handling**: Supports three emission categories (traffic exhaust, road dust, wood combustion) with weighted contributions for overlapping species (e.g., OC and BC split between traffic and wood combustion).
- **Hourly Emissions**: Processes 24 hourly bands from GeoTIFF files, ensuring temporal alignment with PALM’s LOD2 requirements.
- **Automatic Time-Step Synchronization**: Aligns emission data across species using hourly GeoTIFF bands.
- **NaN Handling**: Uses a fill value of -9999.0 for missing data in numerical arrays.


---

## Input data

The following data is required to create the chemistry driver for the PALM simulation using this tool.

1. Downscaled GRETA Emission inventory
	* Check the repo downscale_emissions_local **(https://git.rz.uni-augsburg.de/vaithisa/downscale_emissions_local.git)** to create your own input data. 
    * Produces hourly, gridded emissions as GeoTIFFs (`emission_{species}_temporal.tif`)

2. Static Data 
	* The static data which you have created using the Geospatial data to describe the topography of the simulation domain. 
    * It is used here to extract the AOI and Grid details for the PALM simulation. 

---

## Usage

1. Configure Paths and Parameters:

   - Edit the __main__ block in palm_salsa_driver.py to set:

        * static_file: Path to the PALM static driver file (e.g., /home/vaithisa/GEO4PALM-main/JOBS/Augsburg_10/OUTPUT/Augsburg_3_static).
        * tiff_dir: Directory containing downscaled GRETA GeoTIFF files (e.g., /home/vaithisa/Downscale_Emissions/Downscale_Winter_10m).
        * output_file (driver): Path for the output NetCDF file (e.g., /home/vaithisa/Palm_SALSA/Augsburg_3_salsa).
        * active_categories: List of emission sectors to process (e.g., ['F_RoadTransport', 'C_OtherStationaryComb', etc.,]).

2. Run Main Script

    * **python palm_salsa_driver.py** 

3. Output

    * NetCDF file generated in the defined **output_file (driver) path**

---

## Authors and acknowledgment

Show your appreciation to those who have contributed to the project.
For details and comments, please contact:
1. Sathish Kumar Vaithiyanadhan (sathish.vaithiyanadhan@uni-a.de)
2. Christoph Knote (christoph.knote@med.uni-augsburg.de)

@ Chair of Model-based Environmental Exposure Science (MBEES), Faculty of Medicine, University of Augsburg, Germany.

---

## License

For open source projects, say how it is licensed.

---