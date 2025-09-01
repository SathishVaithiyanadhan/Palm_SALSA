# Chemistry Driver (LOD2) from GRETA emission inventory for PALM Simulations

This repository provides a modular workflow to generate chemistry drivers for the [PALM modeling system](https://gitlab.palm-model.org/releases/palm_model_system/-/releases) using downscaled GRETA emission inventories. The scripts process spatiotemporal emission data into CF-compliant NetCDF files compatible with PALM's LOD2 chemistry module.

LOD 2: Gridded preprocessed hourly (other temporal intervals will be possible in later versions) emission information that is already temporally disaggregated must be supplied by the user. IMPORTANT: In this mode, the initial date of the simulation has to coincide with the first day for which emission values are available - source: https://palm.muk.uni-hannover.de/trac/wiki/doc/app/chemdesc

---

# Attributes and Dimensions

To setup the attributes and the dimensions for the chemistry driver, based on the latest [Anthropogenic Emissions Model](https://docs.palm-model.com/23.10-rc.1/Guide/LES_Model/Modules/Chemistry/CS_model/#temporal-emission-profiles) documentation. This chemistry driver follows the PIDS for the PALM version 25.04. 

## Key features

The driver integrates gridded emission data (e.g., from the GRETA inventory) with PALM's urban microclimate simulations. It implements:
- **AOI extraction** and grids verification from the input static data (the static driver for the simulation is read, necessary data is extracted.)
- **Multiple sector handling**  from the GRETA emission inventory.
- **Multiple species handling**  from the GRETA emission inventory.
- **Source-category handling** for sector-specific emission modeling.
- **Hourly emissions** as the input (LOD2)
- **Emission Stack height** includes the height of the buildings (LOD2).
- Automatic **time-step** synchronization across species based on the emission input. 
- Properly handling the NAN values.

---

## Input data

The following data is required to create the chemistry driver for the PALM simulation using this tool.

1. Downscaled GRETA Emission inventory
	* Check the repo downscale_emissions_local **(https://git.rz.uni-augsburg.de/vaithisa/downscale_emissions_local.git)** to create your own input data. 

2. Static Data 
	* The static data which you have created using the Geospatial data to describe the topography of the simulation domain. 
    * It is used here to extract the AOI and Grid details for the PALM simulation. 

---

## Usage

1. Configure Paths/Parameters
   - Edit chemistry_driver_config.py:

       * Set emis_geotiff_pth and static_pth 

       * Select the preferred active categories and species

2. Run Main Script

    * **python chemistry_driver_main.py** 

3. Output

    * NetCDF file generated at: **{static_pth}/{static}_chemistry**

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