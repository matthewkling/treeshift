## Data sources

### Forest data

US Forest Service Forest Inventory and Analysis database (FIADB), downloaded June 2025. TreeShift uses the most recent visit to each plot, on the standard four-subplot design. Plots with evidence of planting are excluded, so the data describe naturally established forests. "Non-forest" FIA plot locations that were not visited in the field are used only where LANDFIRE Biophysical Settings indicate non-forest vegetation, so they can stand in for places without trees. This tool uses public FIA plot locations, which are "fuzzed" to protect landowner privacy.

### Climate data

CHELSA version 2.1, at about 1 km resolution, aggregated to the grain size chosen in the app. Historical climate is the 1981-2010 average. Future climate comes from five CMIP6 global climate models (GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL) under three emissions pathways (SSP1-2.6, SSP3-7.0, SSP5-8.5) for 2011-2040, 2041-2070, and 2071-2100. Water balance variables (actual evapotranspiration and climatic water deficit) are calculated from monthly temperature and precipitation using the [`hydro`](https://github.com/matthewkling/hydro) R package.

### Human land use

The connectivity setting "What hinders movement" can weigh human land use alongside habitat quality, or in place of it. Land use comes from the Global Human Modification index for 2022 (version 3; Theobald et al. 2024, [doi:10.5281/zenodo.14502573](https://doi.org/10.5281/zenodo.14502573), CC BY 4.0), which rates each place from 0 (unmodified) to 1 (fully modified) by combining agriculture, development, transportation, energy, and other pressures. It is averaged from 300 m to about 1 km, and then to the grain size chosen in the app. Fully modified land is treated as being as hard to cross as land that is not habitat, and open water (which the index leaves blank) as not habitat. The index describes present-day land use, so it is applied unchanged to future maps.

### Boundaries

The polygon library is built from published sources and stored with the app, so it does not change between visits:

- National Park Service units and wilderness areas of at least 100 km²
- National forests and ranger districts (US Forest Service)
- Bureau of Land Management field offices
- States and counties (US Census Bureau, 2025)
- Level III ecoregions (US Environmental Protection Agency)

The `settings` record lists the version of each dataset behind a result.
