## Data sources

### Forest data

US Forest Service Forest Inventory and Analysis database (FIADB), downloaded June 2025. TreeShift uses the most recent visit to each plot, on the standard four-subplot design. Plots with evidence of planting are excluded, so the data describe naturally established forests. "Non-forest" FIA plot locations that were not visited in the field are used only where LANDFIRE Biophysical Settings indicate non-forest vegetation, so they can stand in for places without trees. This tool uses public FIA plot locations, which are "fuzzed" to protect landowner privacy.

### Climate data

CHELSA version 2.1, at about 1 km resolution, aggregated to the grain size chosen in the app. Recent climate is the 2000-2019 average. Future climate comes from five CMIP6 global climate models (GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL) under three emissions pathways (SSP1-2.6, SSP3-7.0, SSP5-8.5) for 2011-2040, 2041-2070, and 2071-2100. Water balance variables (actual evapotranspiration and climatic water deficit) are calculated from monthly temperature and precipitation using the [`hydro`](https://github.com/matthewkling/hydro) R package.

### Boundaries

The polygon library is built from published sources and stored with the app, so it does not change between visits:

- National Park Service units and wilderness areas of at least 100 km²
- National forests and ranger districts (US Forest Service)
- Bureau of Land Management field offices
- States and counties (US Census Bureau, 2025)
- Level III ecoregions (US Environmental Protection Agency)

The `settings` record lists the version of each dataset behind a result.
