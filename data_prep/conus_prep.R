# The contiguous US on the climate grid: data/conus.tif
#
#   Rscript data_prep/conus_prep.R
#
# Run from the project root, after data_prep/boundaries_prep.R. The app makes
# no estimates outside the contiguous US, where there are no FIA plots: the
# climate grid is masked to this raster before every run (ts_climate_grid()).
#
# Source: the boundary library's state units (the 48 contiguous states and the
# District of Columbia, from the Census Bureau's 2025 cartographic boundaries),
# dissolved. A cell is in when its center is.
#
# Output: data/conus.tif, 1 inside and NA outside, on the grid of
# data/clim_historical.tif.

library(terra)

template <- rast("data/clim_historical.tif")[[1]]
states <- vect("data/boundaries.gpkg", layer = "units",
               query = "SELECT * FROM units WHERE id LIKE 'state:%'")
stopifnot(nrow(states) == 49)

conus <- rasterize(aggregate(states), template, field = 1, background = NA) |>
      setNames("conus")

writeRaster(conus, "data/conus.tif", datatype = "INT1U",
            gdal = c("COMPRESS=DEFLATE"), overwrite = TRUE)
