# A stand-in for data/clim_deltas.tif (data_prep/clim_prep.R): a 4 x 4 grid of
# 1-degree cells over x -110..-106, y 38..42, one <ssp>_<period>_<model>_cwd
# layer per model. `cwd` gives each model's change, a number (everywhere) or
# 16 cell values, row by row from the northwest.
fake_cube <- function(cwd, dir, ssp = "ssp370", period = "2041-2070") {
      layers <- lapply(names(cwd), function(m)
            terra::rast(nrows = 4, ncols = 4, xmin = -110, xmax = -106, ymin = 38, ymax = 42,
                        crs = "EPSG:4326", vals = rep_len(cwd[[m]], 16)))
      cube <- terra::rast(layers)
      names(cube) <- paste(ssp, period, names(cwd), "cwd", sep = "_")
      path <- file.path(dir, "clim_deltas.tif")
      terra::writeRaster(cube, path, overwrite = TRUE)
      list(scenarios = tibble::tibble(ssp = ssp, period = period, model = names(cwd)),
           clim_deltas = path)
}
