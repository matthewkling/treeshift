# Climate ----------------------------------------------------------------------------

#' Future climate file path for a scenario and climate model
ts_scenario_path <- function(data, ssp, period, model) {
      file.path(data$clim_future_dir, paste0(ssp, "_", period, "_", model, ".tif"))
}

#' Historical and future climate, uncropped and at native resolution
#'
#' Restricted to `clim_vars`, with the contiguous US mask (data_prep/conus_prep.R)
#' that `ts_climate_grid()` applies. Plot climates are extracted from these
#' rasters, so they reflect each plot's own location rather than a grid cell.
#'
#' @param model One climate model. Climates are never averaged across models:
#'   a model's response to an averaged climate isn't the average of its
#'   responses to each model's climate.
#' @return list(hst, fut, conus) of SpatRasters.
ts_climate <- function(data, clim_vars, ssp, period, model) {
      if (length(model) != 1) stop("`model` must be a single climate model")
      path <- ts_scenario_path(data, ssp, period, model)
      if (!file.exists(path)) stop("No climate data for ", basename(path))

      list(hst = terra::rast(data$clim_historical)[[clim_vars]],
           fut = terra::rast(path)[[clim_vars]],
           conus = terra::rast(data$conus))
}

#' The modeling grid: climate cropped to the bounding box, masked to the
#' contiguous US, and aggregated
#'
#' There are no FIA plots outside the contiguous US, so nothing is estimated
#' there. The mask is applied at native resolution, so a coarse cell on the
#' border takes the climate of its US part.
#'
#' @param climate `ts_climate()`'s list.
#' @param clim_agg Aggregation factor (1 keeps the native grid).
ts_climate_grid <- function(climate, bbox_x, bbox_y, clim_agg) {
      list(hst = ts_grid_of(climate$hst, climate$conus, bbox_x, bbox_y, clim_agg),
           fut = ts_grid_of(climate$fut, climate$conus, bbox_x, bbox_y, clim_agg))
}

#' One raster on the modeling grid (see `ts_climate_grid()`)
#'
#' @param conus The contiguous US mask, `terra::rast(data$conus)`.
ts_grid_of <- function(r, conus, bbox_x, bbox_y, clim_agg) {
      bb <- terra::ext(c(bbox_x, bbox_y))
      r <- terra::mask(terra::crop(r, bb), terra::crop(conus, bb))
      if (clim_agg == 1) r else terra::aggregate(r, fact = clim_agg, na.rm = TRUE)
}
