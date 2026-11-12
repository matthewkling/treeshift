# Climate ----------------------------------------------------------------------------

#' Future climate file paths for a scenario, one per climate model
ts_scenario_paths <- function(data, ssp, period, models) {
      file.path(data$clim_future_dir, paste0(ssp, "_", period, "_", models, ".tif"))
}

#' Historical and future climate, uncropped and at native resolution
#'
#' Restricted to `clim_vars`. With several climate models the future raster is
#' their layer-wise mean (an ensemble). Plot climates are extracted from these
#' rasters, so they reflect each plot's own location rather than a grid cell.
#'
#' @return list(hst, fut) of SpatRasters.
ts_climate <- function(data, clim_vars, ssp, period, models) {
      paths <- ts_scenario_paths(data, ssp, period, models)
      missing <- !file.exists(paths)
      if (any(missing))
            stop("No climate data for ", paste(basename(paths[missing]), collapse = ", "))

      r_hst <- terra::rast(data$clim_historical)[[clim_vars]]

      # Subset before summing: raster addition is positional, so this keeps
      # layers aligned by name across models. Reduce() makes the single-model
      # case a no-op.
      fut_list <- lapply(paths, function(p) terra::rast(p)[[clim_vars]])
      r_fut <- Reduce(`+`, fut_list) / length(fut_list)
      names(r_fut) <- clim_vars   # `+` can drop layer names

      list(hst = r_hst, fut = r_fut)
}

#' The modeling grid: climate cropped to the bounding box and aggregated
#'
#' @param climate `ts_climate()`'s list.
#' @param clim_agg Aggregation factor (1 keeps the native grid).
ts_climate_grid <- function(climate, bbox_x, bbox_y, clim_agg) {
      bb <- terra::ext(c(bbox_x, bbox_y))
      grid_of <- function(r) {
            r <- terra::crop(r, bb)
            if (clim_agg == 1) r else terra::aggregate(r, fact = clim_agg, na.rm = TRUE)
      }
      list(hst = grid_of(climate$hst), fut = grid_of(climate$fut))
}
