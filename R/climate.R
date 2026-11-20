# Climate ----------------------------------------------------------------------------

#' Future climate file path for a scenario and climate model
ts_scenario_path <- function(data, ssp, period, model) {
      file.path(data$clim_future_dir, paste0(ssp, "_", period, "_", model, ".tif"))
}

#' Mean climate change per model over a bounding box, ranked
#'
#' From the coarse cube of climate changes (`data$clim_deltas`, built in
#' data_prep/clim_prep.R): the unweighted mean change in `var` over the box, for
#' each model available for `ssp` x `period`. A box smaller than one cube cell
#' reads the cell it falls in.
#'
#' @param var A climate variable (see `clim_vars()`), in native units.
#' @return A numeric vector named by model, smallest change first. Ties keep
#'   model order.
ts_storyline_ranking <- function(data, ssp, period, bbox_x, bbox_y, var = "cwd") {
      models <- sort(data$scenarios$model[data$scenarios$ssp == ssp &
                                                data$scenarios$period == period])
      if (!length(models)) stop("No climate models for ", ssp, " ", period)
      layers <- paste(ssp, period, models, var, sep = "_")
      cube <- terra::rast(data$clim_deltas)
      missing <- setdiff(layers, names(cube))
      if (length(missing))
            stop("The climate change summary (", basename(data$clim_deltas), ") lacks ",
                 paste(missing, collapse = ", "), "; rebuild it with data_prep/clim_prep.R")

      r <- terra::crop(cube[[layers]], terra::ext(c(bbox_x, bbox_y)), snap = "out")
      means <- stats::setNames(terra::global(r, "mean", na.rm = TRUE)$mean, models)
      if (anyNA(means)) stop("No climate data in the bounding box")
      means[order(means)]
}

#' The climate model a storyline stands for, over a bounding box
#'
#' Ranks the models available for `ssp` x `period` by mean change in `var` over
#' the box (`ts_storyline_ranking()`): "least" is the smallest change, "most"
#' the largest, "median" rank ceiling(n/2).
#'
#' @param storyline One of `storyline_choices`.
ts_storyline_model <- function(data, storyline, ssp, period, bbox_x, bbox_y, var = "cwd") {
      if (!storyline %in% storyline_choices) stop("Unknown climate storyline '", storyline, "'")
      ranking <- ts_storyline_ranking(data, ssp, period, bbox_x, bbox_y, var)
      n <- length(ranking)
      names(ranking)[switch(storyline, least = 1, median = ceiling(n / 2), most = n)]
}

#' Settings with the storyline (if any) resolved to its climate model
#'
#' `ts_settings()` holds a storyline unresolved, since resolving it needs the
#' data. Settings without a storyline, or already resolved (the app resolves
#' when it applies settings, and a link keeps the model it was made with), are
#' returned unchanged.
ts_resolve_scenario <- function(settings, data) {
      if (is.null(settings$scenario_storyline) || !is.null(settings$scenario_model))
            return(settings)
      settings$scenario_model <- ts_storyline_model(
            data, settings$scenario_storyline, settings$scenario_ssp,
            settings$scenario_period, settings$bbox_x, settings$bbox_y)
      settings
}

#' The settings' climate model, which must have been resolved
ts_scenario_model <- function(settings) {
      if (is.null(settings$scenario_model))
            stop("No climate model: resolve the storyline with ts_resolve_scenario() first")
      settings$scenario_model
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
