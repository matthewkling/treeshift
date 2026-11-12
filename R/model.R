# A whole run ------------------------------------------------------------------------
#
# The same stages the app runs, in order, for scripts and tests. The app calls
# the stages separately (server/), so each reruns only when its own inputs
# change.

#' Run the model
#'
#' @param settings `ts_settings()`.
#' @param data `ts_data()`.
#' @return A `ts_result` list: `settings`, `multi` (whether it is a community
#'   run), `climate` (native-resolution rasters), `grid` (the modeling grid),
#'   `pool_bbox`, `plots`, `pool`, `responses` (the per-entity matrix, for a
#'   community run) and `surfaces` (list(bsl, fut)).
ts_model <- function(settings, data) {
      s <- settings
      multi <- ts_multi(s)

      climate <- ts_climate(data, s$clim_vars, s$scenario_ssp, s$scenario_period,
                            s$scenario_model)
      grid <- ts_climate_grid(climate, s$bbox_x, s$bbox_y, s$clim_agg)
      pb <- ts_pool_bbox(s$bbox_x, s$bbox_y, s$theta_geog)
      plots <- ts_plots(data, pb, s$species)

      if (multi) {
            y <- ts_multi_response(data, plots, pb, s$track, s$fia_var)
            pool <- ts_pool(plots, climate)
            surfaces <- ts_multi_surfaces(grid, pool, y, s$theta_clim, s$theta_geog)
      } else {
            y <- NULL
            ratio <- ts_ratio_pieces(plots, s$track, s$fia_var)
            response <- if (is.null(ratio))
                  ts_response(plots, s$track, s$fia_var, s$forest_type, s$forest_group)
            if (is.null(ratio) && is.null(response))
                  stop("The ", s$track, " track has no variable '", s$fia_var, "'")
            pool <- ts_pool(plots, climate, ts_pool_values(plots, response, ratio))
            surfaces <- ts_surfaces(grid, pool, s$fia_var, s$theta_clim, s$theta_geog)
      }

      structure(list(settings = s, multi = multi, climate = climate, grid = grid,
                     pool_bbox = pb, plots = plots, pool = pool, responses = y,
                     surfaces = surfaces),
                class = "ts_result")
}

#' The layer a run maps: the variable, or for a community run a measure
#'
#' @param result `ts_model()`'s result.
#' @param measure For a community run, the community measure.
#' @param support_floor See `ts_value_layer()` and `ts_community_layer()`.
#' @param top_n For a categorical measure, see `ts_community_categories()`.
ts_layer <- function(result, time = "historical", measure = NULL, support_floor = 0,
                     top_n = 10) {
      if (!result$multi)
            return(ts_value_layer(result$surfaces, result$settings$fia_var, time, support_floor))
      if (is.null(measure)) stop("A community run needs a `measure`")
      if (measure %in% categorical_measures)
            return(ts_community_categories(result$surfaces, measure, time,
                                           support_floor, top_n)$raster)
      ts_community_layer(result$surfaces, measure, time, support_floor)
}
