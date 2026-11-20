# Model settings ------------------------------------------------------------------
#
# The settings that determine a model run: what is modeled (track, entity scope,
# entity, variable), where and at what grain (bounding box, aggregation factor),
# against which climates (variables, scenario) and with which bandwidths. These
# are the app's model settings (helpers/model_settings.R) minus connectivity.
# Display choices -- timeframe, community measure, support floor -- are
# arguments to the layer functions instead (R/layers.R), since one run serves
# all of them.

#' Default model settings
#'
#' The app's starting settings. `forest_type` and `forest_group` have no
#' default; they are needed only on those tracks.
ts_default_settings <- function() {
      list(track = "Species", entity_scope = "single", species = "Pinus edulis",
           forest_type = NULL, forest_group = NULL,
           stage = "any", metric = "occurrence", fia_var = NULL,
           clim_agg = 10, bbox_x = c(-116, -100), bbox_y = c(29, 43),
           clim_vars = c("aet", "cwd"),
           scenario_ssp = "ssp370", scenario_period = "2041-2070",
           scenario_storyline = "median", scenario_model = NULL,
           theta_geog = 30, theta_clim = 0.1)
}

#' Model settings, from the defaults
#'
#' @param ... Settings to change from `ts_default_settings()`, by name.
#'   `fia_var` is resolved from `track`, `metric` and `stage` unless given.
#'   The climate model is either a storyline (`scenario_storyline`, resolved to
#'   a model by `ts_resolve_scenario()`, which needs the data) or a specific
#'   `scenario_model`; giving a model without a storyline means the latter.
#' @return A `ts_settings` list.
ts_settings <- function(...) {
      changes <- list(...)
      d <- ts_default_settings()
      unknown <- setdiff(names(changes), names(d))
      if (length(unknown)) stop("Unknown settings: ", paste(unknown, collapse = ", "))
      s <- d
      if ("scenario_model" %in% names(changes) && !"scenario_storyline" %in% names(changes))
            s["scenario_storyline"] <- list(NULL)
      s[names(changes)] <- changes   # unlike modifyList(), keeps NULLs as given

      if (is.null(s$fia_var)) s$fia_var <- resolve_var(s$track, s$metric, s$stage)
      if (is.null(s$fia_var))
            stop("No variable for metric '", s$metric, "' at stage '", s$stage, "'")
      if (s$track == "Forest type" && !ts_multi(s) && is.null(s$forest_type))
            stop("The Forest type track needs `forest_type`")
      if (s$track == "Forest group" && !ts_multi(s) && is.null(s$forest_group))
            stop("The Forest group track needs `forest_group`")
      if (!is.null(s$scenario_storyline)) {
            if (!isTRUE(s$scenario_storyline %in% storyline_choices))
                  stop("Unknown climate storyline '", s$scenario_storyline, "'")
            if (!is.null(s$scenario_model))
                  stop("Give a climate storyline or a climate model, not both")
      }
      if (!all(s$clim_vars %in% clim_vars()))
            stop("Unknown climate variables: ",
                 paste(setdiff(s$clim_vars, clim_vars()), collapse = ", "))
      structure(s, class = "ts_settings")
}

#' Climate storylines: the models ranked least, median and most by change in
#' climatic water deficit (drying) over the bounding box (`ts_storyline_model()`)
storyline_choices <- c("Least drying" = "least", "Median drying" = "median", "Most drying" = "most")

#' Whether settings call for a community (all-entity) run
#'
#' All trees is already an aggregate, so it never runs per entity.
ts_multi <- function(settings) {
      identical(settings$entity_scope, "multi") && !identical(settings$track, "All trees")
}
