# Model settings: the inputs that rerun a model. They take effect when the
# user presses "Run model" (server/model_settings.R), so changing several in
# a row doesn't compute a model for each intermediate combination. Everything
# else -- timeframe, community measure, support floor, polygon, focal site,
# map layers, cartography -- applies right away.
#
# Server code reads these through `model$<id>`, the applied values, never
# `input$<id>`: a direct read would mix the pending settings into results
# computed from the applied ones (tests/testthat/test-model_settings.R).

model_setting_ids <- c("track", "entity_scope", "species", "forest_type", "forest_group",
                       "stage", "metric", "fia_var",
                       "clim_agg", "bbox_x", "bbox_y",
                       "clim_vars", "scenario_ssp", "scenario_period", "scenario_model",
                       "theta_geog", "theta_clim",
                       "cs_sources", "cs_targets", "cs_conductance_layer", "cs_loss_layer",
                       "cs_conductance_scale", "cs_loss_scale")

# How each is named in the "what changed" line under the button
model_setting_labels <- c(track = "entity", entity_scope = "entity", species = "species",
                          forest_type = "forest type", forest_group = "forest group",
                          stage = "life stage", metric = "variable", fia_var = "variable",
                          clim_agg = "grain size", bbox_x = "bounding box", bbox_y = "bounding box",
                          clim_vars = "climate variables", scenario_ssp = "emissions pathway",
                          scenario_period = "time period", scenario_model = "climate models",
                          theta_geog = "geographic bandwidth", theta_clim = "climate bandwidth",
                          cs_sources = "connectivity settings", cs_targets = "connectivity settings",
                          cs_conductance_layer = "connectivity settings",
                          cs_loss_layer = "connectivity settings",
                          cs_conductance_scale = "connectivity settings",
                          cs_loss_scale = "connectivity settings")

# Multi-selects, kept in menu order so that reordering a selection isn't a
# change
model_setting_menus <- function(ctx) list(clim_vars = clim_vars(), scenario_model = ctx$models)

# Whether two values of a setting are the same. Numbers compare by value (a
# slider's integer equals the link's double); multi-selects compare as sets.
same_setting <- function(a, b) {
      if (!length(a) || !length(b)) return(length(a) == length(b))
      if (is.numeric(a) && is.numeric(b))
            return(length(a) == length(b) && isTRUE(all.equal(as.numeric(a), as.numeric(b))))
      length(a) == length(b) && setequal(as.character(a), as.character(b))
}

# The settings that apply given the others, as the link schema decides
# (helpers/url_state.R). `s` needs the model settings plus cs_enable.
model_settings_in_use <- function(s) {
      single <- !url_multi(s)
      cs <- grep("^cs_", model_setting_ids, value = TRUE)
      model_setting_ids[c(
            track = TRUE, entity_scope = !identical(s$track, "All trees"),
            species = single && identical(s$track, "Species"),
            forest_type = single && identical(s$track, "Forest type"),
            forest_group = single && identical(s$track, "Forest group"),
            stage = url_has_stage(s), metric = TRUE, fia_var = TRUE,
            clim_agg = TRUE, bbox_x = TRUE, bbox_y = TRUE,
            clim_vars = TRUE, scenario_ssp = TRUE, scenario_period = TRUE, scenario_model = TRUE,
            theta_geog = TRUE, theta_clim = TRUE,
            stats::setNames(rep(url_cs_on(s), length(cs)), cs))[model_setting_ids]]
}

# The settings in use that differ between the controls and the applied
# values. A change to a setting that doesn't apply -- a forest type picked
# while mapping a species -- isn't pending.
pending_settings <- function(current, applied) {
      ids <- model_settings_in_use(current)
      ids[!vapply(ids, function(id) same_setting(current[[id]], applied[[id]]), logical(1))]
}

# "Species and geographic bandwidth changed."
pending_phrase <- function(ids) {
      ids <- ids[!duplicated(model_setting_labels[ids])]
      p <- stale_phrase(ids, model_setting_labels)
      paste0(toupper(substr(p, 1, 1)), substring(p, 2), " changed.")
}
