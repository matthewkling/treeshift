# Model settings (helpers/model_settings.R): the applied values, the Update
# map button that applies the controls' values, and the signs that the
# controls have moved ahead of the map.

# The applied values. They start as the link's settings, which the controls
# were built from (app.R), so a linked page computes its results straight
# away. reactiveValues only invalidates a field whose value actually changes,
# so applying recomputes only what depends on the settings that changed.
model <- do.call(reactiveValues, link$state[model_setting_ids])

model_menus <- model_setting_menus(url_ctx)

# The controls' values, by setting. One menu (scenario_choice) sets two: a
# storyline, whose model is resolved only when it applies, or a specific model.
controls_now <- reactive({
      v <- lapply(setNames(nm = model_setting_ids), function(id) input[[id]])
      choice <- input$scenario_choice
      if (is.null(choice)) {
            v[c("scenario_storyline", "scenario_model")] <-
                  isolate(list(model$scenario_storyline, model$scenario_model))
      } else if (choice %in% storyline_choices) {
            v[c("scenario_storyline", "scenario_model")] <- list(choice, isolate(model$scenario_model))
      } else {
            v[c("scenario_storyline", "scenario_model")] <- list(NULL, choice)
      }
      v
})

# Connectivity is solved only while it is on the map (server/connectivity.R),
# so its settings cost nothing while it isn't: they wait for Run model only
# while it is shown, and otherwise apply straight away (below).
connectivity_shown <- reactive(identical(input$stat, "connectivity"))

pending <- reactive({
      p <- pending_settings(controls_now(), reactiveValuesToList(model))
      if (connectivity_shown()) p else setdiff(p, connectivity_setting_ids)
})

# Apply the given values (by default, the controls'). Unchanged settings are
# left alone, so they don't invalidate anything. A storyline is resolved to
# its model here, from the applied bounding box and scenario, so nothing
# downstream sees a new bounding box with the old storyline's model.
apply_model_settings <- function(values = isolate(controls_now())) {
      applied <- isolate(reactiveValuesToList(model))
      for (id in names(values)) {
            v <- values[[id]]
            # A NULL storyline is a setting (a specific model); any other NULL
            # is a control not yet built
            if (is.null(v) && id != "scenario_storyline") next
            if (same_setting(v, applied[[id]])) next
            if (id %in% names(model_menus)) v <- url_menu_order(v, model_menus[[id]])
            model[[id]] <- v
            history_step <<- TRUE   # the link this produces gets its own Back step
      }
      s <- isolate(reactiveValuesToList(model))
      if (!is.null(s$scenario_storyline)) {
            # NULL if the bounding box has no climate data; scenario_file()
            # (server/climate.R) then says so
            m <- tryCatch(ts_storyline_model(core_data, s$scenario_storyline, s$scenario_ssp,
                                             s$scenario_period, s$bbox_x, s$bbox_y),
                          error = function(e) NULL)
            if (!same_setting(m, s$scenario_model)) model$scenario_model <- m
      }
}

observeEvent(input$update_model, apply_model_settings())

observe({
      req(!connectivity_shown())
      apply_model_settings(lapply(setNames(nm = connectivity_setting_ids), function(id) input[[id]]))
})

observe({
      n <- length(pending())
      updateActionButton(session, "update_model",
                         label = if (n) "Run model" else "Model is up to date",
                         disabled = n == 0)
})

output$pending_note <- renderUI({
      if (length(pending())) div(class = "ts-pending-note", pending_phrase(pending()))
})

# A light wash over the map while settings are pending: the map stays
# readable and usable (clicks pass through), but plainly not current
output$map_pending <- renderUI({
      if (length(pending())) div(class = "ts-map-stale")
})
