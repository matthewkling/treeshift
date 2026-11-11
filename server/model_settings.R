# Model settings (helpers/model_settings.R): the applied values, the Update
# map button that applies the controls' values, and the signs that the
# controls have moved ahead of the map.

# The applied values. They start as the link's settings, which the controls
# were built from (app.R), so a linked page computes its results straight
# away. reactiveValues only invalidates a field whose value actually changes,
# so applying recomputes only what depends on the settings that changed.
model <- do.call(reactiveValues, link$state[model_setting_ids])

model_menus <- model_setting_menus(url_ctx)

controls_now <- reactive(lapply(setNames(nm = model_setting_ids), function(id) input[[id]]))

pending <- reactive(pending_settings(c(controls_now(), list(cs_enable = input$cs_enable)),
                                     reactiveValuesToList(model)))

# Apply the given values (by default, the controls'). Unchanged settings are
# left alone, so they don't invalidate anything.
apply_model_settings <- function(values = isolate(controls_now())) {
      applied <- isolate(reactiveValuesToList(model))
      for (id in names(values)) {
            v <- values[[id]]
            if (is.null(v) || same_setting(v, applied[[id]])) next
            if (id %in% names(model_menus)) v <- url_menu_order(v, model_menus[[id]])
            model[[id]] <- v
      }
}

observeEvent(input$update_model, apply_model_settings())

observe({
      n <- length(pending())
      updateActionButton(session, "update_model",
                         label = if (n) "Update map" else "Map is up to date",
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
