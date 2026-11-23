## MODULE AVAILABILITY -----------------------------------------
##
## A module that can't show anything until something outside it is set (a
## community run, a focal site) shows one notice in place of its controls and
## figures (module_body(), ui/right_sidebar.R), rather than an error in each
## plot. Hidden outputs aren't computed. Conditions a module's own controls
## can fix, like a scope with no polygon yet, are still reported in place,
## since hiding the controls would hide the fix.
##
## Each entry gives the reason the module is unavailable, or NULL when it's
## available. Its renderers req() the same entry, so nothing computes in the
## flush before the browser hides them.

no_site_msg <- "No focal site yet. Click the map to set one."

module_blockers <- list(
      comm = reactive({
            if (multi_mode()) return(NULL)
            if (identical(model$track, "All trees"))
                  "This module needs a community run. Under 'Entity' in the left sidebar, choose species, forest types or forest groups, then 'community composition'."
            else "This module needs a community run. Under 'Entity' in the left sidebar, choose 'community composition'."
      }),
      sensitivity = reactive({
            if (multi_mode()) "This module is available for single species, forest types and forest groups, not community runs."
            else if (!isTruthy(site$x) || !isTruthy(site$y)) no_site_msg
      }),
      # Held-out accuracy can only be scored for richness among the community
      # measures
      eval = reactive({
            if (multi_mode() && !identical(input$comm_measure, "richness"))
                  "Accuracy can only be measured for richness. Switch the community measure to richness, or switch back to a single species."
      }),
      analogs = reactive(if (!isTruthy(site$x) || !isTruthy(site$y)) no_site_msg)
)

for (.id in names(module_blockers)) local({
      id <- .id
      blocker <- module_blockers[[id]]
      output[[paste0(id, "_notice")]] <- renderUI({
            msg <- blocker()
            if (!is.null(msg)) div(class = "ts-module-notice", icon("circle-info"), span(msg))
      })
      # Read by the body's conditionalPanel, so it must update while the
      # body is hidden
      output[[paste0(id, "_ready")]] <- reactive(is.null(blocker()))
      outputOptions(output, paste0(id, "_ready"), suspendWhenHidden = FALSE)
})
