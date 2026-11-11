# REACTIVE UI COMPONENTS ---------------------------

## Variable selection is a (track x metric x stage) cascade that resolves to
## a single canonical `fia_var` string, which the rest of the app switches
## on. These helpers define the valid choices and the resolution.
## Repopulate variable choices when the track or life stage changes.
## Selection is preserved when it survives the new choice set, so moving
## between two stages that both offer "occurrence" does not bounce the
## user back to the first option.
observeEvent(list(input$track, input$stage), {
      mc <- metric_choices(input$track, input$stage)
      cur <- if (isTruthy(input$metric)) input$metric else ""
      updateSelectInput(session, "metric", choices = mc,
                        selected = if (cur %in% mc) cur else mc[1])
}, ignoreInit = TRUE)

## Entity-scope option text follows the track ("single species" vs "single
## forest type"); only the labels change, the values stay single/multi so
## nothing downstream is affected.
##
## observeEvent rather than observe: the handler body is isolated, so
## reading input$entity_scope here does not create a dependency on the
## input this writes to. A bare observe() would re-fire on its own write.
observeEvent(input$track, {
      req(input$track != "All trees")
      updateSelectInput(session, "entity_scope", choices = entity_scope_choices(input$track),
                        selected = if (isTruthy(input$entity_scope)) input$entity_scope else "single")
}, ignoreInit = TRUE)

## Resolve (track, metric, stage) -> hidden canonical fia_var.
observeEvent(list(input$track, input$metric, input$stage), {
      req(input$metric)
      v <- resolve_var(input$track, input$metric, input$stage)
      req(!is.null(v))
      updateSelectInput(session, "fia_var", selected = v)
})

## Seed the climate-model chart's x/y axis variables from the user's
## left-sidebar climate-variable selection. Done once, on the first
## populated value of input$clim_vars, so it provides a sensible default
## (the first two analog variables) without clobbering a deliberate axis
## choice the user makes later. The dropdowns still offer all six
## variables -- only the initial default tracks the sidebar.
.model_delta_axes_seeded <- reactiveVal(FALSE)
observeEvent(input$clim_vars, {
      if (isTRUE(.model_delta_axes_seeded())) return()
      cv <- input$clim_vars
      req(length(cv) >= 1)
      vx <- cv[1]
      vy <- if (length(cv) >= 2) cv[2] else cv[1]
      updateSelectInput(session, "model_delta_x", selected = vx)
      updateSelectInput(session, "model_delta_y", selected = vy)
      .model_delta_axes_seeded(TRUE)
})

## Keep the climate-relationships module's variable menus in sync with the
## climate variables selected in the left sidebar. Its plot indexes into
## clim()$hst, which only holds the selected variables, so the options must
## be exactly that set. Preserve the current pick when it is still
## available, otherwise fall back to the first selected variable.
## cr_v1 offers every display variable, not just the modeled ones, so its
## choices no longer depend on input$clim_vars and the observer that used
## to narrow them is gone. The UI declares the full set directly.

## Spatial-scope status lines. These exist so the left panel states its
## own preconditions: a right-sidebar panel that renders empty because no
## site or polygon is set now has a visible cause rather than looking
## broken. Both are uiOutput rather than textOutput so the "not set" case
## can carry its own instruction.
output$site_status <- renderUI({
      if (is.null(site$x) || is.null(site$y)) {
            helpText("None selected \u2014 click the map to set one.")
      } else {
            helpText(sprintf("lat %.3f, lon %.3f", site$y, site$x))
      }
})


## Dynamic snap-to-range button (label reflects focal entity)
output$snap_focal_btn <- renderUI({
      switch(input$track,
             "Species"      = actionButton("focal_domain", "Snap to species range"),
             "Forest type"  = actionButton("focal_domain", "Snap to forest type range"),
             "Forest group" = actionButton("focal_domain", "Snap to forest group range"),
             "All trees"    = tags$button(class = "btn btn-default",
                                          disabled = NA,
                                          style = "opacity: 0.5; cursor: not-allowed;",
                                          title = "No focal entity in this mode",
                                          "Snap to focal range"))
})

## Snap domain to focal species/forest-type/forest-group range
observeEvent(input$focal_domain, {
      pad <- 1
      bb <- switch(input$track,
                   "Species" = fia_plots_species_all() %>%
                         filter(species == input$species),
                   "Forest type" = fia_plots_species_all() %>%
                         filter(forest_type == as.integer(input$forest_type)),
                   "Forest group" = {
                         ft_in_group <- forest_type_codes %>%
                               filter(forest_group == as.integer(input$forest_group)) %>%
                               pull(forest_type)
                         fia_plots_species_all() %>%
                               filter(forest_type %in% ft_in_group)
                   },
                   NULL)
      req(!is.null(bb), nrow(bb) > 0)
      bb <- bb %>%
            summarize(xmin = floor(min(x)) - pad,
                      xmax = ceiling(max(x)) + pad,
                      ymin = floor(min(y)) - pad,
                      ymax = ceiling(max(y)) + pad)
      # bbox_x / bbox_y are sliderInputs, so updateSliderInput(value=) is
      # the matching call. These previously used updateSelectInput(selected=)
      # and worked, but that pairing is not the documented one and is not
      # something to rely on.
      updateSliderInput(session, "bbox_x", value = c(bb$xmin, bb$xmax))
      updateSliderInput(session, "bbox_y", value = c(bb$ymin, bb$ymax))
})

## Snap domain to the focal polygon. Covers the buffer too, since the buffered
## summaries need it inside the box; whole degrees, like the species snap.
observeEvent(input$polygon_snap, {
      req(user_polygon_buffered())
      bb <- sf::st_bbox(user_polygon_buffered())
      updateSliderInput(session, "bbox_x", value = c(floor(bb[["xmin"]]), ceiling(bb[["xmax"]])))
      updateSliderInput(session, "bbox_y", value = c(floor(bb[["ymin"]]), ceiling(bb[["ymax"]])))
})

## Snap domain to current map view
observeEvent(input$window_domain, {
      bb <- input$map_bounds
      updateSliderInput(session, "bbox_x", value = c(bb$west, bb$east))
      updateSliderInput(session, "bbox_y", value = c(bb$south, bb$north))
})

# for circuitscape rasters, time period must be "delta"
observe({
      if(input$stat == "connectivity"){
            updateSelectInput(session, "time", selected = "delta")
      }
})

## Show or hide the connectivity entries in the left-sidebar "Map layers"
## menus, driven by the "Enable connectivity" switch in the right-sidebar
## Connectivity panel. Connectivity is a niche analysis, so its options are
## kept out of the core menus until the user asks for them.
##
## updateSelectInput() resets the selection to the first choice whenever
## `choices` changes, so `selected` is always passed explicitly. On enable,
## the current pick is preserved; on disable, a connectivity pick is no
## longer offered and falls back to the menu default.
##
## ignoreInit: the UI is built with the connectivity-off choices and the
## switch defaults to FALSE, so the two already agree at startup.
## Explains the withdrawal rather than leaving a switch that has silently
## stopped doing anything.
output$cs_note <- renderUI({
      req(multi_mode())
      div(class = "ts-notice",
          "Connectivity models the movement of one species or forest type between its
                 current and future range, so it does not apply to community composition. For
                 connectivity of the forest as a whole, set the entity to 'all trees'.")
})

## Also fires on multi_mode(): connectivity is a single-entity analysis.
## Its electrodes are built from one entity's historical and future
## surfaces and then solved, so under community composition it would
## silently describe whichever species the hidden picker still holds --
## producing something that looks like community connectivity and is not.
##
## There is no community version worth building here. The one that would
## be well defined is connectivity of an aggregate surface, and the All
## trees track already provides exactly that.
observeEvent(list(input$cs_enable, multi_mode()), {
      on <- isTRUE(input$cs_enable) && !multi_mode()

      stat_choices <- if(on) stat_choices_cs else stat_choices_base
      vector_choices <- if(on) vector_choices_cs else vector_choices_base

      keep <- function(cur, choices, default){
            if(isTruthy(cur) && cur %in% choices) cur else default
      }

      updateSelectInput(session, "stat", choices = stat_choices,
                        selected = keep(input$stat, stat_choices, "FIA variable"))
      updateSelectInput(session, "vector", choices = vector_choices,
                        selected = keep(input$vector, vector_choices, "occupied FIA plots"))
}, ignoreInit = TRUE)

# Swap eval_stat choices to the metric set cv_performance() produces for
# the current response type. Binary (presence-type) responses get
# ranking/threshold metrics; continuous responses get error metrics.
#
# Keyed on response_is_binary() rather than on fia_var directly: in
# community mode fia_var is still "presence probability" -- that is the
# per-entity variable being summed -- but the response being scored is
# richness, a count. Testing fia_var offered AUC and TSS for it.
#
# observeEvent so the handler is isolated. It also runs at startup: the UI is
# declared with the binary set, and a linked page can open on a continuous
# variable.
observeEvent(response_is_binary(), {
      if (response_is_binary()) {
            updateSelectInput(session, "eval_stat",
                              choices = c("AUC", "TSS", "Brier"),
                              selected = "AUC")
      } else {
            updateSelectInput(session, "eval_stat",
                              choices = c("R-squared", "RMSE"),
                              selected = "R-squared")
      }
})

# highlight focal site analogs: activate correct vector layer; zoom to site
observeEvent(input$highlight_site, {
      updateSelectInput(session, "vector", selected = "focal site analogs")

      if(!is.null(site$y)){
            circle <- geo_circle(cbind(site$x, site$y), input$theta_geog * 3) %>%
                  as.data.frame()

            leafletProxy("map") %>%
                  fitBounds(lng1 = min(circle$lon),
                            lat1 = min(circle$lat),
                            lng2 = max(circle$lon),
                            lat2 = max(circle$lat))
      }
})

# add user-provided polygon (and its buffered version)
observe({
      proxy <- leafletProxy("map") %>%
            clearGroup("user_polygon")

      # Re-draw after any raster change by taking dependency on display_layer
      display_layer()

      # Visibility toggle (Cartography panel). When off, the group is
      # already cleared above; req() stops before redrawing. Toggling
      # back on re-fires this observer and redraws. Mirrors show_site.
      req(input$show_polygon)

      req(user_polygon())
      poly <- user_polygon()

      proxy %>%
            addPolygons(
                  data = user_polygon_buffered(),
                  color = lighten_color(polygon_color()),
                  weight = 2,
                  dashArray = "5, 5",
                  fill = FALSE,
                  group = "user_polygon"
            ) %>%
            addPolygons(
                  data = poly,
                  color = polygon_color(),
                  weight = 2,
                  fillOpacity = 0.1,
                  fillColor = polygon_color(),
                  group = "user_polygon",
                  label = ~name
            )
})
