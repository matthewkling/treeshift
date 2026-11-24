## CLIMATE DATA -----------------------------------------

# Raw (uncropped, unaggregated) historical and future climate rasters,
# restricted to the selected climate variables. The future raster depends
# on the scenario selection; the historical raster is shared across
# scenarios. Used as the source for both the modeling grid (via clim())
# and for plot-level climate extraction (in pool(), site_analogs(), and
# eval_var_set()).

# Active scenario filename, reconstructed from the three scenario inputs.
# A missing file surfaces a user-visible message via validate(), so
# downstream outputs (the map, the focal-site plots, etc.) display it
# instead of erroring.
scenario_file <- reactive({
      req(isTruthy(model$scenario_ssp), isTruthy(model$scenario_period))
      shiny::validate(need(isTruthy(model$scenario_model),
            "No climate model could be chosen: the bounding box has no climate data. Move it over the contiguous US."))
      path <- ts_scenario_path(core_data, model$scenario_ssp, model$scenario_period,
                               model$scenario_model)
      shiny::validate(need(file.exists(path),
            "Selected scenario combination is not available. Try a different emissions pathway, era, or climate model."))
      path
})

clim_pool <- reactive({
      # scenario_file() first, so a missing scenario shows its message rather
      # than ts_climate()'s error
      scenario_file()
      ts_climate(core_data, model$clim_vars, model$scenario_ssp, model$scenario_period,
                 model$scenario_model)
})

clim <- reactive(ts_climate_grid(clim_pool(), model$bbox_x, model$bbox_y, model$clim_agg))


## DISPLAY-ONLY COVARIATES ------------------------------------
##
## The x-axis variable for the climate-relationships plot, resolved
## independently of clim_pool(). Nothing here feeds the analog model: the
## point is that a user can plot against a variable the model was not
## fitted on, and the isolation is what makes that safe.
##
## The cost is one extra crop/aggregate per variable change. Cheap relative
## to a prediction, and it only recomputes when the axis variable, the
## bbox, the grain, or the scenario changes.
cr_grid <- reactive({
      v <- input$cr_v1
      req(isTruthy(v), v %in% display_vars())
      conus <- rast(core_data$conus)
      grid_of <- function(r) ts_grid_of(r, conus, model$bbox_x, model$bbox_y, model$clim_agg)

      if (!is.null(display_extra_vars[[v]])) {
            # A static covariate: the same raster for both timeframes,
            # since elevation does not change between them.
            g <- grid_of(rast(display_extra_vars[[v]]$file))
            names(g) <- v
            return(list(hst = g, fut = g))
      }

      list(hst = grid_of(rast("data/clim_historical.tif")[[v]]),
           fut = grid_of(rast(scenario_file())[[v]]))
})

## Same variable at the FIA plot coordinates, row-aligned with plot_occ().
cr_plot_values <- reactive({
      v <- input$cr_v1
      req(isTruthy(v), v %in% display_vars())
      xy <- as.data.frame(plot_occ()[, c("x", "y")])

      src <- if (!is.null(display_extra_vars[[v]])) {
            rast(display_extra_vars[[v]]$file)
      } else if (isTRUE(input$cr_clim == "future")) {
            rast(scenario_file())[[v]]
      } else {
            rast("data/clim_historical.tif")[[v]]
      }
      as.numeric(terra::extract(src, xy, ID = FALSE, method = "bilinear")[, 1])
})

## Whether a focal site is set, for the site color control's
## conditionalPanel (ui/left_sidebar.R), so it must update while hidden
output$site_set <- reactive(!is.null(site$x))
outputOptions(output, "site_set", suspendWhenHidden = FALSE)

## Resolved once so both the outline and its buffer ring agree, and so a
## typed value that R cannot parse falls back instead of erroring inside
## leaflet.
valid_color <- function(v, fallback) {
      if (!isTruthy(v)) return(fallback)
      ok <- tryCatch({ grDevices::col2rgb(v); TRUE }, error = function(e) FALSE)
      if (ok) v else fallback
}
polygon_color <- reactive(valid_color(input$polygon_color, "darkred"))
site_color    <- reactive(valid_color(input$site_color, "#4c32a8"))

## Restrict a raster to the selected area. Applied to the COMBINED raster
## rather than to each layer, because cropping them separately would leave
## the response and the covariate on different extents.
cr_scope_mask <- function(r) {
      sc <- if (isTruthy(input$cr_scope)) input$cr_scope else "domain"
      if (identical(sc, "domain")) return(r)

      if (sc %in% c("polygon", "buffered")) {
            req(isTruthy(polygon_active()))
            pv <- terra::vect(if (identical(sc, "polygon")) user_polygon()
                              else user_polygon_buffered())
            return(terra::mask(terra::crop(r, pv), pv))
      }

      # Focal site: the neighborhood the model draws analogs from, not the
      # single cell -- one cell would give one point and no relationship.
      req(!is.null(site$x), !is.null(site$y))
      b <- terra::buffer(terra::vect(cbind(site$x, site$y), crs = "EPSG:4326"),
                         width = model$theta_geog * 3 * 1000)
      terra::mask(terra::crop(r, b), b)
}

## Is the x-axis variable one the model was actually fitted on? A
## relationship against anything else is a property of the landscape rather
## than of the model, and the plot says so rather than leaving it to be
## inferred.
cr_var_in_model <- reactive({
      isTruthy(input$cr_v1) && input$cr_v1 %in% model$clim_vars
})
