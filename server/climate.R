## CLIMATE DATA -----------------------------------------

# Raw (uncropped, unaggregated) historical and future climate rasters,
# restricted to the selected climate variables. The future raster depends
# on the scenario selection; the historical raster is shared across
# scenarios. Used as the source for both the modeling grid (via clim())
# and for plot-level climate extraction (in pool(), site_analogs(), and
# eval_var_set()).

# Active scenario filename(s), reconstructed from the three scenario
# inputs. scenario_model may hold multiple GCMs; this returns one path
# per selected model (the future raster is later averaged across them
# to form a model ensemble). Strict: if ANY selected combination is
# missing a file on disk, surfaces a user-visible message via validate()
# naming the missing models, rather than silently averaging a subset.
# Downstream outputs (the map, the focal-site plots, etc.) display the
# message instead of erroring.
scenario_file <- reactive({
      req(isTruthy(input$scenario_ssp),
          isTruthy(input$scenario_period),
          isTruthy(input$scenario_model))
      f <- paste0(input$scenario_ssp, "_",
                  input$scenario_period, "_",
                  input$scenario_model, ".tif")
      paths <- file.path("data/clim_future", f)
      exists <- file.exists(paths)

      shiny::validate(need(
            all(exists),
            if (length(paths) == 1) {
                  "Selected scenario combination is not available. Try a different emissions pathway, era, or climate model."
            } else {
                  paste0("No data for these model(s) at the selected pathway/era: ",
                         paste(toupper(input$scenario_model[!exists]), collapse = ", "),
                         ". Deselect them or choose a different pathway/era.")
            }))
      paths
})

clim_pool <- reactive({
      # Resolve scenario_file() at the top level so a req()/validate()
      # short-circuit (e.g., no selected (ssp, period, model) combination
      # exists on disk) propagates as a clean reactive cancellation
      # rather than an rast() argument-evaluation error.
      sf <- scenario_file()
      r_hst <- rast("data/clim_historical.tif")[[input$clim_vars]]

      # scenario_file() may return multiple paths (one per selected GCM).
      # Read each, subset to the active climate variables, and average
      # layer-wise to form a model ensemble. Subsetting *before* the sum
      # guarantees layers align by position == name across models, since
      # raster addition is positional. Reduce() (vs. do.call("+", ...))
      # is used so the single-model case is a clean no-op rather than an
      # error from calling `+` with one argument.
      fut_list <- lapply(sf, function(p) rast(p)[[input$clim_vars]])
      r_fut <- Reduce(`+`, fut_list) / length(fut_list)
      names(r_fut) <- input$clim_vars   # Reduce()/`+` can drop layer names

      list(hst = r_hst, fut = r_fut)
})

clim <- reactive({
      cp <- clim_pool()
      bb_ext <- ext(c(input$bbox_x, input$bbox_y))
      clim_hst <- cp$hst %>%
            crop(bb_ext) %>%
            aggregate(input$clim_agg, na.rm = TRUE)
      clim_fut <- cp$fut %>%
            crop(bb_ext) %>%
            aggregate(input$clim_agg, na.rm = TRUE)
      list(hst = clim_hst, fut = clim_fut)
})


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
      bb_ext <- ext(c(input$bbox_x, input$bbox_y))

      grid_of <- function(r) r %>% crop(bb_ext) %>%
            aggregate(input$clim_agg, na.rm = TRUE)

      if (!is.null(display_extra_vars[[v]])) {
            # A static covariate: the same raster for both timeframes,
            # since elevation does not change between them.
            g <- grid_of(rast(display_extra_vars[[v]]$file))
            names(g) <- v
            return(list(hst = g, fut = g))
      }

      # Averaged across the selected GCMs, matching clim_pool()'s ensemble.
      sf <- scenario_file()
      fut_list <- lapply(sf, function(pth) rast(pth)[[v]])
      r_fut <- Reduce(`+`, fut_list) / length(fut_list)
      names(r_fut) <- v

      list(hst = grid_of(rast("data/clim_historical.tif")[[v]]),
           fut = grid_of(r_fut))
})

## Same variable at the FIA plot coordinates, row-aligned with plot_occ().
cr_plot_values <- reactive({
      v <- input$cr_v1
      req(isTruthy(v), v %in% display_vars())
      xy <- as.data.frame(plot_occ()[, c("x", "y")])

      src <- if (!is.null(display_extra_vars[[v]])) {
            rast(display_extra_vars[[v]]$file)
      } else if (isTRUE(input$cr_clim == "future")) {
            sf <- scenario_file()
            fl <- lapply(sf, function(pth) rast(pth)[[v]])
            Reduce(`+`, fl) / length(fl)
      } else {
            rast("data/clim_historical.tif")[[v]]
      }
      as.numeric(terra::extract(src, xy, ID = FALSE, method = "bilinear")[, 1])
})

## Outline color for the focal polygon. selectizeInput with a swatch
## renderer rather than a text box: the swatches make the choice visual
## without adding a color-picker dependency, and `create = TRUE` still
## accepts any R color name or hex typed straight in.
output$polygon_color_ui <- renderUI({
      req(isTruthy(polygon_active()))
      # isolate: this renderUI CREATES polygon_color, so reading it
      # reactively would rebuild the widget on every change -- resetting
      # the dropdown and stealing focus mid-selection. The read is only to
      # carry the selection across a polygon swap.
      overlay_color_input("polygon_color", "Focal polygon color",
                          isolate(if (isTruthy(input$polygon_color)) input$polygon_color
                                  else "darkred"))
})

## Matching control for the focal site, shown only once one is set.
output$site_color_ui <- renderUI({
      req(!is.null(site$x), !is.null(site$y))
      overlay_color_input("site_color", "Focal site color",
                          isolate(if (isTruthy(input$site_color)) input$site_color
                                  else "#4c32a8"))
})

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
                         width = input$theta_geog * 3 * 1000)
      terra::mask(terra::crop(r, b), b)
}

## Is the x-axis variable one the model was actually fitted on? A
## relationship against anything else is a property of the landscape rather
## than of the model, and the plot says so rather than leaving it to be
## inferred.
cr_var_in_model <- reactive({
      isTruthy(input$cr_v1) && input$cr_v1 %in% input$clim_vars
})
