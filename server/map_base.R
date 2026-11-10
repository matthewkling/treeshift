# MAP LAYERS ---------------------------------------

### Basemap -----------------------------------------------
# No reactive dependencies — renders once at startup. All subsequent
# map updates go through leafletProxy in dedicated observers.
output$map <- renderLeaflet({
      leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
            addProviderTiles(
                  providers$OpenTopoMap,
                  options = providerTileOptions(className = "basemap-gray")
            ) %>%

            # zoom to FIA extent
            fitBounds(lng1 = fia_xy_bounds$xmin,
                      lat1 = fia_xy_bounds$ymin,
                      lng2 = fia_xy_bounds$xmax,
                      lat2 = fia_xy_bounds$ymax) %>%

            addScaleBar(
                  position = "bottomleft",
                  options = scaleBarOptions(
                        metric = TRUE, # Show kilometers/meters
                        imperial = FALSE, # Hide miles/feet
                        updateWhenIdle = TRUE # Update only after the user stops zooming/panning
                  )
            )
})


### Focal site -----------------------------
# Update site coordinates on click only
observeEvent(input$map_click, {
      site$x <- input$map_click$lng
      site$y <- input$map_click$lat
})

# Clear the focal site. Sets coords back to NULL; downstream observers
# and reactives short-circuit on this, so the map marker, focal-site
# rings, density-plot vertical line, and analog plots all clear.
observeEvent(input$site_clear, {
      site$x <- NULL
      site$y <- NULL
})

# Draw focal site whenever relevant inputs change
observe({
      # Dependencies that should trigger redraw
      input$bbox_x
      input$bbox_y
      input$theta_geog
      input$show_site
      site_color()
      site$x
      site$y

      leafletProxy("map") %>%
            clearGroup("focal_site")

      if(input$show_site && !is.null(site$x)){
            r1 <- geo_circle(cbind(site$x, site$y), input$theta_geog)
            r3 <- geo_circle(cbind(site$x, site$y), input$theta_geog * 3)

            leafletProxy("map") %>%
                  addCircles(lng = site$x, lat = site$y,
                             radius = input$theta_geog / 10,
                             opacity = 1, fillOpacity = .5, color = site_color(),
                             group = "focal_site") %>%
                  addPolygons(lng = r1[,1], lat = r1[,2],
                              color = site_color(), weight = 2, fill = FALSE,
                              dashArray = "12, 5",
                              group = "focal_site") %>%
                  addPolygons(lng = r3[,1], lat = r3[,2],
                              color = site_color(), weight = 2, fill = FALSE,
                              group = "focal_site")
      }
})

# raster layer of euclidean climatic distance to focal site
site_diff <- reactive({
      req(site$x, site$y)

      s <- switch(input$analog_direction,
                  "Contemporary analogs" = site_analogs()$hst_site,
                  "Reverse analogs" = site_analogs()$fut_site,
                  "Forward analogs" = site_analogs()$hst_site)

      r <- switch(input$analog_direction,
                  "Contemporary analogs" = clim()$hst,
                  "Reverse analogs" = clim()$hst,
                  "Forward analogs" = clim()$fut)

      for(i in names(r)) r[[i]] <- abs(r[[i]] - s[[i]]) ^ 2
      r <- sqrt(app(r, sum))
      return(r)
})


### Raster -------------------------------
## When the point layer is the focal site's analogs, which analogs it shows
## depends on analog_direction -- a control in a right-sidebar panel, with
## nothing on the map to say which setting is active. This is that line.
## NULL when the qualifier does not apply, so callers can drop it.
analog_vector_note <- reactive({
      if (!isTruthy(input$vector) || input$vector != "focal site analogs") return(NULL)
      if (!isTruthy(input$analog_direction)) return(NULL)
      paste0("(", tolower(input$analog_direction), " of focal site)")
})

## Plain-text form of the legend title, for ggplot axis labels.
##
## raster_title() emits HTML because the leaflet legend renders it. Callers
## used to strip the "<br/>" it once contained; when that was rewritten to
## <div> blocks the strippers silently stopped matching and started printing
## raw markup as an axis label. Stripping tags generically means the two
## cannot drift apart again.
raster_title_text <- reactive({
      trimws(gsub("\\s+", " ", gsub("<[^>]*>", " ", raster_title())))
})

## "delta" is the internal value; the control and the legend both say
## "change", which also covers the paired measures -- they are change
## quantities without being deltas of anything.
time_label <- reactive({
      if (identical(input$time, "delta")) "change" else input$time
})

## What the evaluation is scoring, for axis and legend labels: the
## community measure in community mode, the variable otherwise.
eval_response_label <- reactive({
      if (multi_mode() && isTruthy(input$comm_measure)) comm_label()
      else tolower(input$fia_var)
})

## Display name for the community measure, with the life stage folded in
## rather than repeated as a separate variable line. "sapling richness"
## says everything "richness (sapling presence)" said, on one line.
comm_label <- reactive({
      req(isTruthy(input$comm_measure))
      stage <- switch(input$fia_var,
                      "seedling presence" = "seedling ",
                      "sapling presence"  = "sapling ",
                      "tree presence"     = "adult ",
                      "")
      paste0(stage, input$comm_measure)
})
