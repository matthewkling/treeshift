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
            ) %>%

            # The cursor's position, for the value readout
            # (server/map_readout.R): at most every 120 ms while it moves,
            # always ending with where it stopped, and NULL when it leaves
            # the map. Also returns to the view saved before Back, Forward or
            # a shortcut reloaded the page (ui/head.R), once, in place of the
            # FIA extent.
            htmlwidgets::onRender("
function(el, x) {
  var map = this, last = 0, timer = null, pending = null;
  try {
    var view = JSON.parse(sessionStorage.getItem('ts-map-view'));
    sessionStorage.removeItem('ts-map-view');
    if (view) map.setView([view.lat, view.lng], view.zoom, {animate: false});
  } catch (e) {}
  function send() {
    last = Date.now();
    Shiny.setInputValue('map_hover', pending);
  }
  map.on('mousemove', function(e) {
    pending = {lat: e.latlng.lat, lng: e.latlng.lng};
    var wait = 120 - (Date.now() - last);
    if (wait <= 0) send();
    else if (!timer) timer = setTimeout(function() { timer = null; send(); }, wait);
  });
  map.on('mouseout', function() {
    clearTimeout(timer); timer = null; pending = null;
    Shiny.setInputValue('map_hover', null);
  });
}")
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
      model$bbox_x
      model$bbox_y
      model$theta_geog
      input$show_site
      site_color()
      site$x
      site$y

      leafletProxy("map") %>%
            clearGroup("focal_site")

      if(input$show_site && !is.null(site$x)){
            r1 <- geo_circle(cbind(site$x, site$y), model$theta_geog)
            r3 <- geo_circle(cbind(site$x, site$y), model$theta_geog * 3)

            leafletProxy("map") %>%
                  addCircles(lng = site$x, lat = site$y,
                             radius = model$theta_geog / 10,
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
      paste(tolower(input$analog_direction), "of focal site")
})

## The legend title as one line, for plot axes
raster_title_text <- reactive(title_phrase(raster_title_parts()))

## What the FIA plot values and the evaluation describe, as a noun phrase:
## the community measure in community mode ("sapling richness"), the variable
## otherwise ("basal area")
response_noun <- reactive({
      if (multi_mode()) {
            req(isTruthy(input$comm_measure))
            comm_noun(input$comm_measure, model$fia_var, model$track)$noun
      } else var_noun(model$fia_var)
})
