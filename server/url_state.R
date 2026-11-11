# URL state (helpers/url_state.R): read the link this session was opened with,
# say what in it couldn't be used, and keep the address bar's link current.

# The link's settings. The UI was already built from them (app.R); the server
# reads them too, for what it sets up itself: the focal site and the polygon
# menu (data.R), and the support floor's units (community.R).
link <- url_decode(isolate(session$clientData$url_search), url_ctx)

if (length(link$problems)) {
      showNotification(
            tagList(tags$strong("Some settings in this link couldn't be used:"),
                    tags$ul(lapply(link$problems, tags$li))),
            duration = NULL, type = "warning")
}

# The current settings, as the codec's state
current_state <- reactive({
      req(input$track, input$metric, input$time)
      list(track = input$track, entity_scope = input$entity_scope,
           species = input$species, forest_type = input$forest_type,
           forest_group = input$forest_group, stage = input$stage,
           metric = input$metric, fia_var = input$fia_var, time = input$time,
           comm_measure = input$comm_measure, comm_top_n = input$comm_top_n,
           support_floor = input$support_floor,
           clim_agg = input$clim_agg, bbox_x = input$bbox_x, bbox_y = input$bbox_y,
           polygon_unit = if (is.null(polygon_active())) "" else polygon_active(),
           site = if (is.null(site$x)) NULL else c(site$x, site$y),
           clim_vars = input$clim_vars, scenario_ssp = input$scenario_ssp,
           scenario_period = input$scenario_period, scenario_model = input$scenario_model,
           theta_geog = input$theta_geog, theta_clim = input$theta_clim,
           cs_enable = input$cs_enable, cs_sources = input$cs_sources,
           cs_targets = input$cs_targets, cs_conductance_layer = input$cs_conductance_layer,
           cs_loss_layer = input$cs_loss_layer, cs_conductance_scale = input$cs_conductance_scale,
           cs_loss_scale = input$cs_loss_scale,
           stat = input$stat, cs_stat = input$cs_stat, vector = input$vector,
           analog_direction = input$analog_direction)
})

# The query string for the current settings. Debounced, so dragging a slider
# rewrites the address bar once it settles rather than on every step.
link_query <- debounce(reactive(url_encode(current_state(), url_ctx)), 400)

# A link ending in #about opens on the About tab, so documentation can be
# linked to directly. The tab isn't a setting, so it lives in the hash rather
# than the schema's query string.
if (startsWith(isolate(session$clientData$url_hash), "#about")) {
      nav_select("main_nav", "about")
}

# "replace", so settings changes don't pile up in the browser's back history.
# Replacing the query also replaces the hash, so the About tab's is rewritten
# with it.
observe({
      updateQueryString(paste0(link_query(), if (identical(input$main_nav, "about")) "#about"),
                        mode = "replace")
})

# The full link, for the settings record
link_url <- reactive({
      cd <- session$clientData
      paste0(cd$url_protocol, "//", cd$url_hostname,
             if (nzchar(cd$url_port)) paste0(":", cd$url_port),
             cd$url_pathname, link_query())
})

# The copy itself happens in the browser (ui/left_sidebar.R), which reports
# back whether the clipboard accepted it
observeEvent(input$copy_link_done, {
      if (identical(input$copy_link_done, "ok")) {
            showNotification("Link copied. It reopens these settings in this version of TreeShift.",
                             type = "message")
      } else {
            showNotification("The link couldn't be copied automatically. Copy it from your browser's address bar instead.",
                             type = "warning")
      }
})
