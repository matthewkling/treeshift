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

# The current settings, as the codec's state: the applied model settings
# (server/model_settings.R) and the current values of everything else. The
# link and the settings record describe what's on the map, never settings
# still waiting for the Run model button.
current_state <- reactive({
      req(model$track, model$metric, input$time)
      c(reactiveValuesToList(model)[model_setting_ids],
        list(time = input$time,
             comm_measure = input$comm_measure, comm_top_n = input$comm_top_n,
             support_floor = input$support_floor,
             polygon_unit = if (is.null(polygon_active())) "" else polygon_active(),
             site = if (is.null(site$x)) NULL else c(site$x, site$y),
             cs_enable = input$cs_enable,
             stat = input$stat, evidence_stat = input$evidence_stat, cs_stat = input$cs_stat, vector = input$vector,
             analog_direction = input$analog_direction))
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

# Browser history: each model update (server/model_settings.R sets
# `history_step`) adds an entry, so Back returns to the previous model run.
# Other changes -- display settings, polygon, site -- replace the current
# entry rather than piling up. Replacing the query also replaces the hash, so
# the About tab's is rewritten with it.
history_step <- FALSE
shown_query <- NULL
observe({
      q <- link_query()
      changed <- !identical(q, shown_query)
      push <- changed && history_step && !is.null(shown_query)
      if (changed) {
            history_step <<- FALSE
            # Each distinct link shown is a view (server/usage_log.R);
            # after_run marks the ones a model run produced
            log_event("view", query = q, after_run = push)
      }
      shown_query <<- q
      updateQueryString(paste0(q, if (identical(input$main_nav, "about")) "#about"),
                        mode = if (push) "push" else "replace")
})

# Back and Forward (reported by ui/head.R) reload the page at that entry's
# link, which rebuilds the controls and results from its settings. Moving
# between the About page's sections changes only the hash, not the settings,
# and doesn't reload.
observeEvent(input$history_nav, {
      if (!identical(input$history_nav, shown_query)) session$reload()
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
