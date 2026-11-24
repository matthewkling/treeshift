# Shortcuts (helpers/shortcuts.R): a list of common questions, each of which
# sets up the map to answer it.
#
# A shortcut that changes only display settings applies them in place. Any
# other reopens the page at the new link, the way Back and Forward do
# (server/url_state.R): the page is built from the link's settings and computes
# its results straight away, so the entity, variable and timeframe menus never
# have to be stepped through one at a time. The link is added to the browser's
# history, so Back undoes a shortcut.

output$shortcuts <- renderUI(shortcut_list(current_state()))

observeEvent(input$shortcut, {
      s <- current_state()
      new <- shortcut_state(input$shortcut, s)
      if (is.character(new)) return()

      changed <- names(new)[!vapply(names(new), function(f) identical(new[[f]], s[[f]]), logical(1))]
      if (length(changed) == 0) return()
      if (all(changed %in% shortcut_display_ids)) {
            for (id in changed) updateSelectInput(session, id, selected = new[[id]])
      } else {
            updateQueryString(url_encode(new, url_ctx), mode = "push")
            # The browser keeps what the link doesn't hold (ui/head.R)
            session$sendCustomMessage("ts_reload", TRUE)
      }
})
