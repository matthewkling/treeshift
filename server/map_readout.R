## MAP READOUT ------------------------------------------
##
## The displayed layer's value under the cursor, in a box at the map's bottom
## left. The map sends the cursor's position as input$map_hover (output$map's
## onRender, server/map_base.R), or NULL when the cursor leaves the map; the
## box then shows the focal site's value instead, which is also what a tap on
## a touch screen, where there's no hover, gives.
##
## Values are read from display_layer(), the layer the map draws, so the
## readout can't disagree with the colors.

## The two timeframes behind a change map, so a readout can show where the
## change starts and ends. NULL for layers whose change isn't a difference of
## two maps: the paired community measures, connectivity's range shift, and
## the evidence layers' combined value.
timeframe_layers <- reactive({
      if (!identical(input$time, "delta") || !identical(rast_stat(), "FIA variable")) return(NULL)
      tf <- c("historical", "future")
      if (multi_mode()) {
            req(isTruthy(input$comm_measure))
            if (input$comm_measure %in% paired_measures) return(NULL)
            lapply(tf, function(t) ts_community_layer(multi_aim(), input$comm_measure, t,
                                                      input$support_floor))
      } else {
            req(aim(), isTruthy(model$fia_var))
            lapply(tf, function(t) ts_value_layer(aim(), model$fia_var, t, input$support_floor))
      }
})

# A value to two significant figures, so the readout doesn't claim more
# precision than the estimates have
fmt_readout <- function(v, signed = FALSE) {
      s <- format(signif(v, 2), big.mark = ",", scientific = FALSE, trim = TRUE)
      if (signed && v > 0) paste0("+", s) else s
}

# Why a cell inside the bounding box has no value, in the order the causes
# apply: no climate data (outside the United States, or water), no analog
# plots for the timeframe(s) shown (a climate unlike any within range), or
# masked by the support floor.
no_estimate_reason <- function(xy) {
      if (is.na(terra::extract(clim()$hst[[1]], xy)[1, 1]))
            return("No estimate: no climate data here")

      if (rast_stat() %in% c("FIA variable", "spread")) {
            surf <- if (multi_mode()) multi_aim() else aim()
            support <- function(tf) terra::extract(surf[[tf]][["sum_weights"]], xy)[1, 1]
            tfs <- switch(input$time, "historical" = "bsl", "future" = "fut", c("bsl", "fut"))
            none <- tfs[vapply(tfs, function(tf) !isTRUE(support(tf) > 0), logical(1))]
            if (length(none)) {
                  which <- if (length(none) == 2) "this cell's historical or projected future climate"
                  else if (none == "bsl") "this cell's historical climate"
                  else "this cell's projected future climate"
                  return(paste("No estimate: no FIA plots within range have a climate similar enough to", which))
            }
      }

      if (masked_layer()) "No estimate: too little here to measure (below the support floor)"
      else "No estimate here"
}

# The readout's lines at a point: a value, or why there isn't one
readout_lines <- function(x, y) {
      r <- display_layer()
      e <- terra::ext(r)
      if (x < e$xmin || x > e$xmax || y < e$ymin || y > e$ymax)
            return("Outside the model bounding box")

      xy <- cbind(x, y)
      v <- terra::extract(r, xy)[1, 1]
      if (length(v) != 1 || is.na(v)) return(no_estimate_reason(xy))

      if (categorical_layer()) {
            cc <- comm_categories()
            return(if (is.numeric(v)) cc$labels[match(v, cc$codes)] else as.character(v))
      }

      units <- if (identical(rast_stat(), "FIA variable") && !multi_mode())
            var_units(model$fia_var)
      with_units <- function(s) if (is.null(units)) s else paste(s, units)

      tl <- timeframe_layers()
      if (!is.null(tl)) {
            ends <- vapply(tl, function(l) terra::extract(l, xy)[1, 1], numeric(1))
            if (all(is.finite(ends)))
                  return(c(sprintf("historical %s → future %s",
                                   fmt_readout(ends[1]), fmt_readout(ends[2])),
                           with_units(paste("change", fmt_readout(v, signed = TRUE)))))
            return(with_units(paste("change", fmt_readout(v, signed = TRUE))))
      }

      with_units(fmt_readout(v))
}

output$map_readout <- renderUI({
      h <- input$map_hover
      if (!is.null(h)) {
            where <- NULL
            x <- h$lng
            y <- h$lat
      } else if (isTruthy(site$x) && isTruthy(site$y)) {
            where <- "Focal site"
            x <- site$x
            y <- site$y
      } else return(NULL)

      lines <- readout_lines(x, y)
      div(class = "ts-map-readout",
          if (!is.null(where)) div(class = "ts-map-readout-where", where),
          lapply(lines, div))
})
