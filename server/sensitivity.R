## SENSITIVITY -----------------------------------------
##
## How much the focal site's estimate changes when one choice at a time is
## varied from the applied settings (R/sensitivity.R). Computed only while the
## panel is open, so clicking sites with it closed costs nothing.

site_sensitivity <- reactive({
      req(identical(input$modules_acc, "Model sensitivity"))
      shiny::validate(need(isTruthy(site$x) && isTruthy(site$y),
                           "Click the map to choose a focal site."))
      shiny::validate(need(!multi_mode(),
                           "Model sensitivity is available for single species, forest types and forest groups."))
      req(isTruthy(model$fia_var))
      scenario_file()   # a missing scenario shows its message
      settings <- lapply(setNames(nm = surfaces_setting_ids), function(id) model[[id]])
      ts_site_sensitivity(settings, core_data, c(site$x, site$y))
})

# Point sizes: coarser classes than the map's support legend, so the legend
# fits on one row below the plot. Classes include their lower break.
sensitivity_support_breaks <- c(0, 1, 10, 100, Inf)
sensitivity_support_labels <- c("<1", "1–10", "10–100", "100+")

# Row titles, top to bottom
sensitivity_rows <- c(gcm = "Climate model", ssp = "Emissions pathway",
                      theta_geog = "Geographic bandwidth", theta_clim = "Climate bandwidth",
                      clim_vars = "Climate variables", sample = "FIA plot sample")

# The plotted points: one per alternative, with its value at the map's
# timeframe, its label, its color on the map's scale and its support
sensitivity_points <- reactive({
      d <- site_sensitivity()
      time <- input$time

      # The support floor blanks a regeneration fraction where its species is
      # too unlikely to occur, as on the map
      if (grepl("fraction", model$fia_var) && isTruthy(input$support_floor)) {
            below <- function(occ) !is.na(occ) & occ < input$support_floor
            blank <- switch(time, historical = below(d$occ_hst), future = below(d$occ_fut),
                            delta = below(d$occ_hst) | below(d$occ_fut))
            shiny::validate(need(!blank[d$row == "sample"], paste(
                  "The estimate is blank at this site: the species' estimated occurrence is below",
                  "the minimum set under 'Map layers', so its regeneration fraction isn't shown there.")))
            d[blank, c("hst", "fut", "delta")] <- NA
      }

      col <- c(historical = "hst", future = "fut", delta = "delta")[[time]]
      d$value <- d[[col]]
      d$lower <- d[[paste0(col, "_lower")]]
      d$upper <- d[[paste0(col, "_upper")]]
      # Support for a change is the harmonic mean of its two ends', as on the
      # map's support layer: a change is only as well supported as its weaker end
      d$support <- switch(time, historical = d$support_hst, future = d$support_fut,
                          delta = 2 / (1 / d$support_hst + 1 / d$support_fut))
      d$support_class <- cut(d$support, sensitivity_support_breaks,
                             labels = sensitivity_support_labels,
                             right = FALSE, include.lowest = TRUE)
      shiny::validate(need(!is.na(d$value[d$row == "sample"]), paste(
            "No FIA plots fall within range of this site's",
            if (time == "delta") "historical or future" else time,
            "climate, so there is no estimate here.")))

      # Climate scenarios don't change a historical estimate
      if (time == "historical") d$value[d$row %in% c("gcm", "ssp")] <- NA

      # Labels: models (tagged with the storyline that picks them), pathways,
      # bandwidths with units, and the variable swapped in
      storyline_of <- tryCatch({
            r <- ts_storyline_ranking(core_data, model$scenario_ssp, model$scenario_period,
                                      model$bbox_x, model$bbox_y)
            n <- length(r)
            stats::setNames(c("least drying", "median drying", "most drying"),
                            names(r)[c(1, ceiling(n / 2), n)])
      }, error = function(e) character())
      lab <- d$alternative
      g <- d$row == "gcm"
      lab[g] <- toupper(d$alternative[g])
      tagged <- g & d$alternative %in% names(storyline_of)
      lab[tagged] <- paste0(lab[tagged], " (", storyline_of[d$alternative[tagged]], ")")
      p <- d$row == "ssp"
      lab[p] <- scenario_meta$ssp_label[match(d$alternative[p], scenario_meta$ssp)]
      lab[d$row == "theta_geog"] <- paste(d$alternative[d$row == "theta_geog"], "km")
      lab[d$row == "theta_clim"] <- paste(d$alternative[d$row == "theta_clim"], "z")
      cv <- d$row == "clim_vars"
      lab[cv] <- vapply(strsplit(d$alternative[cv], ", "), function(v) {
            out <- setdiff(model$clim_vars, v)
            if (!length(out)) paste(v, collapse = ", ") else paste(out, "→", setdiff(v, model$clim_vars))
      }, "")
      lab[d$row == "sample"] <- "current settings"
      d$label <- lab

      d$color <- variable_scale()$color(d$value)
      d$row <- factor(sensitivity_rows[d$row], levels = rev(sensitivity_rows))
      d
})

# The legend goes beside the plot when it is wide, and below it when it is
# narrow, where the plot grows to make room
sensitivity_wide <- function() plot_is_wide(session, "sensitivity_plot")

# "change in Pinus edulis presence probability", "Pinus edulis historical
# presence probability"
sensitivity_title <- function() {
      what <- paste(tolower(model$fia_var), fia_var_units(model$fia_var))
      if (input$time == "delta") paste("change in", focal_label(), what)
      else paste(focal_label(), input$time, what)
}

output$sensitivity_plot <- renderPlot(height = function() if (sensitivity_wide()) 380 else 460, {
      d <- sensitivity_points()
      wide <- sensitivity_wide()
      ref <- d$value[d$row == sensitivity_rows[["sample"]]]
      # Largest first, so smaller points land on top of them
      pts <- d[!is.na(d$value), ]
      pts <- pts[order(-pts$support), ]

      # Each row's range, as a line through its points
      spans <- pts[pts$row != sensitivity_rows[["sample"]], ] %>%
            group_by(row) %>%
            summarise(lo = min(value), hi = max(value), .groups = "drop")

      # Rows with no points: say why
      notes <- d %>% group_by(row) %>% filter(all(is.na(value))) %>% slice(1) %>% ungroup()
      notes$text <- if (input$time == "historical")
            ifelse(notes$row %in% sensitivity_rows[c("gcm", "ssp")],
                   "Doesn't apply to the historical timeframe", "No plots in range")
      else rep("No plots in range", nrow(notes))

      # The current settings' outline, on its points and as the reference
      # line; a point's stroke converted to the same width as a line
      current_stroke <- 1.3
      current_linewidth <- current_stroke * .stroke / 2 / .pt

      ggplot(d, aes(value, row)) +
            geom_blank() +
            # For a change, zero is always in view, so it's plain whether the
            # alternatives agree on its direction
            {if (input$time == "delta") list(geom_vline(xintercept = 0, color = "gray50", linewidth = .5),
                                             expand_limits(x = 0))} +
            geom_vline(xintercept = ref, color = "black", linewidth = current_linewidth) +
            geom_segment(data = spans, aes(x = lo, xend = hi, y = row, yend = row),
                         color = "gray50", linewidth = .5) +
            geom_errorbar(data = d[!is.na(d$lower), ], aes(xmin = lower, xmax = upper),
                          width = .25, color = "gray50", linewidth = .5, orientation = "y") +
            # Every point here, the current one included, so every class
            # present has a legend key; the current point is outlined on top
            geom_point(data = pts, aes(fill = color, size = support_class),
                       shape = 21, color = "gray40", stroke = .4) +
            geom_point(data = pts[pts$current, ], aes(fill = color, size = support_class),
                       shape = 21, color = "black", stroke = current_stroke, show.legend = FALSE) +
            {if (nrow(notes)) geom_text(data = notes, aes(x = -Inf, label = text), hjust = -0.05,
                                        size = 3, color = "gray40")} +
            scale_fill_identity() +
            # Sizes fixed per class, so a class keeps its size whichever
            # classes are present; the legend shows only those present
            scale_size_manual(values = stats::setNames(seq(4, 12, length.out = length(sensitivity_support_labels)),
                                                       sensitivity_support_labels),
                              name = "Analog support",
                              guide = if (wide) guide_legend(ncol = 1)
                                      else guide_legend(nrow = 1, title.position = "top")) +
            scale_y_discrete(expand = expansion(add = .6)) +
            coord_cartesian(clip = "off") +
            theme(legend.position = if (wide) "right" else "bottom",
                  legend.key.spacing.x = unit(2, "pt"),
                  panel.grid = element_line(color = "gray96"),
                  panel.grid.major.y = element_blank()) +
            # The row titles take about 40% of the width
            labs(x = wrap_to_plot(sensitivity_title(), session, "sensitivity_plot", share = 0.6),
                 y = NULL)
})

# Hover readout: every point under the cursor, with its setting and value
output$sensitivity_hover <- renderUI({
      h <- input$sensitivity_hover
      req(h)
      d <- sensitivity_points()
      near <- nearPoints(d[!is.na(d$value), ], h, xvar = "value", yvar = "row",
                         threshold = 12, maxpoints = 8)
      req(nrow(near) > 0)
      fmt <- function(x) format(signif(x, 3), big.mark = ",", trim = TRUE)
      lines <- sprintf("%s: %s%s (support %s)", near$label, fmt(near$value),
                       ifelse(is.na(near$lower), "",
                              sprintf(", 5–95%% range %s to %s", fmt(near$lower), fmt(near$upper))),
                       fmt(near$support))
      div(class = "ts-plot-tooltip",
          style = sprintf("left:%.0fpx; top:%.0fpx;", h$coords_css$x + 12, h$coords_css$y + 12),
          tags$strong(as.character(near$row[1])), br(),
          HTML(paste(htmltools::htmlEscape(lines), collapse = "<br>")))
})

# The evidence behind the current estimate, as the map's evidence layers
# measure it
output$sensitivity_evidence <- renderUI({
      d <- site_sensitivity()
      now <- d[d$row == "sample", ]
      fmt <- function(x) format(signif(x, 2), trim = TRUE)
      helpText(sprintf(paste(
            "At your current settings, the historical estimate rests on support of %s (closeness %s)",
            "and the future estimate on support of %s (closeness %s). Point size shows each alternative's support."),
            fmt(now$support_hst), fmt(now$closeness_hst), fmt(now$support_fut), fmt(now$closeness_fut)))
})
