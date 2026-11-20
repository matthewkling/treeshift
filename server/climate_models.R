## CLIMATE FUTURES --------------------------------------------

# Projected climate change (future minus historical, in native units) for
# every scenario, period and model, averaged over the applied model bounding
# box. Read from the coarse cube of climate changes (data_prep/clim_prep.R),
# the same one the climate storylines rank models on (ts_storyline_model()),
# so the chart always shows the ranking a storyline used.

# Label for the zero-change "historical" anchor that roots every model
# trajectory at the origin
.hist_era_label <- "1981-2010 (historical)"

model_deltas <- reactive({
      shiny::validate(need(nrow(scenario_meta) > 0,
                           "No future climate files found in data/clim_future/."))
      cube <- terra::crop(rast(core_data$clim_deltas),
                          ext(c(model$bbox_x, model$bbox_y)), snap = "out")
      means <- terra::global(cube, "mean", na.rm = TRUE)$mean
      # Layers are named <ssp>_<period>_<model>_<var>
      parts <- stringr::str_split_fixed(names(cube), "_", 4)
      futures <- tibble::tibble(scenario = parts[, 1], era = parts[, 2], model = parts[, 3],
                                var = parts[, 4], value = means) %>%
            tidyr::pivot_wider(names_from = var, values_from = value) %>%
            # Only the scenarios the app offers
            inner_join(scenario_meta %>%
                             select(scenario = ssp, era = period, model,
                                    scenario_label = ssp_label, model_label),
                       by = c("scenario", "era", "model"))

      vars <- intersect(clim_vars(), names(futures))
      anchors <- futures %>%
            distinct(scenario, scenario_label, model, model_label) %>%
            mutate(era = .hist_era_label) %>%
            bind_cols(tibble::as_tibble(stats::setNames(as.list(rep(0, length(vars))), vars))[rep(1, nrow(.)), ])

      bind_rows(anchors, futures) %>%
            # Historical first, then future periods chronologically
            mutate(era = factor(era, levels = c(.hist_era_label, sort(unique(futures$era)))))
})

output$model_delta_plot <- renderPlot({
      d <- model_deltas()
      vx <- input$model_delta_x
      vy <- input$model_delta_y
      req(vx %in% names(d), vy %in% names(d))
      d$dx <- d[[vx]]
      d$dy <- d[[vy]]
      if (!isTRUE(input$model_delta_all_ssp)) d <- d[d$scenario == model$scenario_ssp, ]

      # One glyph per era, cycling if there are more than four
      eras <- levels(d$era)
      shape_vals <- stats::setNames(rep(c(15, 18, 16, 17), length.out = length(eras)), eras)
      size_vals <- stats::setNames(seq(1.5, 3.5, length.out = length(eras)), eras)

      # The applied model, pathway and period, marked by a larger black glyph
      # behind its point
      sel <- d[d$model %in% model$scenario_model & d$scenario == model$scenario_ssp &
                     d$era == model$scenario_period, ]

      # Custom convex hull polygon layer
      StatChull <- ggproto("StatChull", Stat,
                           compute_group = function(data, scales) {
                                 data[chull(data$x, data$y), , drop = FALSE]
                           },
                           required_aes = c("x", "y")
      )
      stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                             position = "identity", na.rm = FALSE, show.legend = NA,
                             inherit.aes = TRUE, ...) {
            layer(
                  stat = StatChull, data = data, mapping = mapping, geom = geom,
                  position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                  params = list(na.rm = na.rm, ...)
            )
      }

      ggplot(d, aes(dx, dy, color = model_label, group = interaction(scenario, model))) +
            geom_hline(yintercept = 0, color = "gray70", linewidth = .3) +
            geom_vline(xintercept = 0, color = "gray70", linewidth = .3) +
            stat_chull(aes(group = era), fill = "gray70", alpha = .5, color = NA, linewidth = .25) +
            geom_path() +
            geom_point(data = sel, aes(shape = era), size = size_vals[as.character(sel$era)] + 2,
                       color = "black", show.legend = FALSE) +
            geom_point(aes(shape = era, size = era)) +
            facet_wrap(~ scenario_label) +
            # Explicit breaks: the highlight layer would otherwise put its era first
            scale_shape_manual(values = shape_vals, breaks = eras) +
            scale_size_manual(values = size_vals) +
            coord_cartesian(clip = "off") +
            theme(legend.position = if (plot_is_wide(session, "model_delta_plot")) "right" else "bottom",
                  legend.direction = "vertical",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(alpha = .5)) +
            strip_style +
            labs(x = paste0("\u0394 ", vx, " (", unit_lab(vx), ")"),
                 y = paste0("\u0394 ", vy, " (", unit_lab(vy), ")"),
                 color = "Climate model", shape = "Era", size = "Era")
})


# The analog pool (R/aim.R): plot climates, the search indexes, and the
# outcome columns -- the single response or a pooled ratio's pieces, plus
# samp_method_cd for plot styling. A community run passes its per-entity
# responses to the surfaces separately.
pool <- reactive({
      po <- plot_occ()
      rp <- ratio_pieces()
      ts_pool(po, clim_pool(), ts_pool_values(po, if (is.null(rp)) response_values(), rp))
})
