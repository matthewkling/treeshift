# CLIMATE RESPONSE PLOT --------------------------------
output$climate_response <- renderPlot({
      req(display_layer(), clim())

      req(isTruthy(input$cr_v1), input$cr_v1 %in% display_vars())
      v1 <- input$cr_v1

      cr_scope <- if (isTruthy(input$cr_scope)) input$cr_scope else "domain"
      scope_set <- switch(cr_scope,
                          "polygon"  = ,
                          "buffered" = isTruthy(polygon_active()),
                          "site"     = !is.null(site$x) && !is.null(site$y),
                          TRUE)
      shiny::validate(need(scope_set, eval_scope_missing_msg(cr_scope)))
      if (cr_scope %in% c("polygon", "buffered")) {
            shiny::validate(need(polygon_in_domain(cr_scope), polygon_outside_msg))
      }
      scope_label <- names(eval_scope_names)[eval_scope_names == cr_scope]

      # Climate for the x-axis. Historical by default; 'future' lets you
      # pair a future prediction with the future climate that produced it,
      # rather than reading it against current climate.
      use_fut <- isTRUE(input$cr_clim == "future")

      if(input$cr_data == "modeled grid cells"){
            cr <- if (use_fut) cr_grid()$fut else cr_grid()$hst
            d <- rast(c(y = display_layer(), clim1 = cr)) %>%
                  cr_scope_mask() %>%
                  as.data.frame()
      }

      if(input$cr_data == "FIA plots"){
            # display_values(), not response_values(): in community mode the
            # per-plot counterpart of the community measure -- observed
            # richness, say -- so this path describes the same quantity the
            # grid-cell path does. response_values() would silently plot the
            # focal species the entity picker still holds while hidden.
            #
            # It also covers the pooled-ratio variables, whose pool()$values
            # carries occ/num/denom rather than a `response` column.
            dv <- display_values()
            shiny::validate(need(
                  !is.null(dv),
                  "This community measure has no per-plot counterpart. Switch the data to 'modeled grid cells'."))
            cv <- cr_plot_values()
            req(length(cv) == length(dv))
            d <- data.frame(y = dv, clim1 = cv)

            # Same membership test the evaluation module uses, so "focal
            # polygon" selects the same plots in both.
            keep <- eval_scope_membership(plot_occ()[, c("x", "y")],
                                          which = cr_scope)[[cr_scope]]
            req(!is.null(keep), length(keep) == nrow(d))
            d <- d[keep, , drop = FALSE]
      }

      # Climate rasters are stored standardized, so back-convert to native
      # units (z * sd + mean) before plotting -- this puts the x-axis in
      # mm / degrees C. Falls back to leaving a variable standardized if it
      # has no scale entry.
      back_transform <- function(x, v) {
            s <- clim_scales[[v]]
            if (is.null(s) || !is.finite(s$sd)) return(x)
            x * s$sd + s$mean
      }
      d$clim1 <- back_transform(d$clim1, v1)

      # A categorical layer inverts the plot rather than adapting it. There
      # is no curve to fit through entity identities, but the question
      # transposes cleanly: instead of "how does the value change along this
      # gradient", it becomes "where along this gradient does each entity
      # win". Horizontal violins answer that -- one distribution of the
      # climate variable per entity, read down the page.
      if (categorical_layer()) {
            cc <- comm_categories()
            d$label <- cc$labels[match(d$y, cc$codes)]
            d <- d[!is.na(d$label) & is.finite(d$clim1), , drop = FALSE]

            shiny::validate(need(nrow(d) > 0,
                                 "No values to plot for this layer over the current area."))

            # Ordered by median position along the climate gradient, so the
            # violins step across the panel and the axis itself reads as the
            # gradient. Frequency order left the entities scattered, which
            # made the sequence -- the point of the plot -- something the
            # reader had to reconstruct.
            #
            # Median rather than mean: these distributions are often skewed
            # or bimodal for a species spanning two climate zones, and the
            # mean can land between the modes, in a range the entity does
            # not actually occupy.
            #
            # "other" stays pinned to the bottom regardless. It is a pooled
            # remainder with no meaningful position on the gradient, and
            # letting it sort among the entities would imply it had one.
            med <- d %>%
                  filter(label != "other") %>%
                  group_by(label) %>%
                  summarize(m = stats::median(clim1, na.rm = TRUE), .groups = "drop") %>%
                  arrange(m)
            lev <- c(if ("other" %in% d$label) "other", as.character(med$label))
            d$label <- factor(d$label, levels = lev)

            p_main <- ggplot(d, aes(clim1, label, fill = label)) +
                  geom_violin(scale = "width", color = NA) +
                  # Same entity colors as the map and the composition panel,
                  # so an entity is recognizable across all of them.
                  scale_fill_manual(values = stats::setNames(cc$colors, cc$labels),
                                    guide = "none") +
                  theme_bw() +
                  labs(x = paste0(input$cr_clim, " ", v1, " (", unit_lab(v1), ")"),
                       y = NULL,
                       subtitle = scope_label,
                       caption = if (!cr_var_in_model())
                             paste0(v1, " is not one of the variables the model was fitted on; this relationship comes from its correlation with those that were."))

            p_dens <- ggplot(d, aes(clim1)) +
                  geom_density(fill = "#2b6a4d", color = NA, alpha = .35) +
                  theme_void()

            return(p_dens / p_main + plot_layout(heights = c(1, 8)))
      }

      # Smooth on a link scale so fitted curves respect each variable's
      # natural bounds: logit for [0, 1] variables, log for non-negative
      # variables, identity for differences. Uses a GLM with a natural-spline
      # basis; quasi-families accept the continuous grid-cell values without
      # integer/proportion warnings. The FIA-plots dataset is the raw
      # historical per-plot variable, so its family depends on the variable
      # alone; only the modeled-grid dataset honors the raster stat and time.
      # In community mode the plotted quantity is the community measure,
      # not the per-entity variable: richness and diversity are unbounded
      # counts however they were built, while evenness, dominance, turnover
      # and the two change fractions live on [0, 1].
      bounded01 <- if (multi_mode())
            isTruthy(input$comm_measure) &&
            input$comm_measure %in% c("evenness", "dominance", "turnover",
                                      "colonization demand", "loss fraction")
      else grepl("presence|fraction", model$fia_var) ||
            model$fia_var == "proportion basal area"
      fam <- if (input$cr_data == "FIA plots") {
            if (bounded01) quasibinomial() else quasipoisson()
      } else if (input$stat == "connectivity" &&
                 isTRUE(input$cs_stat == "current direction")) {
            gaussian()                      # circular (0-360 deg): a log link is meaningless
      } else if (input$stat != "FIA variable") {
            quasipoisson()                  # ESS / similarity / connectivity magnitudes: non-negative
      } else if (input$time == "delta") {
            gaussian()                      # future - historical difference (can be negative)
      } else if (bounded01) {
            quasibinomial()                 # bounded [0, 1]
      } else {
            quasipoisson()                  # counts / basal area: non-negative
      }

      # Base the y-axis on the fitted lines, not the SE ribbons -- the
      # ribbons can blow up where a group's smooth extrapolates into a
      # sparse tail. Refit the same per-group smooth, take the range of the
      # fitted values, and clip the display there with coord_cartesian (the
      # geom_smooth fit itself still uses all the data; only the view is
      # bounded).
      dd <- d[is.finite(d$y) & is.finite(d$clim1), ]
      shiny::validate(need(
            nrow(dd) > 0,
            "No values to plot for this layer over the current area."))
      fits <- tryCatch(fitted(glm(y ~ splines::ns(clim1, 4), family = fam, data = dd)),
                       error = function(e) numeric(0))
      # range() of an all-NA vector returns c(Inf, -Inf), and handing that
      # to coord_cartesian() blanks the panel with no error -- the plot
      # renders, the axes render, and nothing is drawn.
      ylim <- if (length(fits) && any(is.finite(fits)))
            range(fits[is.finite(fits)]) else NULL
      if (!is.null(ylim) && (!all(is.finite(ylim)) || diff(ylim) == 0)) ylim <- NULL

      p_main <- ggplot(d, aes(clim1, y)) +
            geom_smooth(method = "glm",
                        formula = y ~ splines::ns(x, 4),
                        method.args = list(family = fam),
                        color = "#2b6a4d", fill = "#2b6a4d") +
            coord_cartesian(ylim = ylim) +
            # Label the climate axes with their timeframe (so historical vs
            # future is never ambiguous) and their native units.
            labs(x = paste0(input$cr_clim, " ", v1, " (", unit_lab(v1), ")"),
                 subtitle = scope_label,
                 # FIA-plots data is raw historical observation, so label it
                 # by the variable rather than the (stat/time-aware) raster title.
                 y = if (input$cr_data == "FIA plots")
                       paste("observed", if (multi_mode()) comm_label() else model$fia_var)
                 else raster_title_text(),
                 caption = if (!cr_var_in_model() && input$cr_data != "FIA plots")
                       paste0(v1, " is not one of the variables the model was fitted on; this relationship comes from its correlation with those that were."))

      # Marginal density along the x variable: shows where the data
      # actually is, which is what explains any wandering at the ends of the
      # smooth. It carries its own y-scale, so it does not interact with the
      # fitted-line ylim on the main panel.
      p_dens <- ggplot(d, aes(clim1)) +
            geom_density(fill = "#2b6a4d", color = NA, alpha = .35) +
            theme_void()

      # patchwork aligns the two panels on a shared x-axis; both are built
      # from the same clim1.
      p_dens / p_main + plot_layout(heights = c(1, 8))
})
