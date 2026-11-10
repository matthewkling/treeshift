## FOCAL SITE PLOTS ----------------------------

eras <- reactive({
      switch(input$analog_direction,
             "Reverse analogs" = c("historical", "future"),
             "Forward analogs" = c("future", "historical"),
             "Contemporary analogs" = c("historical", "historical"))
})

### climate space --------------------------------
output$clim_plot <- renderPlot({
      req(site_analogs())
      a <- site_analogs()

      v1 <- input$clim_vars[1]
      v2 <- tail(input$clim_vars, 1) # if only one var selected, v2 == v1
      scl1 <- clim_scales[[v1]]
      scl2 <- clim_scales[[v2]]

      if(input$analog_direction == "Reverse analogs"){
            d <- a$rev_analogs
            d$v1 <- d[[paste0(v1, "_hst")]]
            d$v2 <- d[[paste0(v2, "_hst")]]
            s <- a$fut_site
      }else if(input$analog_direction == "Forward analogs"){
            d <- a$fwd_analogs
            d$v1 <- d[[paste0(v1, "_fut")]]
            d$v2 <- d[[paste0(v2, "_fut")]]
            s <- a$hst_site
      }else if(input$analog_direction == "Contemporary analogs"){
            d <- a$bsl_analogs
            d$v1 <- d[[paste0(v1, "_hst")]]
            d$v2 <- d[[paste0(v2, "_hst")]]
            s <- a$hst_site
      }else{
            stop("bug")
      }

      d$v1 <- d$v1 * scl1$sd + scl1$mean
      d$v2 <- d$v2 * scl2$sd + scl2$mean

      # display_values(), not response_values(): the analog SET does not
      # depend on the entity -- the same plots are near and climatically
      # similar whichever species is selected -- but the point coloring
      # does, and in community mode it should be the community measure, as
      # on the map markers. Indexed via analog_index because pool()$values
      # has no `response` column for the pooled-ratio variables.
      dv <- display_values()
      shiny::validate(need(
            !is.null(dv),
            "This community measure has no per-plot value, so the analog plots cannot be colored by it."))
      req(length(dv) >= max(d$analog_index, na.rm = TRUE))
      d$value <- dv[d$analog_index]
      d <- arrange(d, value) %>% mutate(panel = "Climate space")


      # search perimeters
      perims <- rbind(clim_circle(c(s[[v1]], s[[v2]]), input$theta_clim,
                                  c(scl1$mean, scl2$mean), c(scl1$sd, scl2$sd)) %>%
                            mutate(r = "b"),
                      clim_circle(c(s[[v1]], s[[v2]]), input$theta_clim * 3,
                                  c(scl1$mean, scl2$mean), c(scl1$sd, scl2$sd)) %>%
                            mutate(r = "a"))

      fs <- analog_fill(d)
      d <- fs$d

      p <- ggplot() +
            geom_path(data = perims,
                      aes(v1, v2, group = r, linetype = r),
                      color = site_color(), linewidth = .5, alpha = .5) +
            geom_point(data = filter(d, samp_method_cd == 2),
                       aes(v1, v2, fill = value),
                       color = "gray70", size = 1.5) +
            geom_point(data = filter(d, samp_method_cd == 1),
                       aes(v1, v2, fill = value,
                           shape = paste0("analog FIA plot,\n", eras()[1], " climate")),
                       color = "black", size = 3) +
            geom_point(aes(x = s[[v1]] * scl1$sd + scl1$mean,
                           y = s[[v2]] * scl2$sd + scl2$mean),
                       shape = 22,
                       color = site_color(), fill = site_color(), size = 4, alpha = .5) +
            fs$scale +
            scale_shape_manual(values = 21) +
            scale_linetype(guide = "none") +
            facet_wrap(~ panel) +
            theme(legend.position = "top") +
            strip_style +
            labs(y = v1,
                 x = v2,
                 shape = NULL)

      if(input$analog_direction != "Contemporary analogs"){
            h <- a$hst_site
            f <- a$fut_site
            p <- p +
                  annotate(geom = "segment",
                           x = h[[v1]] * scl1$sd + scl1$mean,
                           y = h[[v2]] * scl2$sd + scl2$mean,
                           xend = f[[v1]] * scl1$sd + scl1$mean,
                           yend = f[[v2]] * scl2$sd + scl2$mean,
                           arrow = arrow(angle = 15, type = "closed"),
                           color = site_color(), alpha = .5)
      }

      if(nrow(filter(d, env_dist <= input$theta_clim * 3)) == 0){
            p <- p +
                  annotate(geom = "text", color = "red",
                           label = "No FIA plots found\nwithin geographic & climatic range\nof focal site",
                           x = mean(perims$v1), y = mean(perims$v2))
      }

      p
})


### distance space --------------------------------


output$dist_plot <- renderPlot({
      req(site_analogs())

      d <- switch(input$analog_direction,
                  "Reverse analogs" = site_analogs()$rev_analogs,
                  "Forward analogs" = site_analogs()$fwd_analogs,
                  "Contemporary analogs" = site_analogs()$bsl_analogs)
      # display_values(), not response_values(): the analog SET does not
      # depend on the entity -- the same plots are near and climatically
      # similar whichever species is selected -- but the point coloring
      # does, and in community mode it should be the community measure, as
      # on the map markers. Indexed via analog_index because pool()$values
      # has no `response` column for the pooled-ratio variables.
      dv <- display_values()
      shiny::validate(need(
            !is.null(dv),
            "This community measure has no per-plot value, so the analog plots cannot be colored by it."))
      req(length(dv) >= max(d$analog_index, na.rm = TRUE))
      d$value <- dv[d$analog_index]
      d <- arrange(d, value) %>% mutate(panel = "Distance space")

      fs <- analog_fill(d)
      d <- fs$d

      ymax <- if(input$filter_clim) input$theta_clim * 3.3 else NA

      # Both axes here ARE the two distances, so the weight is a function
      # of position and can be drawn as a field. That says what the scheme
      # is regardless of which points are shown -- point size alone goes
      # silent when the displayed plots all fall outside the cutoff, which
      # is exactly when someone has turned the filter off to look at them.
      #
      xmax_g <- input$theta_geog * 3.3
      ymax_g <- if (is.na(ymax)) max(d$env_dist, na.rm = TRUE) else ymax
      wfield <- expand.grid(
            geog_dist = seq(0, xmax_g, length.out = 140),
            env_dist  = seq(0, ymax_g, length.out = 140))
      wfield$w <- analog_weight(wfield$geog_dist, wfield$env_dist)

      ggplot(d, aes(geog_dist, env_dist, fill = value)) +
            # Now the only representation of the weighting, so it spans a
            # wide value range rather than a light tint. Stopping short of
            # pure black and pure white at the ends keeps the black-outlined
            # points legible over the darkest corner and leaves the panel
            # background distinguishable from a zero-weight region.
            #
            # alpha rather than fill, because fill already carries the
            # response value and a second fill scale would need another
            # package.
            geom_raster(data = wfield, aes(geog_dist, env_dist, alpha = w),
                        fill = "black", inherit.aes = FALSE) +
            scale_alpha_continuous(range = c(0.05, 0.95), limits = c(0, 1),
                                   guide = "none") +
            geom_vline(xintercept = input$theta_geog, linetype = "longdash", color = site_color(), alpha = .5) +
            geom_vline(xintercept = input$theta_geog * 3, color = site_color(), alpha = .5) +
            geom_hline(yintercept = input$theta_clim, linetype = "longdash", color = site_color(), alpha = .5) +
            geom_hline(yintercept = input$theta_clim * 3, color = site_color(), alpha = .5) +
            scale_x_continuous(expand = c(0, 0), limits = c(0, input$theta_geog * 3.3)) +
            scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +

            geom_point(data = filter(d, samp_method_cd == 2), color = "gray70", size = 1.5) +
            geom_point(data = filter(d, samp_method_cd == 1), shape = 21, color = "black", size = 3) +

            geom_point(data = tibble::tibble(geog_dist = 0, env_dist = 0,
                                             shape_label = paste0("selected site,\n",
                                                                  eras()[2], " climate")),
                       aes(x = geog_dist, y = env_dist, shape = shape_label),
                       color = site_color(), fill = site_color(),
                       size = 10, alpha = .5, inherit.aes = FALSE) +

            fs$scale +
            scale_shape_manual(values = 22) +
            facet_wrap(~ panel) +
            theme(legend.position = "top") +
            strip_style +
            labs(y = "Climate difference from focal site (z-score)",
                 x = "Geographic distance to focal site (km)",
                 shape = NULL)
})

### forest composition --------------------------------
## Composition profile: what the community IS, historical and future.
##
## Informationally this is the same two numbers per entity as the
## scatterplot -- (hst, fut) versus (hst, change) -- but the encoding
## carries something neither of the plots above does: because the bands are
## stacked, the total height of each column is the COMMUNITY SUM. So
## per-entity change and aggregate change are legible in one figure, and a
## community that reshuffles without changing size looks different from one
## that shrinks.
##
## Built from the model stack rather than from the analog plots, so it
## follows the "Summarize over" control like the rest of the panel. The
## previous version was site-only and silently ignored it.
##
## Forward analogs are gone: "where does this site's climate go" is a
## question about the analog machinery, not about composition, and it has
## no counterpart in the AIM surfaces.
output$comm_composition <- renderPlot({
      shiny::validate(need(is.null(comm_profile_ready()), comm_profile_ready()))
      cb <- comm_bands()
      cols <- cb$cols

      # Stack positions; "other" is last in the band order so it ends up
      # on top rather than displacing the dominant entities' baseline.
      bands <- cb$bands %>%
            mutate(hst_top = cumsum(hst), hst_bot = hst_top - hst,
                   fut_top = cumsum(fut), fut_bot = fut_top - fut)

      # Smoothstep between the two columns, which is what gives the
      # ribbons their shape. Straight edges read as a chart of two bars;
      # the curve reads as flow.
      tt <- seq(0, 1, length.out = 60)
      sm <- tt * tt * (3 - 2 * tt)
      rib <- purrr::map_dfr(seq_len(nrow(bands)), function(i) {
            b <- bands[i, ]
            tibble::tibble(entity = b$entity, t = tt,
                           ymin = b$hst_bot + (b$fut_bot - b$hst_bot) * sm,
                           ymax = b$hst_top + (b$fut_top - b$hst_top) * sm)
      })
      rib$entity <- factor(rib$entity, levels = bands$entity)

      # Label side. Left is the default because the left column IS the
      # current-day composition, and reading it as a list is the
      # first-order thing this plot is for. Entities that grow move to the
      # right, where their band is thicker and the label has room -- an
      # arriving entity has almost no historical band to sit against.
      lab <- bands %>%
            mutate(side = ifelse(fut > hst, "right", "left"),
                   y_anchor = ifelse(side == "right",
                                     (fut_bot + fut_top) / 2,
                                     (hst_bot + hst_top) / 2))

      # Every entity gets a label, so they have to be pushed apart rather
      # than dropped. One pass per side: sort by anchor, then walk upward
      # enforcing a minimum gap. The gap is sized so the labels always fit
      # in the panel, and a leader line ties each back to its band.
      y_max <- max(bands$hst_top, bands$fut_top, 0, na.rm = TRUE)
      spread <- function(y) {
            if (!length(y)) return(y)
            gap <- y_max / max(26, length(y) + 1)
            o <- order(y)
            ys <- y[o]
            for (i in seq_along(ys)[-1]) ys[i] <- max(ys[i], ys[i - 1] + gap)
            # A long run can push the top label past the panel; slide the
            # whole stack down by the overshoot rather than let it clip.
            over <- max(0, ys[length(ys)] - y_max)
            ys <- ys - over
            out <- numeric(length(y))
            out[o] <- ys
            out
      }
      lab <- lab %>%
            group_by(side) %>%
            mutate(y_label = spread(y_anchor)) %>%
            ungroup() %>%
            mutate(x_anchor = ifelse(side == "right", 1, 0),
                   x_label  = ifelse(side == "right", 1.03, -0.03),
                   hjust    = ifelse(side == "right", 0, 1))

      ggplot(rib, aes(t, ymin = ymin, ymax = ymax, fill = entity, group = entity)) +
            geom_ribbon(color = "white", linewidth = .2) +
            # Leader lines, in the band's own color, so a displaced label
            # is still unambiguously attached to its band.
            geom_segment(data = lab, inherit.aes = FALSE,
                         aes(x = x_anchor, xend = x_label,
                             y = y_anchor, yend = y_label, color = entity),
                         linewidth = .3, alpha = .6) +
            geom_text(data = lab, inherit.aes = FALSE,
                      aes(x = x_label, y = y_label, label = entity,
                          hjust = hjust, color = entity),
                      size = 3.1) +
            scale_fill_manual(values = cols, guide = "none") +
            scale_color_manual(values = cols, guide = "none") +
            scale_x_continuous(breaks = c(0, 1),
                               labels = c("historical", "future"),
                               expand = ggplot2::expansion(mult = c(0.42, 0.42))) +
            theme_bw() +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor = element_blank()) +
            labs(x = NULL, y = tolower(input$fia_var))
})
