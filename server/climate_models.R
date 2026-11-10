## CLIMATE MODEL SELECTOR -------------------------------------

# Per-(scenario x era x model) climate deltas (future - historical) for
# the model-chooser module. Independent of input$clim_vars, the scenario
# selection, and the AIM pipeline: reads ALL climate variables straight
# from the raw rasters for *every* future file on disk, so the chart
# shows the full ensemble as time trajectories rather than a single
# selection snapshot. Extracting every variable (rather than just the
# current axis pair) lets the user re-pick chart axes with no recompute.
#
# Deltas are back-converted from standardized (z-score) units -- how the
# rasters are stored -- to native climate units by multiplying by each
# variable's SD; see clim_units / clim_scales.
#
# The region over which each delta is summarized is user-controlled
# (input$model_delta_region): the modeling-domain bbox, the focal site,
# or the uploaded study-area polygon.
#
# Cost: one raw future-raster read per future file -- the full ensemble,
# so this is heavy. It is therefore gated behind the "Update chart"
# button via bindEvent(): nothing recomputes until the user asks, even
# if the summary region changes in the meantime. Axis-variable changes
# never recompute -- all variables are already extracted, so the plot
# just reselects columns. bindCache() underneath still earns its place:
# re-clicking Update after returning to a previously-seen region is a
# cache hit rather than a re-read. (reactive + bindCache + bindEvent is
# the composable equivalent of an eventReactive that also caches.)

# Label used for the synthetic zero-delta "historical" anchor that roots
# every model trajectory at the origin.
.hist_era_label <- "1991-2010 (historical)"

model_deltas <- reactive({
      region <- input$model_delta_region
      req(isTruthy(region))

      shiny::validate(need(nrow(scenario_meta) > 0,
                           "No future climate files found in data/clim_future/."))

      if (region %in% c("cell", "site")) {
            shiny::validate(need(!is.null(site$x) && !is.null(site$y),
                                 "Click the map to select a focal site, or choose a different summary region."))
      } else if (region %in% c("polygon", "buffered")) {
            shiny::validate(need(isTruthy(polygon_upload_active()),
                                 "Upload a focal polygon, or choose a different summary region."))
      }

      # Extract ALL climate variables (not just the current axis pair),
      # so the user can re-pick chart axes without triggering a recompute
      # -- the plot just reselects columns from the cached result. One
      # extra cost is reading 6 raster layers per file instead of 2,
      # which is minor next to the per-file I/O.
      all_vars <- clim_vars()
      hst <- rast("data/clim_historical.tif")[[all_vars]]

      # Region geometry resolved once, outside the per-file loop. The
      # buffered polygon and the site neighborhood differ from their
      # unbuffered counterparts only slightly for a climate delta, but they
      # are offered anyway: a wide buffer around a small polygon can span a
      # real elevation gradient, and a user should not have to learn which
      # module offers which areas.
      bb_ext <- ext(c(input$bbox_x, input$bbox_y))
      region_vect <- switch(
            region,
            "polygon"  = terra::vect(st_transform(user_polygon(), crs(hst))),
            "buffered" = terra::vect(st_transform(user_polygon_buffered(), crs(hst))),
            "site"     = terra::buffer(
                  terra::vect(cbind(site$x, site$y), crs = "EPSG:4326"),
                  width = input$theta_geog * 3 * 1000),
            NULL)

      # Summarize a (future - historical) delta SpatRaster to a named
      # numeric vector, one value per climate variable.
      summarize_delta <- function(delta) {
            if (region == "cell") {
                  v <- terra::extract(delta, cbind(site$x, site$y),
                                      method = "bilinear")
                  stats::setNames(as.numeric(v[1, ]), names(delta))
            } else if (!is.null(region_vect)) {
                  v <- terra::extract(delta, region_vect, fun = "mean",
                                      na.rm = TRUE, ID = FALSE)
                  stats::setNames(as.numeric(v[1, ]), names(delta))
            } else {
                  terra::global(terra::crop(delta, bb_ext), "mean",
                                na.rm = TRUE)[, 1] %>%
                        stats::setNames(names(delta))
            }
      }

      # Per-variable SD for back-converting standardized deltas to native
      # units. mean cancels in a difference, so only SD is needed. A
      # missing scale entry falls back to 1 (delta stays standardized).
      var_sd <- vapply(all_vars, function(v) {
            s <- clim_scales[[v]]$sd
            if (is.null(s) || !is.finite(s)) 1 else s
      }, numeric(1))

      rows <- lapply(seq_len(nrow(scenario_meta)), function(i) {
            m <- scenario_meta[i, ]
            fut <- rast(file.path("data/clim_future", m$file))[[all_vars]]
            d <- summarize_delta(fut - hst)
            # Back-convert each delta from z-score units to native units.
            d <- d[all_vars] * var_sd
            tibble::tibble(scenario = m$ssp, scenario_label = m$ssp_label,
                           era = m$period, model = m$model,
                           model_label = m$model_label) %>%
                  bind_cols(tibble::as_tibble(as.list(d)))
      })
      futures <- dplyr::bind_rows(rows)

      # Synthetic historical anchor: one zero-delta row per
      # (scenario x model) trajectory, so every path starts at the
      # origin. Era label is the constant defined above.
      zero_delta <- stats::setNames(as.list(rep(0, length(all_vars))),
                                    all_vars)
      anchors <- futures %>%
            distinct(scenario, scenario_label, model, model_label) %>%
            mutate(era = .hist_era_label) %>%
            bind_cols(tibble::as_tibble(zero_delta)[rep(1, nrow(.)), ])

      bind_rows(anchors, futures) %>%
            # Era ordering: synthetic historical first, then future
            # periods chronologically (string sort works for "YYYY-YYYY").
            mutate(era = factor(era,
                                levels = c(.hist_era_label,
                                           sort(unique(futures$era)))))
}) %>%
      # theta_geog is part of the key because the buffered polygon and the
      # site neighborhood are both sized by it -- without it, changing the
      # bandwidth would return a cached delta for the old geometry.
      bindCache(input$model_delta_region,
                input$bbox_x, input$bbox_y,
                site$x, site$y,
                input$theta_geog,
                polygon_upload_active()) %>%
      bindEvent(input$model_delta_update)

output$model_delta_plot <- renderPlot({
      d <- model_deltas()

      # Before the first "Update chart" click, model_deltas() is NULL
      # (bindEvent gate). Show an instruction rather than a blank panel.
      if (is.null(d) || nrow(d) == 0) {
            return(
                  ggplot() +
                        annotate("text", x = 0, y = 0,
                                 label = "Choose a summary region, then click \u201cUpdate chart\u201d.",
                                 color = "gray40", size = 5) +
                        theme_void()
            )
      }

      # Axes: the per-variable delta columns selected by the dropdowns.
      # Changing axes re-renders from the cached data with no recompute.
      vx <- input$model_delta_x
      vy <- input$model_delta_y
      req(vx %in% names(d), vy %in% names(d))
      d$dx <- d[[vx]]
      d$dy <- d[[vy]]

      # Axis labels use the shared unit_lab() helper; deltas are already
      # back-converted from z-scores in model_deltas.

      # Selected models (current GCM ensemble) are drawn at full opacity;
      # non-selected trajectories are faded to 0.5 so the active ensemble
      # stands out without needing in-plot text labels.
      d$selected <- d$model %in% input$scenario_model

      # Shape palette: one glyph per era, generated to fit however many
      # eras exist (filled square / diamond / circle / triangle, then
      # cycling) rather than a hardcoded 4-value scale.
      eras <- levels(d$era)
      shape_vals <- stats::setNames(
            rep(c(15, 18, 16, 17), length.out = length(eras)), eras)

      ggplot(d, aes(dx, dy,
                    color = model_label,
                    group = interaction(scenario, model))) +
            geom_hline(yintercept = 0, color = "gray70", linewidth = .3) +
            geom_vline(xintercept = 0, color = "gray70", linewidth = .3) +
            geom_path(aes(alpha = selected)) +
            geom_point(aes(shape = era, size = era, alpha = selected)) +
            facet_wrap(~ scenario_label) +
            scale_shape_manual(values = shape_vals) +
            scale_size_manual(values = stats::setNames(
                  seq(1.5, 3.5, length.out = length(eras)), eras)) +
            scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5),
                               guide = "none") +
            coord_cartesian(clip = "off") +
            theme(legend.position = "top",
                  legend.direction = "vertical") +
            strip_style +
            labs(x = paste0("\u0394 ", vx, " (", unit_lab(vx), ")"),
                 y = paste0("\u0394 ", vy, " (", unit_lab(vy), ")"),
                 color = "Climate model", shape = "Era", size = "Era")
})


pool <- reactive({
      req(plot_occ())
      xy <- plot_occ()[, c("x", "y")]
      # Extract climate from the raw (uncropped, unaggregated) rasters at
      # each plot's actual location. This gives the climate the plot is
      # actually experiencing rather than the spatial average over a
      # coarsened cell, and means pool() doesn't depend on the modeling
      # grid (input$clim_agg) or the bbox crop.
      cp <- clim_pool()
      clim_hst <- terra::extract(cp$hst, xy, ID = FALSE, method = "bilinear")
      clim_fut <- terra::extract(cp$fut, xy, ID = FALSE, method = "bilinear")

      # # code used interactively during dev to determine optimal index res
      # tune_index_res(x = cbind(xy, clim_fut), pool = cbind(xy, clim_hst),
      #                env = kernel("gaussian", input$theta_clim, input$theta_clim*3),
      #                geog = kernel("gaussian", input$theta_geog, input$theta_geog*3),
      #                stat = "sum_weights")

      # Response(s) for AIM smoothing. Pooled-ratio variables pack the
      # numerator/denominator/occurrence pieces so aim() can form the
      # surface as occurrence-gated ratio-of-means; everything else uses a
      # single response column.
      rp <- ratio_pieces()
      vals <- if (is.null(rp)) {
            tibble::tibble(response = response_values(),
                           samp_method_cd = plot_occ()$samp_method_cd)
      } else {
            tibble::tibble(occ = as.integer(rp$denom > 0),
                           num = rp$num,
                           denom = rp$denom,
                           samp_method_cd = plot_occ()$samp_method_cd)
      }

      # The same plots, in the same rows, in AIM coordinates (see
      # aim_surfaces()), as a table for the FFT method and an index for
      # the exact one. index_hst/index_fut stay lon/lat for the focal-site
      # analog search.
      aim_pts <- data.frame(aim_xy(xy), clim_hst)

      list(index_hst = build_analog_index(cbind(xy, clim_hst), coord_type = "lonlat",
                                          geog_res_adj = 4),
           index_fut = build_analog_index(cbind(xy, clim_fut), coord_type = "lonlat",
                                          geog_res_adj = 4),
           aim_pts = aim_pts,
           aim_index = build_analog_index(aim_pts, coord_type = "projected",
                                          geog_res_adj = 4),
           clim_hst = clim_hst,
           clim_fut = clim_fut,
           # Active response(s) for AIM smoothing (built as `vals` above),
           # plus samp_method_cd carried through for plot styling.
           values = vals)
})
