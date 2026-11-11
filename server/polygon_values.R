## POLYGON RASTER VALUES ------------------------

# Equal-area CRS for buffering in meters (NAD83 / Conus Albers).
# Adequate for any CONUS-scale buffer; would need to change for
# other regions.
.buffer_crs <- "EPSG:5070"

# Polygon buffered by theta_geog * 3 (km), the analog search's geographic
# cutoff. Features are unioned before buffering so multi-feature polygons
# produce a single dissolved buffer rather than per-feature buffers.
# Returned in WGS84 to match the user polygon's CRS so map and raster
# operations are consistent.
user_polygon_buffered <- reactive({
      req(user_polygon(), isTruthy(input$theta_geog))
      poly <- user_polygon()
      st_transform(poly, .buffer_crs) %>%
            st_union() %>%
            st_buffer(dist = input$theta_geog * 3 * 1000) %>%   # km -> m
            st_transform(4326)
})

## Raster values within the focal site's analog neighborhood. The site was
## the only area this module treated as a single cell rather than a
## distribution, which made it the odd one out against every other module
## and gave it a vertical line where the others get a curve. The single
## cell is still marked, so "where my site sits within its neighborhood"
## stays readable.
site_raster_values <- reactive({
      req(!is.null(site$x), !is.null(site$y), display_layer(),
          isTruthy(input$theta_geog))
      b <- terra::buffer(terra::vect(cbind(site$x, site$y), crs = "EPSG:4326"),
                         width = input$theta_geog * 3 * 1000)
      terra::extract(display_layer(), b, na.rm = TRUE) %>%
            pull(2) %>%
            na.omit()
})

polygon_raster_values <- reactive({
      req(user_polygon(), display_layer())

      poly <- user_polygon()
      r <- display_layer()

      # reproject polygon to match raster CRS
      poly_proj <- st_transform(poly, crs(r))

      terra::extract(r, vect(poly_proj), na.rm = TRUE) %>%
            pull(2) %>%  # column 1 is ID, column 2 is values
            na.omit()
})

polygon_buffered_raster_values <- reactive({
      req(user_polygon_buffered(), display_layer())

      poly <- user_polygon_buffered()
      r <- display_layer()

      poly_proj <- st_transform(poly, crs(r))

      terra::extract(r, vect(poly_proj), na.rm = TRUE) %>%
            pull(2) %>%
            na.omit()
})

output$polygon_density <- renderPlot({
      req(display_layer())

      # A categorical layer has no distribution to draw a density of --
      # entity codes are labels, not magnitudes. The question transposes
      # instead: what share of each area does each entity account for.
      # Same areas, same colors, a composition bar apiece.
      if (categorical_layer()) {
            cc <- comm_categories()
            share <- function(v) {
                  v <- v[!is.na(v)]
                  if (!length(v)) return(NULL)
                  tibble::tibble(label = cc$labels[match(v, cc$codes)]) %>%
                        filter(!is.na(label)) %>%
                        count(label) %>%
                        mutate(share = n / sum(n))
            }

            areas <- list("Model bounding box" = na.omit(values(display_layer(), mat = FALSE)))
            if (isTruthy(polygon_active())) {
                  areas[["Buffered polygon"]] <- polygon_buffered_raster_values()
                  areas[["Focal polygon"]]    <- polygon_raster_values()
            }
            if (!is.null(site$x) && !is.null(site$y)) {
                  areas[["Focal site neighborhood"]] <- site_raster_values()
                  # The single clicked cell, as its own bar. It is one
                  # entity so the bar is solid, which is the point: it says
                  # which entity holds the site, in the same colors as
                  # every other row.
                  areas[["Focal site"]] <- terra::extract(
                        display_layer(), cbind(site$x, site$y))[1, 1]
            }

            db <- bind_rows(lapply(names(areas), function(nm) {
                  sh <- share(areas[[nm]])
                  if (is.null(sh)) NULL else mutate(sh, source = nm)
            }))
            shiny::validate(need(nrow(db) > 0, "No values to summarize over these areas."))

            db$source <- factor(db$source, levels = names(areas))
            # Vertical bars, so the stacking order runs down the page in
            # the same direction as the legend beside it -- the eye can
            # match a band to its key entry by height alone. Reversed
            # levels put the first legend entry at the top of the stack.
            db$label <- factor(db$label, levels = rev(cc$labels))

            return(
                  ggplot(db, aes(source, share, fill = label)) +
                        geom_col(width = 0.7) +
                        scale_fill_manual(values = stats::setNames(cc$colors, cc$labels),
                                          breaks = cc$labels, name = NULL) +
                        scale_y_continuous(labels = scales::percent, expand = expansion(0, 0)) +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 20, hjust = 1),
                              plot.caption = element_text(hjust = 0),
                              plot.caption.position = "plot") +
                        labs(x = NULL, y = "share of area",
                             subtitle = raster_title_text(),
                             caption = if (isTruthy(polygon_active()) && !polygon_in_domain())
                                   stringr::str_wrap(polygon_outside_msg, 45))
            )
      }

      domain_vals <- na.omit(values(display_layer(), mat = FALSE))

      # Polygon-dependent values: NULL when no polygon is selected. The
      # corresponding curves (focal polygon, buffer) are then simply absent
      # from the plot. Check polygon_active() directly to avoid
      # triggering user_polygon()'s internal req() short-circuit.
      has_polygon <- isTruthy(polygon_active())
      # Outside the box the polygon curve has no values to draw; say why
      # rather than leaving it silently missing. (Its buffer can still reach
      # into the box, and is drawn if so.)
      polygon_outside <- has_polygon && !polygon_in_domain()
      poly_vals   <- if (has_polygon) polygon_raster_values() else NULL
      buffer_vals <- if (has_polygon) polygon_buffered_raster_values() else NULL

      has_site  <- !is.null(site$x) && !is.null(site$y)
      site_vals <- if (has_site) site_raster_values() else NULL

      # Follow the user's overlay colors, so a curve here and the boundary
      # it describes on the map are recognizably the same thing.
      study_color   <- polygon_color()
      buffer_color  <- lighten_color(polygon_color())
      domain_color  <- "gray60"
      # Follows the user's focal-site color, so the marker here and the
      # marker on the map are the same thing to look at.
      focal_color   <- site_color()

      # Focal site value: extract from the displayed raster at the click
      # location. NA if no site selected or if the point falls on a
      # non-modeled cell.
      focal_val <- if (!is.null(site$x) && !is.null(site$y)) {
            v <- terra::extract(display_layer(),
                                cbind(site$x, site$y),
                                method = "bilinear")
            v <- v[1, 1]
            if (is.finite(v)) v else NA_real_
      } else {
            NA_real_
      }

      # 'uniform' rank-transforms all three distributions against the
      # modeling-domain ECDF, so the three curves remain comparable on a
      # common [0, 1] axis. Done in data space (ggplot scale_x_*
      # transforms operate on the visual layer only, which would distort
      # density estimates by squishing them after they're computed).
      x_lab <- raster_title_text()
      if (input$density_transform == "uniform") {
            ecdf_ref <- ecdf(domain_vals)
            if (has_polygon) {
                  poly_vals   <- ecdf_ref(poly_vals)
                  buffer_vals <- ecdf_ref(buffer_vals)
            }
            if (has_site && length(site_vals)) site_vals <- ecdf_ref(site_vals)
            domain_vals <- ecdf_ref(domain_vals)
            if (!is.na(focal_val)) focal_val <- ecdf_ref(focal_val)
            x_lab <- paste0(x_lab, " (uniform rank)")
            d_bounds <- c(0, 1)
      } else {
            d_bounds <- display_bounds()
      }

      # Build the long data frame, including only the curves we have data
      # for.
      parts <- list(
            tibble::tibble(value = domain_vals, source = "Model bounding box")
      )
      if (has_polygon) {
            parts <- c(parts, list(
                  tibble::tibble(value = buffer_vals, source = "Buffered polygon"),
                  tibble::tibble(value = poly_vals,   source = "Focal polygon")
            ))
      }
      if (has_site && length(site_vals)) {
            parts <- c(parts, list(
                  tibble::tibble(value = site_vals, source = "Focal site neighborhood")
            ))
      }
      d <- bind_rows(parts) %>%
            mutate(source = factor(source, levels = c("Model bounding box",
                                                      "Buffered polygon",
                                                      "Focal polygon",
                                                      "Focal site neighborhood",
                                                      "Focal site")))

      pal <- c("Model bounding box"      = domain_color,
               "Buffered polygon"        = buffer_color,
               "Focal polygon"           = study_color,
               "Focal site neighborhood" = lighten_color(focal_color),
               "Focal site"              = focal_color)
      ltypes <- c("Model bounding box"      = "solid",
                  "Buffered polygon"        = "dashed",
                  "Focal polygon"           = "solid",
                  "Focal site neighborhood" = "dashed",
                  "Focal site"              = "solid")

      # Optional focal-site reference line. Skipped (NULL → no-op when
      # added to ggplot) if no focal site is selected or the cell value
      # isn't finite. Built as a data-driven geom_vline so the color and
      # linetype scales pick it up for the legend.
      focal_layer <- if (!is.na(focal_val)) {
            geom_vline(data = tibble::tibble(xint = focal_val,
                                             source = factor("Focal site",
                                                             levels = c("Model bounding box",
                                                                        "Buffered polygon",
                                                                        "Focal polygon",
                                                                        "Focal site neighborhood",
                                                                        "Focal site"))),
                       aes(xintercept = xint, color = source,
                           linetype = source),
                       linewidth = 0.6,
                       alpha = 0.7,
                       inherit.aes = FALSE)
      } else {
            NULL
      }

      # Legend break order: include Focal site only when one is active.
      legend_breaks <- c("Model bounding box", "Buffered polygon", "Focal polygon")
      if (has_site && length(site_vals))
            legend_breaks <- c(legend_breaks, "Focal site neighborhood")
      if (!is.na(focal_val)) legend_breaks <- c(legend_breaks, "Focal site")

      # ECDF only. The paired density and ECDF panels showed the same
      # comparison twice, and the density was the weaker of the two here:
      # what the module is for is reading off how one area sits relative to
      # another at a given value, which is the cumulative curve. Two panels
      # also halved the width available for either.
      p <- ggplot(d, aes(x = value, color = source, linetype = source)) +
            stat_ecdf(aes(group = source),
                      geom = "step", linewidth = 0.7, pad = FALSE) +
            focal_layer +
            scale_y_continuous(
                  breaks = c(0, 0.25, 0.5, 0.75, 1),
                  expand = expansion(0, 0)
            ) +
            scale_color_manual(values = pal,
                               breaks = legend_breaks) +
            scale_linetype_manual(values = ltypes,
                                  breaks = legend_breaks) +
            theme(legend.position = "top",
                  plot.caption = element_text(hjust = 0),
                  plot.caption.position = "plot") +
            labs(x = x_lab,
                 y = "cumulative share of area",
                 color = NULL, linetype = NULL,
                 caption = if (polygon_outside) stringr::str_wrap(polygon_outside_msg, 45))

      # X-axis limits from the active variable's natural bounds. Infinite
      # bounds become NA so ggplot auto-scales that side from data.
      x_limits <- ifelse(is.finite(d_bounds), d_bounds, NA_real_)

      if (input$density_transform == "sqrt") {
            p <- p + scale_x_sqrt(limits = x_limits,
                                  expand = expansion(0, 0))
      } else if (input$density_transform == "log10") {
            # log10 lower limit can't be 0 (undefined); let data set it.
            log_limits <- x_limits
            if (!is.na(log_limits[1]) && log_limits[1] <= 0) log_limits[1] <- NA
            p <- p + scale_x_log10(limits = log_limits,
                                   expand = expansion(0, 0))
      } else {
            p <- p + scale_x_continuous(limits = x_limits,
                                        expand = expansion(0, 0))
      }

      p
})
