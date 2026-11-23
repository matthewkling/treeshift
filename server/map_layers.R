# MAP LAYERS (continued): the raster and vector layers.


## Evidence layers: what the FIA plots behind each cell's estimate amount to.
## Support is the sum of the plots' weights, i.e. the number of exactly
## matching plots the evidence is worth; closeness is the typical weight of
## the plots that drive the estimate (support = effective number of plots x
## closeness); spread is how much those plots' values disagree. Layer names
## in aim()'s output.
evidence_layers <- c(support = "sum_weights", closeness = "effective_weight",
                     spread = "spread")

## Support is drawn in fixed classes rather than a continuous ramp. In
## climatically uniform regions (e.g. the Great Plains) it runs to hundreds of
## plots, which on a ramp stretched to the map's maximum squeezes everywhere
## else into the bottom colors. What matters for reading it is the low end,
## and the classes are log-spaced because an estimate's precision grows with
## the square root of the evidence: going from 2 to 20 plots helps about as
## much as going from 20 to 200. Fixed breaks also mean a color is the same
## amount of support on every map (at the same bandwidths); the top class is
## open-ended. Classes include their lower break.
support_breaks <- c(0, 1, 3, 10, 30, 100, Inf)
support_labels <- c("under 1", "1 to 3", "3 to 10", "10 to 30", "30 to 100", "100 or more")

## Why the spread layer can't be drawn, or NULL when it can. Shown under the
## raster menu; display_layer() draws nothing in these cases.
evidence_note <- reactive({
      if (!identical(rast_stat(), "spread")) return(NULL)
      if (multi_mode()) {
            "Analog spread describes the plots behind one species, forest type or group, so it isn't shown for community measures."
      } else if (!spread_computed()) {
            "Analog spread isn't shown for ratio variables such as fractions and proportions."
      } else if (identical(input$time, "delta")) {
            "Analog spread is shown for the historical or future timeframe, not for change."
      }
})

## Why a connectivity map can't be drawn (its shiny::validate() message), or
## NULL. The map's observer would otherwise drop the message silently.
connectivity_note <- reactive({
      if (!identical(rast_stat(), "connectivity")) return(NULL)
      tryCatch({ connectivity_layer(); NULL },
               shiny.silent.error = function(e) if (nzchar(conditionMessage(e))) conditionMessage(e))
})

output$stat_note <- renderUI({
      msg <- c(evidence_note(), connectivity_note())
      if (!is.null(msg)) div(class = "ts-notice", msg)
})


## Paired measures are their own quantity, not a delta of a per-timeframe
## one, so they must not pick up the diverging zero-centered delta palette
## even when input$time happens to be "delta". Expected gains and losses
## are strictly non-negative; turnover runs 0 (identical) to 1 (no
## overlap). None of them center on zero.
## The *_var reactives say what kind of FIA-variable layer is selected,
## whatever the map shows; the *_layer ones, whether the map shows it.
paired_var <- reactive({
      multi_mode() && isTruthy(input$comm_measure) && input$comm_measure %in% paired_measures
})
paired_layer <- reactive(isTruthy(rast_stat()) && rast_stat() == "FIA variable" && paired_var())

## Categorical layers: the value is an entity index, not a quantity.
## Everything that assumes continuity has to branch on this -- palette,
## reprojection method, scale transform, and the panels that plot value
## distributions.
categorical_var <- reactive({
      multi_mode() && isTruthy(input$comm_measure) && input$comm_measure %in% categorical_measures
})
categorical_layer <- reactive(isTruthy(rast_stat()) && rast_stat() == "FIA variable" && categorical_var())

## Category codes, labels and colors for the categorical layers. Derived
## from the raster's own factor levels so the markers cannot drift from the
## surface they sit on: one mapping, two consumers.
comm_categories <- reactive({
      req(categorical_layer())
      comm_category_result()[c("codes", "labels", "colors")]
})

## Layers bounded on [0, 1] that read better on a diverging palette than a
## sequential one, with the domain pinned to the full range so the
## midpoint stays put rather than floating with the observed max.
##
## Diversity is deliberately NOT here: it is exp(H) in effective-entity
## units, unbounded, and shares units with richness -- keeping the two on
## the same sequential palette is what makes them comparable.
bounded_ratio_var <- reactive({
      if (multi_mode()) return(isTruthy(input$comm_measure) &&
                                     input$comm_measure %in% c("evenness", "dominance"))
      isTruthy(model$fia_var) && grepl("fraction", model$fia_var)
})
bounded_ratio_layer <- reactive(isTruthy(rast_stat()) && rast_stat() == "FIA variable" && bounded_ratio_var())

## Layers whose NAs carry meaning -- a support mask or an entity not found
## within range -- rather than merely marking the edge of the analog
## domain. Filling these with 0 on a delta map would read as "no change"
## instead of "no information", and for a ratio 0 is a real value at one
## end of the scale.
masked_var <- reactive({
      if (multi_mode()) return(isTruthy(input$comm_measure) &&
                                     !identical(input$comm_measure, "richness"))
      isTruthy(model$fia_var) && grepl("fraction", model$fia_var)
})
masked_layer <- reactive(isTruthy(rast_stat()) && rast_stat() == "FIA variable" && masked_var())


display_layer <- reactive({
      req(focal_id())

      # In All mode the community reduction takes the place of the
      # per-entity response. The model diagnostics (support, closeness, focal
      # site similarity) are entity-independent and keep working
      # unchanged, so only the "FIA variable" branch is replaced -- and
      # aim() is not required, so All mode does not also pay for a
      # single-entity run it will not display.
      if (multi_mode() && rast_stat() == "FIA variable") return(composition_layer())

      req(aim())  # only run once we have data

      if(rast_stat() %in% names(evidence_layers)){
            req(is.null(evidence_note()))
            r <- ts_evidence_layer(aim(), evidence_layers[[rast_stat()]], input$time)
      }else if(rast_stat() == "focal site similarity"){
            r <- site_diff()
      }else if(rast_stat() == "connectivity"){
            r <- connectivity_layer()
      }else{ # rast_stat() == "FIA variable"
            req(isTruthy(model$fia_var))
            r <- ts_value_layer(aim(), model$fia_var, input$time, input$support_floor)
      }

      r
})

## The displayed layer's title parts (helpers/titles.R)
raster_title_parts <- reactive({
      layer_title(stat = rast_stat(), time = input$time, track = model$track,
                  focal = focal_label(), fia_var = model$fia_var,
                  comm_measure = input$comm_measure, multi = multi_mode(),
                  direction = input$analog_direction)
})

## The raster's legend title: subject and quantity, then muted details
raster_title <- reactive({
      parts <- raster_title_parts()
      lines <- c(parts$subject, parts$quantity, parts$notes)

      # When the layer is categorical the marker legend is suppressed --
      # the points reuse these same entity colors, so a second copy would
      # repeat it verbatim. That leaves this legend serving both layers,
      # which means it also has to say what the points are.
      if (categorical_layer() && !is.null(analog_vector_note()))
            lines <- c(lines, paste("points:", analog_vector_note()))

      legend_title(lines, n_main = length(parts$subject) + 1)
})

## How a raster layer is drawn: the values drawn (transformed, if a scale
## transform is chosen), the palette and its domain, and `color()`, which gives
## a value in the layer's own units its color on that scale, clamped to the
## ends of the scale. `stat` is the kind of layer (a rast_stat() value). Shared
## by the map and the Model sensitivity panel (server/sensitivity.R), so the panel's
## points take the colors the legend shows.
raster_scale_of <- function(r, stat) {
      fia <- stat == "FIA variable"
      categorical <- fia && categorical_var()
      raw_layer <- r

      # Scale transforms are meaningless on category codes, and
      # setValues() would drop the factor levels the legend reads from.
      # Fixed classes (support, connectivity) already do the job of a
      # transform.
      bins <- switch(stat,
                     # Reversed, so well-supported cells take the dark, cool end
                     support = list(breaks = support_breaks, labels = support_labels,
                                    palette = rev(scales::viridis_pal()(length(support_labels)))),
                     connectivity = list(breaks = cs_breaks, labels = cs_break_labels,
                                         palette = cs_break_colors))
      binned <- !is.null(bins)
      # Maps a transformed value back to the variable's units, for the legend
      # labels: colors follow the transformed scale, labels stay in the
      # original units. NULL while the drawn raster is untransformed.
      # `to_scale` puts any other value on the transformed scale.
      untransform <- NULL
      to_scale <- identity
      if (!categorical && !binned && input$rast_transform != "identity") {
            raw <- r[]
            trans <- switch(input$rast_transform,
                            "log10" = log10,
                            "sqrt" = sqrt,
                            "uniform" = function(x) scales::rescale(rank(x, na.last = "keep")))
            to_scale <- switch(input$rast_transform,
                               "uniform" = function(x) stats::ecdf(raw)(x),
                               trans)
            # Ranks have no formula to undo: a value's rescaled rank is its
            # percentile, so the label is the data's value at that percentile
            untransform <- switch(input$rast_transform,
                                  "log10" = function(x) 10^x,
                                  "sqrt" = function(x) x^2,
                                  "uniform" = function(x) stats::quantile(raw, x, na.rm = TRUE, names = FALSE))
            r <- setValues(r, trans(raw))
            r[!is.finite(r[])] <- NA
      }

      ## palettes ##
      proj_method <- "bilinear"
      mm <- minmax(r)
      dom <- if(input$rast_transform == "identity") c(0, mm[2]) else mm

      choose_palette <- function(x) if(input$rast_palette == "default") x else input$rast_palette

      if(categorical){
            # Nearest neighbour: bilinear interpolation of category codes
            # invents indices that mean nothing. Same reason current
            # direction uses it.
            proj_method <- "ngb"
            cc <- comm_categories()
            # `levels`, not `domain`: colorFactor() builds its levels with
            # sort(unique(domain)), so passing the codes as a domain pairs
            # color k with the k-th SMALLEST entity index while the legend
            # labels it by frequency rank. Every name shifts. Supplying
            # levels explicitly keeps the given order.
            rast_pal <- colorFactor(palette = cc$colors, levels = cc$codes,
                                    na.color = "transparent")
            dom <- cc$codes
      } else if(input$time == "delta" && fia && !paired_var()){
            # units are deltas, drawn untransformed
            r <- raw_layer
            untransform <- NULL
            to_scale <- identity
            # Fill non-analog NAs -- except where the NAs are a deliberate
            # mask (species not found within range, or below the support
            # floor). Filling those would read as "no change" on the
            # diverging scale rather than "no information". Only within the
            # grid's climate: outside the contiguous US nothing is estimated.
            if(!masked_var()) r[is.na(r) & !is.na(clim()$hst[[1]])] <- 0
            mm <- minmax(r)
            dom <- max(abs(mm)) * c(-1, 1)
            rast_pal <- colorNumeric(
                  palette = c("darkred", "orange", "gray90", "dodgerblue", "darkblue") %>% choose_palette(),
                  domain = dom,
                  na.color = "transparent"
            )
      } else if(binned){
            rast_pal <- colorBin(
                  palette = bins$palette %>% choose_palette(),
                  domain = NULL,
                  bins = bins$breaks,
                  na.color = "transparent"
            )
            dom <- range(bins$breaks[is.finite(bins$breaks)])
      } else if(stat %in% names(evidence_layers)){
            # Closeness has a fixed 0-1 scale, so its colors mean the same
            # thing on every map. Viridis sets these diagnostics apart from
            # the variable's own gray-to-green palette.
            if(stat == "closeness" && input$rast_transform == "identity") dom <- c(0, 1)
            rast_pal <- colorNumeric(
                  palette = "viridis" %>% choose_palette(),
                  domain = dom,
                  na.color = "transparent"
            )
      } else if(fia) {
            # historical or future estimate: domain is unit except for total BA
            if(bounded_ratio_var()){
                  # Bounded [0, 1] ratios -- the regeneration fractions
                  # (relict <-> colonizing), and community evenness and
                  # dominance -- default to RdBu. Fixing the domain to the
                  # full [0, 1] range keeps the midpoint put rather than
                  # floating with the observed max, so the same color means
                  # the same value across species and scenarios.
                  if(input$rast_transform == "identity") dom <- c(0, 1)
                  rast_pal <- colorNumeric(
                        palette = "RdBu" %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent"
                  )
            } else {
                  rast_pal <- colorNumeric(
                        palette = c("gray90", green) %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent"
                  )
            }
      } else if(stat == "focal site similarity") {
            rast_pal <- colorNumeric(
                  palette = "magma" %>% choose_palette(),
                  domain = dom,
                  na.color = "transparent",
                  reverse = TRUE
            )
      } else {
            rast_pal <- colorNumeric(
                  palette = c("gray90", "darkorchid4") %>% choose_palette(),
                  domain = dom,
                  na.color = "transparent"
            )
      }

      color <- function(x) {
            v <- to_scale(x)
            if (!categorical) v <- pmin(pmax(v, min(dom)), max(dom))
            rast_pal(v)
      }
      list(r = r, pal = rast_pal, dom = dom, untransform = untransform, bins = bins,
           binned = binned, proj_method = proj_method, color = color)
}

# The map's raster layer, as drawn
raster_scale <- reactive({
      r <- display_layer()
      req(r)
      raster_scale_of(r, rast_stat())
})

# The FIA variable's scale at the map's timeframe: the map's own while the map
# shows the variable, otherwise the one it would have
variable_scale <- reactive({
      if (identical(rast_stat(), "FIA variable")) return(raster_scale())
      req(aim(), isTruthy(model$fia_var))
      raster_scale_of(ts_value_layer(aim(), model$fia_var, input$time, input$support_floor),
                      "FIA variable")
})

observe({
      sc <- raster_scale()
      leafletProxy("map") %>%
            clearImages() %>%
            addRasterImage(
                  sc$r,
                  colors  = sc$pal,
                  opacity = input$rast_opacity,
                  method = sc$proj_method,
                  group = "AIM results",
                  maxBytes = 10 * 1024 * 1024 # increase limit to 10 MB
            )
})


## The FIA plot markers' values and colors, shared by the markers and their
## legend
point_vectors <- c("all FIA plots", "occupied FIA plots", "focal site analogs")
point_scale <- reactive({
      pd <- vector_points()

      # Regeneration fractions get the variable's full [0, 1] range so
      # 0.5 stays on the diverging palette's midpoint. Not in All
      # mode, where the fractions have no per-entity form and the
      # displayed value is a count.
      cat_layer <- categorical_layer()
      cc <- if (cat_layer) comm_categories() else NULL

      dom <- if(cat_layer) cc$codes
      else if(response_is_binary()) 0:1
      else if(bounded_ratio_layer()) c(0, 1)
      else seq_domain(pd$display)

      # Markers follow the raster palette only when the raster is
      # itself showing FIA response on a sequential scale (see
      # marker_follows_raster_palette); otherwise the canonical
      # gray90 -> green ramp, except for the regeneration fractions,
      # which default to the same diverging RdBu as their raster (on
      # the sequential ramp their 0 = relict end renders near-white).
      # colorNumeric accepts either a color vector or a palette name.
      marker_palette <- if(input$rast_palette != "default" &&
                           marker_follows_raster_palette()){
            input$rast_palette
      } else if(bounded_ratio_layer()){
            "RdBu"
      } else {
            c("gray90", green)
      }
      # Categorical layers key on entity codes, so the markers need
      # the raster's own colorFactor rather than a numeric ramp --
      # and its breaks are the categories themselves.
      plot_pal <- if(cat_layer){
            # levels, not domain -- see the raster palette above.
            colorFactor(palette = cc$colors, levels = cc$codes)
      } else {
            colorNumeric(palette = marker_palette, domain = dom)
      }
      list(pd = pd, pal = plot_pal, dom = dom, categorical = cat_layer)
})


### Legend ----------------------------
## One panel for the raster and the FIA plot markers (helpers/cartography.R).
## NULL for a section with nothing to key.

raster_legend <- reactive({
      if (!isTRUE(input$rast_opacity > 0)) return(NULL)
      sc <- raster_scale()
      op <- input$rast_opacity
      key <- if (categorical_layer()) {
            cc <- comm_categories()
            legend_swatches(cc$colors, cc$labels, op)
      } else if (sc$binned) {
            # Plain class labels ("100 or more", not "100 - Inf"), colored
            # from a value inside each class
            b <- head(sc$bins$breaks, -1)
            legend_swatches(sc$pal(b + c(diff(b) / 2, 1)), sc$bins$labels, op)
      } else legend_bar(sc$pal, sc$dom, sc$untransform, op,
                        cuts = if (!is.null(sc$untransform) && input$rast_transform == "log10")
                              log_ticks(sc$dom))
      legend_section(raster_title(), key)
})

point_legend <- reactive({
      if (!isTruthy(input$vector) || !input$vector %in% point_vectors) return(NULL)
      ps <- point_scale()
      # A categorical layer's markers take the raster's entity colors, which
      # its legend already keys
      if (ps$categorical) return(NULL)
      breaks <- pretty(seq(ps$dom[1], ps$dom[2], length.out = 100), n = 5)
      title <- legend_title(c("FIA plots", paste("observed historical", response_noun()),
                              analog_vector_note()),
                            n_main = 2)
      legend_section(title, legend_circles(ps$pal(breaks), legend_num(breaks)))
})

observe({
      # A section that can't be drawn right now (its layer is waiting on a
      # req(), say) is left out rather than holding up the other
      part <- function(f) tryCatch(f(), shiny.silent.error = function(e) NULL,
                                   validation = function(e) NULL)
      html <- c(part(raster_legend), part(point_legend))
      proxy <- leafletProxy("map") %>% removeControl("map_legend")
      if (length(html)) proxy %>%
            addControl(HTML(paste(html, collapse = "")), position = "topright",
                       layerId = "map_legend", className = "ts-legend leaflet-control")
})


### Vectors ----------------------------
observe({
      req(plot_occ())

      proxy <- leafletProxy("map")

      ## clear old layers ##
      proxy %>%
            # clearMarkers() %>%
            clearGroup("FIA plots") %>%
            clearGroup("current_flow") %>%
            clearGroup("current_field")

      if(input$vector %in% point_vectors){

            ps <- point_scale()
            pd <- ps$pd
            plot_pal <- ps$pal

            proxy <- proxy %>%
                  addCircles(
                        data = filter(pd, samp_method_cd == 2),
                        lng = ~x,
                        lat = ~y,
                        radius = point_radius()/2, # meters
                        stroke = FALSE, fillColor = "black",
                        group = "FIA plots"
                  ) %>%
                  addCircles(
                        # Drop plots with no defined value: the regeneration
                        # fractions are NA wherever the species is absent as
                        # both the stage and adults, and colorNumeric would
                        # otherwise paint them in its na.color (a mid gray
                        # that reads darker than a true zero). No-op for the
                        # other variables, which are 0 rather than NA when
                        # the focal entity is absent.
                        data = filter(pd, samp_method_cd == 1, !is.na(display)),
                        lng = ~x,
                        lat = ~y,
                        radius = point_radius(), # meters
                        fillColor = ~plot_pal(display),
                        fillOpacity = 1,
                        stroke = TRUE, color = "black", weight = 1,
                        group = "FIA plots"
                  )
      }

      # Net flow of the range shift, scaled so the strongest is 1
      flow <- function() {
            f <- cs_shift()[[c("flow_x", "flow_y")]]
            f / max(global(sqrt(f[[1]]^2 + f[[2]]^2), "max", na.rm = TRUE)$max)
      }

      if(input$vector == "current animation"){
            f <- flow()
            v_data <- format_velocity_data(f[[1]], f[[2]])
            proxy <- proxy %>%
                  addVelocity(
                        content = v_data,
                        options = velocityOptions(
                              colorScale = c("#000000", "#000000"),
                              lineWidth = 2,
                              velocityScale = 10^input$vector_scale * 1,
                              particleMultiplier = 0.005,
                              displayValues = FALSE,
                              velocityType = "Current",
                              emptyString = "No current data"
                        ),
                        group = "current_flow"
                  )
      }

      if(input$vector == "current field"){
            f <- flow()
            spokes <- make_spokes(f[[1]], f[[2]],
                                  subsample = 11 - model$clim_agg,
                                  arrow_scale = 10^input$vector_scale * .1)

            line_pal <- colorNumeric(palette = c("gray90", "black"), domain = sqrt(spokes$lines$mag))
            point_pal <- colorNumeric(palette = c("gray90", "black"), domain = sqrt(spokes$points$mag))

            proxy <- proxy %>%
                  addGlPolylines(data = spokes$lines, color = line_pal(sqrt(spokes$lines$mag)), weight = .25,
                                 group = "current_field") %>%
                  addCircles(
                        data = spokes$points %>% bind_cols(st_coordinates(spokes$points)),
                        lng = ~X, lat = ~Y,
                        radius = point_radius()/2 / 10^input$vector_scale,
                        fillOpacity = 1, stroke = FALSE, fillColor = point_pal(sqrt(spokes$points$mag)),
                        group = "current_field"
                  )
      }

})
