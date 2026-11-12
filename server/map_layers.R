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

output$stat_note <- renderUI({
      msg <- evidence_note()
      if (!is.null(msg)) div(class = "ts-notice", msg)
})


## Paired measures are their own quantity, not a delta of a per-timeframe
## one, so they must not pick up the diverging zero-centered delta palette
## even when input$time happens to be "delta". Expected gains and losses
## are strictly non-negative; turnover runs 0 (identical) to 1 (no
## overlap). None of them center on zero.
paired_layer <- reactive({
      isTruthy(rast_stat()) && rast_stat() == "FIA variable" &&
            multi_mode() && isTruthy(input$comm_measure) &&
            input$comm_measure %in% paired_measures
})

## Categorical layers: the value is an entity index, not a quantity.
## Everything that assumes continuity has to branch on this -- palette,
## reprojection method, scale transform, and the panels that plot value
## distributions.
categorical_layer <- reactive({
      isTruthy(rast_stat()) && rast_stat() == "FIA variable" &&
            multi_mode() && isTruthy(input$comm_measure) &&
            input$comm_measure %in% categorical_measures
})

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
bounded_ratio_layer <- reactive({
      if (!isTruthy(rast_stat()) || rast_stat() != "FIA variable") return(FALSE)
      if (multi_mode()) return(isTruthy(input$comm_measure) &&
                                     input$comm_measure %in% c("evenness", "dominance"))
      isTruthy(model$fia_var) && grepl("fraction", model$fia_var)
})

## Layers whose NAs carry meaning -- a support mask or an entity not found
## within range -- rather than merely marking the edge of the analog
## domain. Filling these with 0 on a delta map would read as "no change"
## instead of "no information", and for a ratio 0 is a real value at one
## end of the scale.
masked_layer <- reactive({
      if (!isTruthy(rast_stat()) || rast_stat() != "FIA variable") return(FALSE)
      if (multi_mode()) return(isTruthy(input$comm_measure) &&
                                     !identical(input$comm_measure, "richness"))
      isTruthy(model$fia_var) && grepl("fraction", model$fia_var)
})


display_layer <- reactive({
      req(focal_id())

      # In All mode the community reduction takes the place of the
      # per-entity response. The model diagnostics (support, closeness, focal
      # site similarity, connectivity) are entity-independent and keep working
      # unchanged, so only the "FIA variable" branch is replaced -- and
      # aim() is not required, so All mode does not also pay for a
      # single-entity run it will not display.
      if (multi_mode() && rast_stat() == "FIA variable") return(composition_layer())

      req(aim())  # only run once we have data

      if(rast_stat() %in% names(evidence_layers)){
            req(is.null(evidence_note()))
            layer <- evidence_layers[[rast_stat()]]
            # For change, the harmonic mean of the two timeframes, which leans
            # toward the weaker one: a change is only as well supported as its
            # less supported end.
            r <- switch(input$time,
                        "delta" = 2 / (1/aim()$fut[[layer]] + 1/aim()$bsl[[layer]]),
                        "historical" = aim()$bsl[[layer]],
                        "future" = aim()$fut[[layer]])
      }else if(rast_stat() == "focal site similarity"){
            r <- site_diff()
      }else if(rast_stat() == "connectivity"){
            st <- cs_stats[cs_stat_names == input$cs_stat]
            r <- connectivity()[[st]]
      }else{ # rast_stat() == "FIA variable"
            # For the pooled-ratio variables, aim() carries the occurrence
            # layer alongside the response; it is the denominator these
            # fractions are conditional on. Mask per timeframe, so a delta
            # is blank wherever either end is unsupported.
            is_frac <- isTruthy(model$fia_var) && grepl("fraction", model$fia_var)
            floor_at <- function(tf) {
                  x <- tf$response
                  if (is_frac && !is.null(tf$occ))
                        x <- apply_support_floor(x, tf$occ, input$support_floor)
                  x
            }
            r <- switch(input$time,
                        "delta" = floor_at(aim()$fut) - floor_at(aim()$bsl),
                        "historical" = floor_at(aim()$bsl),
                        "future" = floor_at(aim()$fut))
      }

      r
})

raster_title <- reactive({
      ## Legend title. Each semantic unit is its own block, so it starts on
      ## a new line but WRAPS internally rather than being split at a
      ## hard-coded point. The old version broke "presence probability"
      ## after the first word regardless of available width; now that the
      ## title sits in a flexible column beside the color bar, forced
      ## breaks fight the layout instead of helping it.
      lines <- character()
      muted <- logical()
      add <- function(txt, dim = FALSE){
            if (!length(txt) || is.na(txt) || !nzchar(txt)) return(invisible(NULL))
            lines <<- c(lines, txt)
            muted <<- c(muted, dim)
      }

      if(rast_stat() %in% c("support", "closeness")){
            add(names(evidence_choices)[evidence_choices == rast_stat()])
            add(time_label())
            add(if (rast_stat() == "support") "(equivalent number of exactly matching plots)"
                else "(1 = exact match)", dim = TRUE)

      }else if(rast_stat() == "focal site similarity"){
            add(gsub(" analogs", "", input$analog_direction))
            add("climate distance from focal site")

      }else if(multi_mode() && rast_stat() == "FIA variable"){
            add(focal_label())
            add(time_label())
            add(comm_label())

            # Which per-entity variable the measure was computed FROM.
            # Only richness implies it -- richness is defined over
            # occurrence and nothing else, and comm_label() already folds
            # the life stage into its name. Every other measure can be
            # built from any variable, so "diversity" alone names two
            # different maps: the diversity of occurrence probabilities and
            # the diversity of basal area.
            if (!identical(input$comm_measure, "richness"))
                  add(paste0("of ", tolower(model$fia_var)), dim = TRUE)

            # The measures that keep the variable's units. Most reduce to a
            # dimensionless quantity -- diversity and richness are effective
            # entity counts, and evenness, dominance, turnover and the two
            # change fractions are ratios in which the units cancel. The
            # expected gains and losses are sums of per-entity change, so
            # they are still in whatever the variable measures.
            if (input$comm_measure %in% c("expected gains", "expected losses"))
                  add(switch(model$fia_var,
                             "total basal area"    = "(ft2 BA / acre)",
                             "number of trees"     = "(trees / acre)",
                             "number of saplings"  = "(saplings / acre)",
                             "number of seedlings" = "(seedlings / acre)",
                             # Summed change in occurrence probability is a
                             # count of entities, like richness.
                             if (identical(model$track, "Species")) "(species)"
                             else "(classes)"),
                      dim = TRUE)

      }else{
            if(rast_stat() == "spread") add("analog spread")
            add(focal_label())
            add(time_label())
            add(tolower(model$fia_var))
            add(switch(model$fia_var,
                       "total basal area"    = "(ft2 BA / acre)",
                       "number of trees"     = "(trees / acre)",
                       "number of saplings"  = "(saplings / acre)",
                       "number of seedlings" = "(seedlings / acre)",
                       NULL),
                dim = TRUE)
            if(rast_stat() == "connectivity") add(input$cs_stat)
      }

      # When the layer is categorical the marker legend is suppressed --
      # the points reuse these same entity colors, so a second copy would
      # repeat it verbatim. That leaves this legend serving both layers,
      # which means it also has to say what the points are.
      if (categorical_layer() && !is.null(analog_vector_note())) {
            add(paste("points:", analog_vector_note()), dim = TRUE)
      }

      # Support qualifier. Stated in the legend rather than only in a
      # tooltip because the legend is what ends up in a screenshot in
      # someone's planning document -- and "expected richness 8" over a
      # 10 km cell will otherwise be read as eight species in 100 km2,
      # which is wrong by a large factor rather than a subtlety.
      #
      # Applied to spread too: it is the variation among plot values. Not
      # applied to the other model diagnostics: support and closeness
      # describe the weights, focal site similarity is a climate distance,
      # connectivity is a circuit quantity.
      # In All mode this applies whatever the underlying variable: every
      # community reduction describes the composition of a plot-sized
      # footprint, and diversity in particular rises with area regardless
      # of whether it was built from occurrence or basal area.
      if (rast_stat() %in% c("FIA variable", "spread") &&
          (multi_mode() ||
           (isTruthy(model$fia_var) && model$fia_var %in% support_dependent_vars))) {
            add("(per ~0.07 ha plot)", dim = TRUE)
      }

      paste0("<div", ifelse(muted, " style='font-weight:400; opacity:.75;'", ""),
             ">", lines, "</div>", collapse = "")
})

observe({
      r <- display_layer()
      req(r)

      # Scale transforms are meaningless on category codes, and
      # setValues() would drop the factor levels the legend reads from.
      # Support's fixed classes already do the job of a transform.
      binned <- identical(rast_stat(), "support")
      # Maps a transformed value back to the variable's units, for the legend
      # labels: colors follow the transformed scale, labels stay in the
      # original units. NULL while the drawn raster is untransformed.
      untransform <- NULL
      if (!categorical_layer() && !binned && input$rast_transform != "identity") {
            raw <- r[]
            trans <- switch(input$rast_transform,
                            "log10" = log10,
                            "sqrt" = sqrt,
                            "uniform" = function(x) scales::rescale(rank(x, na.last = "keep")))
            # Ranks have no formula to undo: a value's rescaled rank is its
            # percentile, so the label is the data's value at that percentile
            untransform <- switch(input$rast_transform,
                                  "log10" = function(x) 10^x,
                                  "sqrt" = function(x) x^2,
                                  "uniform" = function(x) stats::quantile(raw, x, na.rm = TRUE, names = FALSE))
            r <- setValues(r, trans(raw))
            r[!is.finite(r[])] <- NA
      }

      proxy <- leafletProxy("map") %>%
            clearImages() %>%
            removeControl("raster_legend")

      ## palettes ##
      proj_method <- "bilinear"
      mm <- minmax(r)
      dom <- if(input$rast_transform == "identity") c(0, mm[2]) else mm

      choose_palette <- function(x) if(input$rast_palette == "default") x else input$rast_palette

      if(categorical_layer()){
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
      } else if(input$time == "delta" && rast_stat() == "FIA variable" && !paired_layer()){
            # units are deltas, drawn untransformed
            r <- display_layer()
            untransform <- NULL
            # Fill non-analog NAs -- except where the NAs are a deliberate
            # mask (species not found within range, or below the support
            # floor). Filling those would read as "no change" on the
            # diverging scale rather than "no information".
            if(!masked_layer()) r[is.na(r[])] <- 0
            mm <- minmax(r)
            dom <- max(abs(mm)) * c(-1, 1)
            rast_pal <- colorNumeric(
                  palette = c("darkred", "orange", "gray90", "dodgerblue", "darkblue") %>% choose_palette(),
                  domain = dom,
                  na.color = "transparent"
            )
      } else if(binned){
            rast_pal <- colorBin(
                  palette = "viridis" %>% choose_palette(),
                  domain = NULL,
                  bins = support_breaks,
                  na.color = "transparent"
            )
      } else if(rast_stat() %in% names(evidence_layers)){
            # Closeness has a fixed 0-1 scale, so its colors mean the same
            # thing on every map. Viridis sets these diagnostics apart from
            # the variable's own gray-to-green palette.
            if(rast_stat() == "closeness" && input$rast_transform == "identity") dom <- c(0, 1)
            rast_pal <- colorNumeric(
                  palette = "viridis" %>% choose_palette(),
                  domain = dom,
                  na.color = "transparent"
            )
      } else if(rast_stat() == "FIA variable") {
            # historical or future estimate: domain is unit except for total BA
            if(bounded_ratio_layer()){
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
      } else if(rast_stat() == "connectivity" && input$cs_stat == "current direction") {
            proj_method <- "ngb" # because bilinear interpolation breaks for a circular variable
            dom <- if(input$rast_transform == "identity") c(0, 360) else mm
            rast_pal <- colorNumeric(
                  palette = rainbow(20)[c(1:20, 1)],
                  domain = dom,
                  na.color = "transparent"
            )
      } else if(rast_stat() == "focal site similarity") {
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

      proxy <- proxy %>%
            addRasterImage(
                  r,
                  colors  = rast_pal,
                  opacity = input$rast_opacity,
                  method = proj_method,
                  group = "AIM results",
                  maxBytes = 10 * 1024 * 1024 # increase limit to 10 MB
            )

      if(input$rast_opacity > 0){
            # addLegend() derives its labels from `values` whenever `pal`
            # is supplied, and silently ignores a `labels` argument -- which
            # is why the categorical legend has to use the manual
            # colors/labels form instead. Entity names, not code numbers.
            proxy <- if (categorical_layer()) {
                  cc <- comm_categories()
                  proxy %>%
                        addLegend(
                              position = "topright",
                              colors   = cc$colors,
                              labels   = cc$labels,
                              title    = HTML(raster_title()),
                              layerId  = "raster_legend",
                              opacity  = input$rast_opacity,
                              className = "info legend legend-discrete"
                        )
            } else if (binned) {
                  # Manual form, for plain labels ("100 or more", not
                  # "100 - Inf"); colors read from the palette at a value
                  # inside each class
                  inside <- head(support_breaks, -1) + c(diff(head(support_breaks, -1)) / 2, 1)
                  proxy %>%
                        addLegend(
                              position = "topright",
                              colors   = rast_pal(inside),
                              labels   = support_labels,
                              title    = HTML(raster_title()),
                              layerId  = "raster_legend",
                              opacity  = input$rast_opacity,
                              className = "info legend legend-discrete"
                        )
            } else {
                  # Labels in the variable's units under a transform, each
                  # to 3 significant figures (leaflet's labelFormat() gives
                  # every label the same number of decimals: "1.000", "31.623")
                  lab_format <- if (is.null(untransform)) labelFormat() else
                        function(type, cuts, ...) {
                              format_one <- function(v) format(signif(v, 3), big.mark = ",",
                                                               scientific = FALSE, trim = TRUE)
                              vapply(untransform(cuts), format_one, character(1))
                        }
                  proxy %>%
                        addLegend(
                              position = "topright",
                              pal      = rast_pal,
                              values   = dom,
                              labFormat = lab_format,
                              title    = HTML(raster_title()),
                              layerId  = "raster_legend",
                              opacity  = input$rast_opacity
                        )
            }
      }
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
            clearGroup("current_field") %>%
            removeControl("point_legend")

      if(input$vector %in% c("all FIA plots", "occupied FIA plots", "focal site analogs")){

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
            # Only feeds the marker legend, which categorical layers do
            # not draw -- guarded so pretty() is not handed a code vector.
            breaks <- if(cat_layer) numeric(0)
            else pretty(seq(dom[1], dom[2], length.out = 100), n = 5)

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

            # legend
            breaks <- pretty(seq(dom[1], dom[2], length.out = 100), n = 5)
            legend_labels <- c(breaks, "Non-sampled plot")
            legend_colors <- c(plot_pal(breaks), "gray") # gray, to match the semi-transparent black markers
            legend_radii <- c(rep(6, length(breaks)), 3)
            # Categorical layers: markers use the raster's own colors and
            # categories, so a second legend would repeat it verbatim. The
            # raster legend serves both.
            if (!cat_layer) proxy <- proxy %>%
                  addControl(
                        circleLegend(
                              # In All mode the markers show observed
                              # richness, not the per-entity variable.
                              title  = HTML(paste(c("FIA plots: historical",
                                                    if (multi_mode()) paste("observed", comm_label())
                                                    else model$fia_var,
                                                    analog_vector_note()),
                                                  collapse = "<br/>")),
                              labels = legend_labels,
                              colors = legend_colors,
                              radii = legend_radii
                        ),
                        position = "topright",
                        layerId  = "point_legend",
                        className = "circle-legend leaflet-control"
                  )
      }

      if(input$vector == "current animation"){
            v_data <- format_velocity_data(connectivity()[["current_flow"]],
                                           connectivity()[["current_direction"]])
            proxy <- proxy %>%
                  addVelocity(
                        content = v_data,
                        options = velocityOptions(
                              colorScale = c("#000000", "#000000"),
                              lineWidth = 2,
                              velocityScale = 10^input$vector_scale * 25,
                              particleMultiplier = 0.005,
                              displayValues = FALSE,
                              velocityType = "Current",
                              emptyString = "No current data"
                        ),
                        group = "current_flow"
                  )
      }

      if(input$vector == "current field"){
            mag <- connectivity()[["current_flow"]]
            mag <- mag / minmax(mag)[2]
            spokes <- make_spokes(mag, connectivity()[["current_direction"]],
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
