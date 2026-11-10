# Helpers and state shared by the other server files.


point_radius <- reactive({10^input$vector_scale * 1000})

green <- "darkgreen"

# Active focal identifier and human-readable label
## All-entity ("community composition") mode. Guarded on track as well as
## the control, because input$entity_scope keeps its value while hidden --
## switching to All trees with the scope left on "multi" would otherwise
## try to run a community reduction over a track that is already an
## aggregate.
multi_mode <- reactive({
      isTruthy(input$entity_scope) && input$entity_scope == "multi" &&
            isTruthy(input$track) && input$track != "All trees"
})

focal_id <- reactive({
      if (multi_mode())
            return(paste0("all_", tolower(gsub(" ", "_", input$track))))
      switch(input$track,
             "Species"      = input$species,
             "Forest type"  = input$forest_type,
             "Forest group" = input$forest_group,
             "All trees"    = "all_trees")
})

focal_label <- reactive({
      if (multi_mode()) return(switch(input$track,
                                      "Species"      = "All species",
                                      "Forest type"  = "All forest types",
                                      "Forest group" = "All forest groups"))
      switch(input$track,
             "Species"      = input$species,
             "Forest type"  = {
                   nm <- ft_name(as.integer(input$forest_type))
                   if (is.na(nm)) paste("Forest type", input$forest_type) else nm
             },
             "Forest group" = {
                   nm <- fg_name(as.integer(input$forest_group))
                   if (is.na(nm)) paste("Forest group", input$forest_group) else nm
             },
             "All trees"    = "All trees")
})

# TRUE if the active response is a 0/1 indicator (any presence variable:
# "presence probability", "seedling presence", and the forest type/group
# classification probabilities); FALSE if continuous (basal area / counts).
response_is_binary <- reactive({
      # A community reduction is never binary even though it is built from
      # presence columns: richness is a count.
      if (multi_mode()) return(FALSE)
      isTruthy(input$fia_var) && grepl("presence", input$fia_var)
})

# Sequential color domain c(0, max) for markers/points, robust to the case
# where the values are all NA or all 0 (e.g. a fraction variable for a
# species absent from the current view) -- an empty/degenerate range makes
# colorNumeric error, so fall back to c(0, 1).
seq_domain <- function(x) {
      mx <- suppressWarnings(max(x, na.rm = TRUE))
      c(0, if (is.finite(mx) && mx > 0) mx else 1)
}

# Lower/upper bounds for the currently-displayed raster. Used by the
# density plot to keep KDE mass from leaking past natural limits (e.g.,
# presence probability is bounded [0, 1]; total BA is non-negative).
# Returns c(-Inf, Inf) when no informative bound applies.
display_bounds <- reactive({
      if (input$stat == "ESS") {
            c(0, Inf)
      } else if (input$stat == "focal site similarity") {
            c(0, Inf)
      } else if (input$stat == "connectivity") {
            if (input$cs_stat == "current direction") {
                  c(0, 360)
            } else {
                  # other connectivity stats: non-negative magnitudes
                  c(0, Inf)
            }
      } else if (input$stat == "FIA variable") {
            if (input$time == "delta") {
                  c(-Inf, Inf)
            } else if (response_is_binary()) {
                  c(0, 1)
            } else if (grepl("fraction", input$fia_var)) {
                  c(0, 1)  # regeneration fractions are bounded [0, 1]
            } else {
                  c(0, Inf)
            }
      } else {
            c(-Inf, Inf)
      }
})

# TRUE when the FIA-plot markers (map and focal-site scatterplots) should
# adopt the raster's palette. The markers always show historical FIA
# response on a [0, max] sequential scale; they're commensurable with the
# raster only when the raster is also showing FIA response on a
# sequential scale -- i.e. stat == "FIA variable" and time is historical
# or future (delta flips the raster to a diverging zero-centered scale).
# In all other cases the markers keep their own canonical gray90 -> green
# ramp. See palette discussion: this preserves the raster/marker
# alignment that holds by default without ever giving the markers an
# inappropriate scale.
marker_follows_raster_palette <- reactive({
      # Not for categorical layers: the map markers take the raster's
      # colorFactor directly, and the focal-site scatterplots are
      # continuous, so handing them a qualitative palette name would
      # reach scale_fill_distiller() with a palette it does not know.
      input$stat == "FIA variable" && input$time != "delta" && !categorical_layer()
})

## The panel plots are always categorical, but rast_palette holds a
## continuous ramp whenever the map is showing a continuous measure -- so
## without this they would fall back to the default palette the moment the
## user switched the map off dominant entity. Remembering the last
## qualitative choice keeps the panel on the palette the user picked.
last_cat_palette <- reactiveVal("default")
observeEvent(input$rast_palette, {
      if (input$rast_palette %in% rast_palettes_categorical)
            last_cat_palette(input$rast_palette)
})

## Swap the palette choices to match the layer type. observeEvent so the
## handler's read of input$rast_palette is isolated and cannot re-trigger
## its own write. "default" is in both sets, so it is the natural landing
## point when a selection does not survive the swap.
observeEvent(categorical_layer(), {
      ch <- if (categorical_layer()) rast_palettes_categorical
      else rast_palettes_continuous
      cur <- if (isTruthy(input$rast_palette)) input$rast_palette else "default"
      updateSelectInput(session, "rast_palette", choices = ch,
                        selected = if (cur %in% ch) cur else "default")
}, ignoreInit = TRUE)

# ggplot fill-scale equivalent of the leaflet choose_palette() logic.
# Resolve one of rast_palettes_categorical to n colors. Falls back to the
# default hue ramp for any name it does not recognise, which is what makes
# the choice-set swap safe: a stale continuous name arriving here during a
# flush produces a usable palette rather than an error.
# Drop neutrals: near-black, near-white and low-saturation grays. Gray is
# reserved for the pooled "other" band, so a gray in the palette reads as
# "other"; near-white leaves an entity with no visible swatch at all; and
# black is heavy enough to dominate a map it only holds a slice of.
# Okabe-Ito leads with black and Polychrome 36's second entry is near-white,
# so both need this.
drop_neutral_colors <- function(cols) {
      h <- grDevices::rgb2hsv(grDevices::col2rgb(cols))
      # Saturation alone catches grays AND near-whites (near-white is
      # bright but unsaturated); the value floor catches near-black. An
      # upper value bound would be wrong -- a pure hue like magenta or
      # green sits at v = 1 and is vivid, not washed out.
      keep <- h["s", ] > 0.25 & h["v", ] > 0.35
      if (!any(keep)) return(cols)
      cols[keep]
}

# When more entities are shown than the palette holds, recycled colors are
# shaded rather than repeated exactly. Two entities the same color is
# unavoidable past the palette length, but a darker and a lighter version
# of one hue can still be told apart, where two identical swatches cannot.
recycle_shaded <- function(cols, n) {
      if (n <= length(cols)) return(cols[seq_len(n)])
      k <- length(cols)
      vapply(seq_len(n), function(i) {
            base <- cols[((i - 1) %% k) + 1]
            cycle <- (i - 1) %/% k
            if (cycle == 0) return(base)
            h <- grDevices::rgb2hsv(grDevices::col2rgb(base))
            # Alternate darker / lighter on successive passes.
            darker <- (cycle %% 2) == 1
            v <- if (darker) max(0.20, h["v", 1] * 0.62) else min(0.98, h["v", 1] * 1.30)
            sat <- if (darker) h["s", 1] else h["s", 1] * 0.55
            grDevices::hsv(h["h", 1], min(1, sat), v)
      }, character(1))
}

cat_palette_colors <- function(n, pal) {
      if (!pal %in% rast_palettes_categorical) pal <- "default"
      if (n < 1) return(character(0))

      pick <- function(name, max_n) {
            cols <- drop_neutral_colors(
                  unname(grDevices::palette.colors(max_n, palette = name)))
            recycle_shaded(cols, n)
      }

      switch(pal,
             "Polychrome" = pick("Polychrome 36", 36),
             "Alphabet"   = pick("Alphabet", 26),
             "Tableau"    = pick("Tableau 10", 10),
             "Set 3"      = pick("Set 3", 12),
             "Okabe-Ito"  = pick("Okabe-Ito", 8),
             # Default adapts: Tableau 10 is the nicer set at small n, but
             # it runs out, and past that perceptual separation matters
             # more than individual colors looking good.
             if (n <= 9) pick("Tableau 10", 10) else pick("Polychrome 36", 36))
}

# rast_palette values are heterogeneous: "default" -> the canonical
# gray90 -> green gradient; "viridis"/"magma"/"plasma"/"inferno" ->
# viridis options; the rest ("RdYlBu"/"Spectral"/"RdBu"/"BrBG") ->
# RColorBrewer palettes via scale_fill_distiller. `limits` is passed
# through so the scale matches the marker domain; `guide` defaults to
# "none" since the focal-site plots don't show a fill legend.
marker_fill_scale <- function(limits, guide = "none") {
      p <- input$rast_palette
      if (p == "default" || !marker_follows_raster_palette()) {
            ggplot2::scale_fill_gradientn(colors = c("gray90", green),
                                          limits = limits, guide = guide)
      } else if (p %in% c("viridis", "magma", "plasma", "inferno")) {
            ggplot2::scale_fill_viridis_c(option = p,
                                          limits = limits, guide = guide)
      } else {
            # RColorBrewer diverging/sequential palettes
            ggplot2::scale_fill_distiller(palette = p,
                                          limits = limits, guide = guide)
      }
}


# Relative weight each analog plot carries in a prediction.
#
# analog_search() returns distances but not weights (pairs mode gives
# env_dist and geog_dist only), and the search in this panel is called with
# cutoff-only kernels -- so nothing upstream knows about weighting. It has
# to be recomputed here.
#
# *** THIS FORMULA MUST MATCH THE KERNELS IN aim_kernels():
#     kernel(max = theta * 3, weight = "gaussian", theta = theta)
# Assumed Gaussian with theta as the standard deviation, and the joint
# weight the product across the two dimensions. That reading comes from the
# app's own help text ("bandwidths specify the standard deviations") and
# from the cutoff sitting at three of them. If analogs uses a different
# constant, the point sizes here are subtly wrong while everything still
# looks reasonable -- so this is the line to check against the package. ***
analog_weight <- function(geog_dist, env_dist) {
      w <- exp(-0.5 * (geog_dist / input$theta_geog)^2) *
            exp(-0.5 * (env_dist / input$theta_clim)^2)
      # The cutoff is hard, not a tail: kernel(max = theta * 3) drops a
      # plot entirely rather than down-weighting it. Without this a plot
      # just past the cutoff gets exp(-4.5), which is small but not zero,
      # and is drawn as though it still counts for something.
      w[geog_dist > input$theta_geog * 3 | env_dist > input$theta_clim * 3] <- 0
      w
}

# Value column and fill scale for the analog plots, which color each plot
# by the response at that plot. For a categorical layer that response is an
# entity code, and handing codes to a continuous ramp paints arbitrary
# integers on a gradient -- the points come out gray-to-green and mean
# nothing. Convert to a labeled factor and use the entity colors instead,
# the same ones the map and the composition panel use.
analog_fill <- function(d) {
      if (categorical_layer()) {
            cc <- comm_categories()
            d$value <- factor(cc$labels[match(d$value, cc$codes)],
                              levels = cc$labels)
            # guide "none": the map legend beside these plots already keys
            # the same colors, and a second copy would crowd two small panels.
            list(d = d,
                 scale = ggplot2::scale_fill_manual(
                       values = stats::setNames(cc$colors, cc$labels),
                       drop = FALSE, guide = "none", na.value = "gray85"))
      } else {
            dom <- if (!response_is_binary()) seq_domain(d$value) else 0:1
            list(d = d, scale = marker_fill_scale(limits = dom))
      }
}
