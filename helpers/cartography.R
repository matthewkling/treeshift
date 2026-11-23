# Map styling: overlay colors, raster palettes, legends, and outline shapes.

# Colors offered for the focal polygon outline. A curated set rather than a full
# picker: the polygon is a boundary drawn over a colored raster, so what matters
# is contrast against a busy background, and a dozen strong, well-separated hues
# cover that better than an unconstrained wheel. The selectize `create` option
# still accepts anything typed, so a specific brand color is one keystroke away.
overlay_colors <- c("darkred", "#e05c00", "#c2b200", "#2b6a4d", "#1f78b4",
                    "#4c32a8", "#6a3d9a", "#e7298a", "black", "white", "gray40")

# Swatch renderer for the overlay color selectors. Defined once because both use
# it, and a divergence between them would be invisible until someone compared
# the two dropdowns side by side.
overlay_color_render <- function() {
      swatch <- paste0(
            "'<span style=\"display:inline-block;width:12px;height:12px;",
            "margin-right:8px;vertical-align:middle;border:1px solid #999;",
            "background:' + escape(item.value) + ';\"></span>'")
      I(paste0("{
            option: function(item, escape) {
                  return '<div style=\"padding:4px 8px;\">' + ", swatch, " + escape(item.label) + '</div>';
            },
            item: function(item, escape) {
                  return '<div>' + ", swatch, " + escape(item.label) + '</div>';
            }
      }"))
}

# Blend a color toward white. Used for the buffer ring, which follows the focal
# polygon's color but needs to stay distinguishable from it -- a tint plus the
# dashed line reads as "the same thing, extended" rather than as a second
# unrelated boundary.
lighten_color <- function(col, amount = 0.45) {
      v <- grDevices::col2rgb(col)[, 1]
      grDevices::rgb(t(v + (255 - v) * amount), maxColorValue = 255)
}

# One selector for either overlay. `create = TRUE` accepts any R color name or
# hex typed straight in, so the curated list is a starting point rather than a
# limit.
overlay_color_input <- function(id, label, current) {
      selectizeInput(id, label,
                     choices = union(overlay_colors, current),
                     selected = current,
                     options = list(create = TRUE, render = overlay_color_render()))
}

# Palette choices for the raster legend. The two sets are disjoint by design:
# a sequential or diverging ramp applied to entity identities implies an
# ordering that does not exist, and a qualitative set applied to a continuous
# surface is unreadable. One control swaps between them.
#
# The categorical options come from grDevices::palette.colors() rather than
# hcl.colors(): the hcl qualitative palettes sample a 1D path through HCL space,
# so past about eight categories adjacent colors differ only slightly and the
# result is a discretized gradient. scales::hue_pal() has the same problem --
# a hue sweep at fixed chroma and luminance.
#
# Polychrome 36 and Alphabet are built by maximizing perceptual distance across
# hue, chroma AND lightness together, which is what actually survives n = 10-25.
# Tableau 10 and Set 3 look better at small n but run out; Okabe-Ito is the
# colorblind-safe choice at eight. All are base R since 4.0, no new dependency.
rast_palettes_continuous  <- c("default", "viridis", "magma", "plasma", "inferno",
                               "RdYlBu", "Spectral", "RdBu", "BrBG")
rast_palettes_categorical <- c("default", "Polychrome", "Alphabet", "Tableau",
                               "Set 3", "Okabe-Ito")

## The map's legend panel (server/map_layers.R): one box at top right, with a
## section for the raster and one for the FIA plot markers. Each section is a
## title over its key. Continuous scales run horizontally under the title, so
## the title can use the panel's full width; categories and fixed classes are
## listed down the panel. Styled in ui/head.R (.ts-legend).

# A number for a legend label, to 3 significant figures
legend_num <- function(x) {
      vapply(x, function(v) format(signif(v, 3), big.mark = ",", scientific = FALSE, trim = TRUE),
             character(1))
}

# A title from its lines, the first `n_main` in bold and the rest muted
legend_title <- function(lines, n_main) {
      paste0("<div class='ts-legend-title'>",
             paste0("<div", ifelse(seq_along(lines) > n_main, " class='ts-legend-muted'", ""), ">",
                    htmltools::htmlEscape(lines), "</div>", collapse = ""),
             "</div>")
}

legend_section <- function(title, key) {
      paste0("<div class='ts-legend-section'>", title, key, "</div>")
}

# Tick positions on a log10 scale: whole decades, at most five, counted down
# from the highest so the top of the range is always labeled. NULL when the
# range doesn't span a decade.
log_ticks <- function(dom) {
      lo <- ceiling(min(dom))
      hi <- floor(max(dom))
      if (hi - lo < 1) return(NULL)
      ints <- lo:hi
      rev(seq(max(ints), min(ints), by = -ceiling((length(ints) - 1) / 4)))
}

# A horizontal color bar for a continuous palette over `dom`, labeled in the
# variable's units (`untransform` maps a transformed value back to them), at
# `cuts` on the drawn scale, or at round numbers if NULL
legend_bar <- function(pal, dom, untransform = NULL, opacity = 1, cuts = NULL) {
      stops <- seq(dom[1], dom[2], length.out = 21)
      if (is.null(cuts)) cuts <- pretty(dom, n = 4)
      cuts <- cuts[cuts >= min(dom) & cuts <= max(dom)]
      pos <- if (diff(range(dom)) > 0) (cuts - dom[1]) / diff(dom) * 100 else rep(0, length(cuts))
      labels <- legend_num(if (is.null(untransform)) cuts else untransform(cuts))
      paste0("<div class='ts-legend-bar' style='opacity:", opacity,
             "; background:linear-gradient(to right, ", paste(pal(stops), collapse = ", "), ");'></div>",
             "<div class='ts-legend-ticks'>",
             paste0("<span style='left:", round(pos, 1), "%;'>", labels, "</span>", collapse = ""),
             "</div>")
}

# Color swatches down the panel, one row per category or class
legend_swatches <- function(colors, labels, opacity = 1) {
      paste0("<div class='ts-legend-row'><i style='background:", colors, "; opacity:", opacity, ";'></i>",
             htmltools::htmlEscape(labels), "</div>", collapse = "")
}

# The FIA plot markers' key: outlined circles in a row, values beneath, and
# the non-sampled plots' small dots on a line of their own
legend_circles <- function(colors, labels) {
      paste0("<div class='ts-legend-circles'>",
             paste0("<div><svg width='12' height='12'><circle cx='6' cy='6' r='5' fill='", colors,
                    "' stroke='black' stroke-width='1'></circle></svg><span>", labels, "</span></div>",
                    collapse = ""),
             "</div>",
             "<div class='ts-legend-row'><svg width='12' height='12'><circle cx='6' cy='6' r='2.5' fill='gray'></circle></svg>Non-sampled plot</div>")
}

geo_circle <- function(pts, radius, thresh = .1) {
      radius <- radius * 1000 # m to km
      thresh <- thresh * 1000

      angles <- seq(from = 0, to = 360, by = 2)
      vertices <- geosphere::destPoint(p = pts, b = angles, d = radius)

      if(pts[1,1] > thresh) vertices[,1] <- ifelse(vertices[,1] > 0, vertices[,1], vertices[,1] + 360)
      if(pts[1,1] < -thresh) vertices[,1] <- ifelse(vertices[,1] < 0, vertices[,1], vertices[,1] - 360)

      vertices
}

clim_circle <- function(pt, radius, means, sds) {
      angles <- seq(from = 0, to = 360, by = 2) / 180 * pi
      as.data.frame(cbind(v1 = pt[1] * sds[1] + means[1] + cos(angles) * radius * sds[1],
                          v2 = pt[2] * sds[2] + means[2] + sin(angles) * radius * sds[2]))
}
