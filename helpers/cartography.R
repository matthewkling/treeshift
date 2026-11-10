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

circleLegend <- function(title, labels, colors, radii = 6) {
      # If only one radius is provided, repeat it for all labels
      if(length(radii) == 1) radii <- rep(radii, length(labels))

      items <- Map(function(lbl, col, rad) {
            # Use the 'rad' variable to define SVG dimensions and circle size
            sprintf(
                  "<div style='display:flex; align-items:center; margin:4px 0;'>
         <svg width='%d' height='%d' style='flex:0 0 auto;'>
           <circle cx='%d' cy='%d' r='%d' fill='%s'></circle>
         </svg>
         <span style='margin-left:8px;'>%s</span>
       </div>",
                  12, 12, # Keep container size consistent for alignment
                  6, 6,   # Center the circle in the 12x12 box
                  rad, col, lbl
            )
      }, labels, colors, radii)

      HTML(sprintf(
            "<div class='leaflet-control'
          style='background: rgba(255,255,255,0.75);
                 padding: 10px 12px;
                 margin: 2px;
                 border-radius: 5px;
                 box-shadow: 0 0 15px rgba(0,0,0,0.2);
                 color: #444;
                 line-height: 18px;'>
       <div style='font-weight:600; font-size:14px; margin-bottom:8px;'>%s</div>
       %s
     </div>",
            title,
            paste(items, collapse = "")
      ))
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
