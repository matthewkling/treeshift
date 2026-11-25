# Choice lists for the map-layer menus.

# Connectivity's map classes (server/map_layers.R), as multiples of an even
# landscape's flow. The detail is at the funneled end; everything under half
# is one class. Classes include their lower break.
cs_breaks <- c(0, 0.5, 1, 2, 5, 10, Inf)
cs_break_labels <- c("under ½ (impeded)", "½ to 1", "1 to 2", "2 to 5", "5 to 10",
                     "10 or more (funneled)")
cs_break_colors <- c("#bdbdbd", "#f0f0f0", "#fed976", "#fd8d3c", "#e31a1c", "#800026")

# Connectivity's measures (server/connectivity.R): where movement through the
# habitat is funneled (corridors); how much habitat each habitat cell can
# reach, discounted by how hard it is to reach (connected habitat); and the
# habitat within its dispersal reach, whatever lies between (proximity of
# suitable habitat). Under change, the habitat measures have two directions: from each
# cell of historical habitat to future habitat (can it track its habitat?),
# and from each cell of future habitat to historical habitat (can it be
# colonized?).
conn_choices <- function(time) {
      if (identical(time, "delta"))
            c("corridors",
              "connected future habitat" = "tracking",
              "connected historical habitat" = "colonization",
              "proximity of suitable future habitat" = "reach_tracking",
              "proximity of suitable historical habitat" = "reach_colonization")
      else c("corridors", "connected habitat" = "connectedness",
             "proximity of suitable habitat" = "reach")
}

# The habitat measures' classes (server/map_layers.R) for values `x` in km2:
# on a log scale, as values span orders of magnitude, with breaks at round
# numbers (1-2-5, 1-3 or tenfold steps: the finest that gives at most
# `max_classes`), spanning the values on the map and no more. Values under
# 1 km2 share one class, labeled isolated with `anchor` (connected habitat,
# which accounts for what lies between). Classes include their lower break.
# Red at the low end, so isolated habitat stands out.
cs_area_bins <- function(x, anchor = FALSE, max_classes = 8) {
      x <- x[is.finite(x)]
      if (!length(x)) x <- 0
      breaks <- c(0, 1)
      if (max(x) >= 1) for (steps in list(c(1, 2, 5), c(1, 3), 1)) {
            series <- sort(as.vector(outer(steps, 10^(0:12))))
            lo <- max(series[series <= max(min(x), 1)])
            hi <- min(series[series > max(x)])
            breaks <- c(if (min(x) < 1) 0, series[series >= lo & series <= hi])
            if (length(breaks) - 1 <= max_classes) break
      }
      km2 <- function(v) paste0(scales::comma(v), " km\u00b2")
      a <- head(breaks, -1)
      b <- breaks[-1]
      labels <- ifelse(a == 0, paste("under", km2(b)), paste(scales::comma(a), "to", km2(b)))
      if (anchor && a[1] == 0) labels[1] <- paste(labels[1], "(isolated)")
      ramp <- grDevices::colorRampPalette(scales::brewer_pal(palette = "RdYlBu")(11))
      list(breaks = breaks, labels = labels, palette = ramp(length(a)))
}

# The weight (%) connectivity's conductance gives human land use over habitat
# quality (server/connectivity.R)
cs_landuse_steps <- seq(0, 100, by = 25)

# Connectivity's dispersal reach (km, as text so "Inf" survives a menu)
cs_reach_choices <- c("10 km" = "10", "25 km" = "25", "50 km" = "50", "100 km" = "100",
                      "250 km" = "250", "no limit" = "Inf")

# The variables connectivity is offered for: those that measure how much of
# the entity a place holds, so they can stand for habitat. Ratio variables
# (proportion basal area, regeneration fraction) and richness don't.
cs_metrics <- c("occurrence", "abundance", "basal area")

# Left-sidebar menu choices, without and with the connectivity entries, which
# are offered only where connectivity applies (url_cs_available()). The
# `*_base` sets are the UI defaults; an observer in server/controls.R swaps
# between the two. Listed explicitly rather than appended so "none" stays
# last in the vector menu.
# The evidence layers describe the FIA plots behind each cell's estimate
# (server/map_layers.R). They share one "analog evidence" entry in the raster
# menu, with the measure chosen in a second menu.
evidence_choices <- c("analog support" = "support", "analog closeness" = "closeness",
                      "analog spread" = "spread")
stat_choices_base <- c("FIA variable", "analog evidence" = "evidence", "focal site similarity")
stat_choices_cs   <- c(stat_choices_base, "connectivity")

vector_choices_base <- c("occupied FIA plots", "all FIA plots", "focal site analogs",
                         "none")
vector_choices_cs   <- c("occupied FIA plots", "all FIA plots", "focal site analogs",
                         "range shift flow (animated)" = "current animation",
                         "range shift flow (arrows)" = "current field", "none")

# Single-or-community choices. Only the labels follow the track ("single
# species" vs "single forest type"); the values stay single/multi. All trees
# has no such choice (its menu is hidden), so it gets the species labels.
entity_scope_choices <- function(track) {
      one <- switch(track,
                    "Forest type"  = "single forest type",
                    "Forest group" = "single forest group",
                    "single species")
      stats::setNames(c("single", "multi"), c(one, "community composition"))
}

# Entity and life-stage menus (labels shown, values used throughout the
# server). Shared by the UI and the settings record.
track_choices <- c("species"      = "Species",
                   "forest type"  = "Forest type",
                   "forest group" = "Forest group",
                   "all trees"    = "All trees")
stage_choices <- c("seedlings only"  = "seedling",
                   "saplings only"   = "sapling",
                   "adults only"     = "adult",
                   "all life stages" = "any")
