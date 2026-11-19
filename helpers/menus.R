# Choice lists for the map-layer menus.

# Connectivity's map classes (server/map_layers.R), as multiples of an even
# landscape's flow. The detail is at the funneled end; everything under half
# is one class. Classes include their lower break.
cs_breaks <- c(0, 0.5, 1, 2, 5, 10, Inf)
cs_break_labels <- c("under ½ (impeded)", "½ to 1", "1 to 2", "2 to 5", "5 to 10",
                     "10 or more (funneled)")
cs_break_colors <- c("#bdbdbd", "#f0f0f0", "#fed976", "#fd8d3c", "#e31a1c", "#800026")

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
