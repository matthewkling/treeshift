# Choice lists for the map-layer menus.

cs_stat_names <- c("current flow", "current on target", "current loss", "current pinch points",
                   "current direction", "voltage", "conductance", "loss")
cs_stats <- c("current_flow", "current_destinations", "current_loss", "current_power",
              "current_direction", "voltage", "conductance", "loss")

# Left-sidebar menu choices, in connectivity-off and connectivity-on variants.
# Connectivity is a niche module, so its entries are kept out of the core
# `stat` and `vector` menus unless the user flips "Enable connectivity" in the
# right-sidebar Connectivity panel. The `*_base` sets are the UI defaults; the
# observeEvent(input$cs_enable) handler in the server swaps between the two.
# Listed explicitly rather than appended so "none" stays last in the vector menu.
# The evidence layers describe the FIA plots behind each cell's estimate
# (server/map_layers.R).
evidence_choices <- c("analog support" = "support", "analog closeness" = "closeness",
                      "analog spread" = "spread")
stat_choices_base <- c("FIA variable", evidence_choices, "focal site similarity")
stat_choices_cs   <- c("FIA variable", evidence_choices, "focal site similarity",
                       "connectivity")

vector_choices_base <- c("occupied FIA plots", "all FIA plots", "focal site analogs",
                         "none")
vector_choices_cs   <- c("occupied FIA plots", "all FIA plots", "focal site analogs",
                         "current animation", "current field", "none")

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
