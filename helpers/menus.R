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
stat_choices_base <- c("FIA variable", "effective sample size" = "ESS",
                       "focal site similarity")
stat_choices_cs   <- c("FIA variable", "effective sample size" = "ESS",
                       "focal site similarity", "connectivity")

vector_choices_base <- c("occupied FIA plots", "all FIA plots", "focal site analogs",
                         "none")
vector_choices_cs   <- c("occupied FIA plots", "all FIA plots", "focal site analogs",
                         "current animation", "current field", "none")
