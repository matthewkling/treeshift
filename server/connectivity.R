## CONNECTIVITY -------------------------------------------------
#
# Circuit-theory connectivity of the mapped entity's habitat, by circuitscaper's
# random-dipole method. The timeframe picks one of three maps:
#   historical, future  connectivity within that period's habitat (dp_omni())
#   change              range shift pathways linking historical habitat to
#                       future habitat (dp_directionality(): the undirected map,
#                       plus the net flow the current vector layers draw)
# Each is its own reactive, so only the map being shown is solved, and nothing
# is solved while connectivity isn't on the map, however the model settings
# change in the meantime.
#
# Each map shows normalized current: the flow through a cell relative to the
# flow it would carry on an even landscape with the same habitat and reach.
# Well above 1, movement is funneled (a corridor or pinch point); near 1 it
# spreads as on an even landscape; well below 1 it is impeded.

# Conductance of a cell with no habitat, relative to the mean
cs_floor <- 1e-3

# Habitat below this share of the maximum is not a source
cs_threshold <- 0.01

# Fixed, so a map is reproducible and cacheable
cs_seed <- 1L

# Julia threads for the solves: all the machine's cores, as analogs uses for
# the model. circuitscaper starts Julia with one unless asked; the results
# don't depend on the thread count.
cs_threads <- parallel::detectCores()

# The two periods' habitat, scaled by their common maximum: the periods weigh
# alike, and the habitat threshold is relative, whatever the variable's units.
# Cells outside the contiguous US or the climate data are not part of the
# circuit (NA), so the map stops at the border. Cells within it that the model
# leaves without an estimate (future climates with no analogs) are not
# habitat. Habitat is scaled by the cell's land fraction, so open water is not
# habitat.
cs_habitat <- reactive({
      bsl <- aim()$bsl$response
      fut <- aim()$fut$response
      top <- max(global(c(bsl, fut), "max", na.rm = TRUE)$max)
      shiny::validate(need(is.finite(top) && top > 0,
                           "Connectivity needs habitat, and the model finds none in this area."))
      outside <- is.na(aim()$bsl$sum_weights)
      land <- cs_ghm()$land
      scale <- function(r) {
            r <- r / top
            r[is.na(r)] <- 0
            r <- r * land
            r[outside] <- NA
            r
      }
      list(hist = scale(bsl), fut = scale(fut))
})

# gHM on the model grid (data_prep/ghm_prep.R): `land`, the share of the cell's
# 1 km cells that are land (gHM has no value on water, and the climate grid
# includes lakes and coastal water), and `ghm`, the mean gHM over that land.
cs_ghm <- reactive({
      template <- aim()$bsl$response
      g <- rast(file.path(core_data$dir, "ghm.tif")) |> crop(template, snap = "out")
      to_grid <- function(r) resample(r, template, method = "average")
      land <- to_grid(!is.na(g))
      land[is.na(land)] <- 0
      list(land = land, ghm = to_grid(g))
})

# How readily land use lets a cell be crossed, on the same scale of contrast as
# habitat: from 1 on unmodified land down to cs_floor (as hard to cross as
# non-habitat) on fully modified land, log-linear in gHM, and scaled by the
# land fraction so water is hard to cross. NULL when land use has no weight.
cs_permeability <- reactive({
      if (model$cs_landuse == 0) return(NULL)
      g <- cs_ghm()
      ghm <- g$ghm
      ghm[is.na(ghm)] <- 0
      g$land * cs_floor^ghm
})

# Conductance: a weighted geometric mean of habitat and permeability, with
# weight `landuse` (%) on permeability. At 0 it is habitat alone; at 100 it is
# land use alone, and habitat only places the sources. Never quite zero.
cs_conductance <- function(habitat, permeability = NULL, landuse = 0) {
      w <- landuse / 100
      if (!is.null(permeability)) habitat <- habitat^(1 - w) * permeability^w
      habitat / global(habitat, "mean", na.rm = TRUE)$mean + cs_floor
}

# One map, from the result cache when these settings have been solved before
# (R/cache.R). `era` is "hist", "fut" or "shift".
cs_solve <- function(era) {
      # Belt and braces: the connectivity layers are withdrawn when it
      # doesn't apply, but a stale selection during a flush must not reach
      # a solve
      req(cs_available())
      settings <- lapply(setNames(nm = c(surfaces_setting_ids, connectivity_setting_ids)),
                         function(id) model[[id]])
      logged_run("connectivity", ts_connectivity_key(settings, core_data, era),
                 function() cs_compute(era, cs_habitat(), as.numeric(model$cs_reach),
                                       cs_permeability(), model$cs_landuse),
                 describe = function(r) list(era = era, cells = terra::ncell(r)),
                 message = "Computing connectivity.",
                 detail = "Large or fine-grained maps can take a minute or more.")
}

cs_compute <- function(era, habitat, reach, permeability = NULL, landuse = 0) {
      # Starts Julia on first use; afterwards returns at once
      circuitscaper::cs_setup(threads = cs_threads)
      solve <- function() {
            if (era == "shift") {
                  r <- circuitscaper::dp_directionality(
                        cs_conductance(sqrt(habitat$hist * habitat$fut), permeability, landuse),
                        habitat$hist, habitat$fut,
                        radius = reach, radius_units = "km", resistance_is = "conductances",
                        source_threshold = cs_threshold, seed = cs_seed)
                  c(setNames(r$omni_normalized_current, "normalized_current"), r[[c("flow_x", "flow_y")]])
            } else {
                  r <- circuitscaper::dp_omni(
                        cs_conductance(habitat[[era]], permeability, landuse), reach,
                        source_strength = habitat[[era]],
                        radius_units = "km", resistance_is = "conductances",
                        source_threshold = cs_threshold, seed = cs_seed)
                  r[["normalized_current"]]
            }
      }
      tryCatch(solve(), error = function(e) {
            # dp_*() stop when no two habitat cells are within reach
            if (!grepl("within reach|no pairs", conditionMessage(e), ignore.case = TRUE)) stop(e)
            shiny::validate(need(FALSE, paste(
                  "No habitat is within the dispersal reach of other habitat here.",
                  "Try a longer reach under 'Landscape connectivity'.")))
      })
}

cs_hist <- reactive(cs_solve("hist"))
cs_fut <- reactive(cs_solve("fut"))
cs_shift <- reactive(cs_solve("shift"))

# The map for the timeframe shown, over the whole model domain
connectivity_layer <- reactive({
      r <- switch(input$time, historical = cs_hist(), future = cs_fut(), delta = cs_shift())
      r[["normalized_current"]]
})
