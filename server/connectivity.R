## CONNECTIVITY -------------------------------------------------
#
# Circuit-theory connectivity of the mapped entity's habitat, by circuitscaper's
# random-dipole method. Two measures (conn_choices(), helpers/menus.R), each
# with three maps picked by the timeframe.
#
# Corridors: where movement is funneled.
#   historical, future  connectivity within that period's habitat (dp_omni())
#   change              range shift pathways linking historical habitat to
#                       future habitat (dp_directionality(): the undirected map,
#                       plus the net flow the current vector layers draw)
# Each map shows normalized current: the flow through a cell relative to the
# flow it would carry on an even landscape with the same habitat and reach.
# Well above 1, movement is funneled (a corridor or pinch point); near 1 it
# spreads as on an even landscape; well below 1 it is impeded.
#
# Connected habitat and proximity of suitable habitat: how much habitat each
# habitat cell can reach (dp_accessibility()), as area (km2).
#   historical, future  the rest of that period's habitat within reach
#   change              from historical habitat to future habitat (tracking),
#                       and from future habitat to historical (colonization)
# Proximity is the habitat area within the dispersal reach, whatever
# lies between. Connected habitat counts it in full where it is as easy to
# reach as through continuous, fully suitable habitat on unmodified land, and
# less where it is harder. Both grow with the reach. One solve gives both, so
# switching between them doesn't solve again. Shown on habitat cells only.
#
# Each map is its own reactive, so only the map being shown is solved, and
# nothing is solved while connectivity isn't on the map, however the model
# settings change in the meantime.

# Conductance of a cell with no habitat, relative to the mean
cs_floor <- 1e-3

# Habitat below this share of the maximum is not a source
cs_threshold <- 0.01

# Radius of the site each cell stands for in the habitat measures, in km:
# fixed on the ground, so the measures don't change with the grain size.
# dp_accessibility() takes it in pixels (cs_compute()). Where the grain is
# coarser than this, the site is the cell alone.
cs_site_radius <- 5

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
# Cells outside the circuit stay NA (the mask: terra takes NA^0 as 1). Scaled
# by its mean, so its "ideal" attribute is the conductance of full habitat on
# unmodified land, which connected habitat is measured against.
cs_conductance <- function(habitat, permeability = NULL, landuse = 0) {
      w <- landuse / 100
      if (!is.null(permeability)) habitat <- mask(habitat^(1 - w) * permeability^w, habitat)
      m <- global(habitat, "mean", na.rm = TRUE)$mean
      r <- habitat / m + cs_floor
      attr(r, "ideal") <- 1 / m + cs_floor
      r
}

# One map, from the result cache when these settings have been solved before
# (R/cache.R). `era` is "hist", "fut" or "shift" for corridors, and the same
# with "_iso" for the habitat measures.
cs_solve <- function(era) {
      # Belt and braces: the connectivity layers are withdrawn when it
      # doesn't apply, but a stale selection during a flush must not reach
      # a solve
      req(cs_available())
      # The habitat measures count the habitat within reach, which an
      # unlimited reach doesn't define
      if (endsWith(era, "_iso"))
            shiny::validate(need(is.finite(as.numeric(model$cs_reach)),
                                 "This measure needs a dispersal reach: choose one other than 'no limit' under 'Landscape connectivity'."))
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
            if (endsWith(era, "_iso")) {
                  # The change map runs on the habitat of both periods, as
                  # the range shift pathways do
                  shift <- era == "shift_iso"
                  hab <- if (shift) sqrt(habitat$hist * habitat$fut) else habitat[[sub("_iso$", "", era)]]
                  conductance <- cs_conductance(hab, permeability, landuse)
                  # The side of a cell in km, as the square root of its mean
                  # area: cells on the lon/lat grid aren't square
                  px_km <- sqrt(global(cellSize(hab, unit = "km"), "mean")$mean)
                  r <- circuitscaper::dp_accessibility(
                        conductance, reach,
                        source_strength = if (shift) habitat$hist else hab,
                        destination = if (shift) habitat$fut,
                        radius_units = "km", resistance_is = "conductances",
                        site_radius = cs_site_radius / px_km,
                        source_threshold = cs_threshold, seed = cs_seed,
                        normalize = FALSE,
                        reference_resistance = 1 / attr(conductance, "ideal"))
                  r[[grep("(accessible|target)_area$", names(r))]]
            } else if (era == "shift") {
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
cs_hist_iso <- reactive(cs_solve("hist_iso"))
cs_fut_iso <- reactive(cs_solve("fut_iso"))
cs_shift_iso <- reactive(cs_solve("shift_iso"))

# The map for the measure and timeframe shown, over the whole model domain
connectivity_layer <- reactive({
      # The measure's choices follow the timeframe one flush behind
      # (server/controls.R)
      req(input$conn_stat %in% conn_choices(input$time))
      if (input$conn_stat == "corridors") {
            r <- switch(input$time, historical = cs_hist(), future = cs_fut(), delta = cs_shift())
            return(r[["normalized_current"]])
      }
      # Under change, the direction says whose habitat the map is on, and the
      # layers carry its prefix
      m <- input$conn_stat
      from <- if (m %in% c("tracking", "reach_tracking")) "hist"
              else if (m %in% c("colonization", "reach_colonization")) "fut"
              else switch(input$time, historical = "hist", future = "fut")
      prefix <- if (input$time != "delta") "" else if (from == "hist") "origin_" else "destination_"
      layer <- if (startsWith(m, "reach")) "target_area" else "accessible_area"
      r <- switch(input$time, historical = cs_hist_iso(), future = cs_fut_iso(),
                  delta = cs_shift_iso())[[paste0(prefix, layer)]]
      # Shown for the habitat it describes: elsewhere it would say how much
      # habitat a place would reach if it were habitat, which is easy to misread
      habitat <- cs_habitat()[[from]]
      mask(r, habitat >= cs_threshold, maskvalues = FALSE)
})
