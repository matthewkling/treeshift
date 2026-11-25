# Titles for map layers and the plots that show them. Every title is built from
# the same parts -- what it's about (the subject), what's measured and when
# (the quantity), and muted details such as units -- and drawn either as the
# map legend's stacked lines or as a one-line phrase for a plot axis:
#
#   Pinus edulis: change in basal area (ft² / acre)
#   All species: historical diversity (based on basal area; per ~0.07 ha plot)
#   historical analog support (equivalent number of exactly matching plots)
#
# The subject leads, set off by a colon, so the quantity reads as an ordinary
# noun phrase ("change in basal area") whatever the entity's name is.

# The quantity an FIA variable measures, as a noun phrase
var_noun <- function(fia_var) {
      switch(fia_var,
             "seedling presence"   = "seedling presence probability",
             "sapling presence"    = "sapling presence probability",
             "tree presence"       = "adult presence probability",
             "total basal area"    = "basal area",
             "number of trees"     = "adult tree density",
             "number of saplings"  = "sapling density",
             "number of seedlings" = "seedling density",
             "seedling fraction"   = "seedling regeneration fraction",
             "sapling fraction"    = "sapling regeneration fraction",
             "proportion basal area" = "proportion of basal area",
             fia_var)
}

# Units for an FIA variable measured per area, or NULL for the dimensionless
# ones (probabilities, fractions, counts of entities)
var_units <- function(fia_var) {
      switch(fia_var,
             "total basal area"    = "ft² / acre",
             "number of trees"     = "trees / acre",
             "number of saplings"  = "saplings / acre",
             "number of seedlings" = "seedlings / acre",
             NULL)
}

# The plot-scale qualifier, for values that depend on the area they're
# estimated over (support_dependent_vars, helpers/variables.R). Every
# community measure does: it describes the composition of a plot-sized
# footprint.
plot_scale_note <- "per ~0.07 ha plot"

# A noun phrase qualified by timeframe: "historical basal area", "change in
# basal area"
with_time <- function(noun, time) {
      switch(time,
             "delta" = paste("change in", noun),
             paste(time, noun))
}

# An FIA variable's axis label, with its units: "historical basal area (ft² /
# acre)", or with no timeframe, "basal area (ft² / acre)"
var_label <- function(fia_var, time = NULL) {
      noun <- var_noun(fia_var)
      if (!is.null(time)) noun <- with_time(noun, time)
      u <- var_units(fia_var)
      if (is.null(u)) noun else paste0(noun, " (", u, ")")
}

# What a community is made of, for measures that name its members
track_noun <- function(track, plural = FALSE) {
      switch(track,
             "Species"      = "species",
             "Forest type"  = if (plural) "forest types" else "forest type",
             "Forest group" = if (plural) "forest groups" else "forest group")
}

# A community measure's quantity and details, without timeframe. Richness is
# defined over presence alone, so it takes the life stage into its name;
# every other measure can be built from any variable, so it says which.
comm_noun <- function(measure, fia_var, track) {
      stage <- switch(fia_var,
                      "seedling presence" = "seedling ",
                      "sapling presence"  = "sapling ",
                      "tree presence"     = "adult ",
                      "")
      based_on <- paste("based on", var_noun(fia_var))
      occurrence <- fia_var %in% occurrence_vars
      member <- track_noun(track)

      switch(measure,
             "richness" = list(noun = paste0(stage, "richness"), notes = NULL),
             "dominant entity" = list(noun = paste("dominant", member), notes = based_on),
             "biggest winner"  = list(noun = paste(member, "with the largest gain"), notes = based_on),
             "biggest loser"   = list(noun = paste(member, "with the largest loss"), notes = based_on),
             # Sums of per-entity change: a count of members for a presence
             # variable, otherwise in the variable's units
             "expected gains"  = ,
             "expected losses" = {
                   dir <- sub("expected ", "", measure)
                   if (occurrence)
                         list(noun = paste("expected", member, dir), notes = based_on)
                   else list(noun = paste("expected", dir, "in", var_noun(fia_var)),
                             notes = var_units(fia_var))
             },
             list(noun = measure, notes = based_on))
}

# The parts of a layer's title. `stat` is the kind of layer (rast_stat());
# `focal` the entity's display name (focal_label()); `direction`
# input$analog_direction; `conn` the connectivity measure and `reach` its
# dispersal reach in km, as model$cs_reach holds it.
layer_title <- function(stat, time, track = NULL, focal = NULL, fia_var = NULL,
                        comm_measure = NULL, multi = FALSE, direction = NULL,
                        conn = "corridors", reach = NULL) {
      subject <- NULL
      notes <- NULL

      if (stat %in% c("support", "closeness")) {
            noun <- paste("analog", stat)
            # The "change" map of an evidence layer is not a change: it's the
            # harmonic mean of the two timeframes (ts_evidence_layer())
            quantity <- if (time == "delta") paste(noun, "across both timeframes")
            else with_time(noun, time)
            notes <- if (stat == "support") "equivalent number of exactly matching plots"
            else "1 = exact match"

      } else if (stat == "focal site similarity") {
            # Which climates are compared depends on the analog direction
            # (site_diff())
            quantity <- switch(direction,
                               "Reverse analogs" = "historical climate distance from the focal site's future climate",
                               "Forward analogs" = "future climate distance from the focal site's historical climate",
                               "historical climate distance from the focal site")
            notes <- "standardized units"

      } else if (stat == "connectivity" && !identical(conn, "corridors")) {
            subject <- focal
            quantity <- switch(conn,
                               "connectedness"      = paste("connected", time, "habitat"),
                               "tracking"           = "future habitat connected to historical habitat",
                               "colonization"       = "historical habitat connected to future habitat",
                               "reach"              = paste("proximity of suitable", time, "habitat"),
                               "reach_tracking"     = "proximity of suitable future habitat to historical habitat",
                               "reach_colonization" = "proximity of suitable historical habitat to future habitat")
            # These measures need a finite reach (cs_solve())
            within <- if (is.finite(as.numeric(reach %||% Inf))) paste("within", reach, "km")
                      else "within dispersal reach"
            notes <- paste0("habitat area ", within, ", ",
                            if (startsWith(conn, "reach")) "whatever lies between"
                            else "discounted by difficulty of access")

      } else if (stat == "connectivity") {
            subject <- focal
            quantity <- switch(time,
                               "historical" = "connectivity corridors within historical habitat",
                               "future"     = "connectivity corridors within future habitat",
                               "delta"      = "range shift corridors, historical to future habitat")
            notes <- "movement compared with an even landscape"

      } else if (multi) {
            subject <- focal
            cn <- comm_noun(comm_measure, fia_var, track)
            # Paired measures already compare the two timeframes
            quantity <- if (comm_measure %in% paired_measures) cn$noun else with_time(cn$noun, time)
            notes <- c(cn$notes, plot_scale_note)

      } else {
            subject <- focal
            noun <- var_noun(fia_var)
            quantity <- if (stat == "spread") with_time(paste("analog spread in", noun), time)
            else with_time(noun, time)
            notes <- c(var_units(fia_var),
                       if (fia_var %in% support_dependent_vars) plot_scale_note)
      }

      list(subject = subject, quantity = quantity, notes = notes)
}

# A title as one line, for a plot axis
title_phrase <- function(parts) {
      paste0(if (length(parts$subject)) paste0(parts$subject, ": "),
             parts$quantity,
             if (length(parts$notes)) paste0(" (", paste(parts$notes, collapse = "; "), ")"))
}
