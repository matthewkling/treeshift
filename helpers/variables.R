# Climate and FIA variable catalog: available variables, their units, and how
# a (track, metric, stage) choice resolves to an FIA variable.


# Native (un-standardized) units for each climate variable, used to label the
# climate-model selector chart after deltas are back-converted from z-scores.
# The stored climate rasters are standardized; a delta back-converts to native
# units by multiplying the z-score delta by the variable's SD (the mean
# cancels in a difference). NOTE: ppt is log10-transformed before standardizing,
# so its delta is in log10-space, not mm -- labelled accordingly.
clim_units <- c(tmean  = "\u00b0C",
                tmincm = "\u00b0C",
                tmaxwm = "\u00b0C",
                aet    = "mm",
                cwd    = "mm",
                ppt    = "log10 mm")

# Variables the climate-relationships plot can put on its x axis. Deliberately
# wider than the set the model is fitted on: a relationship against a variable
# the model did not use is still a real property of the landscape, and often the
# more interesting one -- richness peaking at mid elevation says something even
# though elevation is not a predictor.
#
# Kept entirely separate from the modeling path. Widening clim_pool() instead
# would be tidier-looking and much riskier: pool()'s analog index is built from
# whatever climate columns it is handed, so a missed subset would silently build
# a six-dimensional index instead of two and change every prediction with
# nothing on screen to say so.
#
# Extra (non-climate) rasters are declared here with their file and units. Any
# single-layer raster on the same grid slots in; back_transform() already leaves
# variables with no clim_scales entry untouched, which is right for a raster
# stored in native units.
display_extra_vars <- list(
      elevation = list(file = "data/elevation.tif", units = "m")
)
display_extra_vars <- display_extra_vars[
      vapply(display_extra_vars, function(x) file.exists(x$file), logical(1))]

# Available x-axis variables, and which of them the model can be fitted on.
display_vars <- function() c(clim_vars(), names(display_extra_vars))

# Native-unit label for a climate variable. Climate rasters are stored in
# standardized (z-score) units, so any plot showing native units must first
# back-convert via clim_scales (z * sd + mean). Falls back to "standardized"
# for a variable with no unit entry.
unit_lab <- function(v) {
      if (!is.null(display_extra_vars[[v]])) return(display_extra_vars[[v]]$units)
      u <- unname(clim_units[v])
      if (is.null(u) || is.na(u)) "standardized" else u
}

## Variables available for a track AND life stage. Stage is now the
## upstream control (see the UI comment), so this is where the
## stage-specific gaps in the FIA data get expressed:
##
##   basal area              only under "any" -- FIA's basal_area is not
##                           broken out by size class
##   proportion basal area   same
##   abundance               needs a specific stage; there is no combined
##                           stem count across seedlings, saplings and trees
##   regeneration fraction   seedling or sapling only; it is a ratio
##                           AGAINST adults, so "adult" and "any" are
##                           meaningless
metric_choices <- function(track, stage = "any") {
      if (track %in% c("Forest type", "Forest group")) return("occurrence")
      any_stage <- identical(stage, "any")

      c("occurrence",
        if (!any_stage) "abundance",
        if (any_stage) "basal area",
        if (identical(track, "All trees")) "species richness",
        if (any_stage && identical(track, "Species")) "proportion basal area",
        if (stage %in% c("seedling", "sapling")) "regeneration fraction")
}

# There is no stage_choices() counterpart: since the dependency was inverted,
# stage is the upstream control and its options are the same four everywhere.
# The observer that used to repopulate them from `metric` is gone -- it was the
# old backwards direction, and it also overwrote the labeled choice vector
# ("all life stages") with bare values and reset the user's selection.

# Variables whose value depends on the SUPPORT it was estimated over, and so
# need the plot-scale qualifier in the legend.
#
# Densities are not here: a kernel-weighted mean of per-acre basal area is the
# expected per-acre basal area whatever area you evaluate it over, so support
# affects their variance but not their expectation.
#
# Occurrence is different in kind, because presence is a threshold on an
# underlying density -- P(>=1 stem in 0.07 ha) and P(>=1 stem in a 10 km cell)
# are different numbers with no unit conversion between them.
#
# Two non-obvious members. "proportion basal area" is a ratio of means, which
# would be scale-free on its own, but ts_surfaces() multiplies it by the
# occurrence layer, and that multiplier carries the support dependence. The
# regeneration fractions are ratios of occurrence-based counts, and the support
# does not cancel: both components saturate toward 1 as support grows, dragging
# the ratio toward 0.5.
support_dependent_vars <- c("presence probability", "seedling presence",
                            "sapling presence", "tree presence",
                            "proportion basal area",
                            "seedling fraction", "sapling fraction",
                            "species richness", "seedling species richness",
                            "sapling species richness", "adult species richness")
