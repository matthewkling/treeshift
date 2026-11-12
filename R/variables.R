# Variable catalog: climate variables, how a (track, metric, stage) choice
# resolves to an FIA variable, and the community measures.

# The climate variables the rasters carry
clim_vars <- function() c("tmean", "ppt", "aet", "cwd", "tmincm", "tmaxwm")

# Canonical fia_var for a (track, metric, stage) combination. Returns NULL for
# combinations metric_choices() would not have offered -- "abundance" has no
# "any" stage, for instance -- so callers must guard rather than assume a value.
resolve_var <- function(track, metric, stage) {
      if (track %in% c("Forest type", "Forest group")) return("presence probability")
      switch(metric,
             "occurrence" = switch(stage,
                                   "seedling" = "seedling presence",
                                   "sapling"  = "sapling presence",
                                   "adult"    = "tree presence",
                                   "any"      = "presence probability"),
             "abundance" = switch(stage,
                                  "seedling" = "number of seedlings",
                                  "sapling"  = "number of saplings",
                                  "adult"    = "number of trees"),
             "basal area"            = "total basal area",
             "species richness" = switch(stage,
                                         "seedling" = "seedling species richness",
                                         "sapling"  = "sapling species richness",
                                         "adult"    = "adult species richness",
                                         "any"      = "species richness"),
             "proportion basal area" = "proportion basal area",
             "regeneration fraction" = switch(stage,
                                              "seedling" = "seedling fraction",
                                              "sapling"  = "sapling fraction"))
}

# Occurrence-type variables. Richness is a sum of occurrence probabilities and
# is only defined over these: summing "total basal area" across entities gives
# community basal area, which is a perfectly good surface but is not richness.
occurrence_vars <- c("presence probability", "seedling presence",
                     "sapling presence", "tree presence")

# FIA variables that generalize to a per-entity column in the multi-entity
# response matrix. The ratio variables (seedling fraction, sapling fraction)
# carry the occurrence-gated numerator/denominator structure that
# ratio_pieces() hands to aim() as three coupled columns; generalizing those
# multiplies the column count by three and is deferred.
#
# Note "proportion basal area" IS here: its denominator is the plot total,
# shared across entities, so it is a plain per-entity column divided by a
# common vector rather than a pooled ratio.
multi_vars <- c("presence probability", "seedling presence", "sapling presence",
                "tree presence", "total basal area", "proportion basal area",
                "number of trees", "number of saplings", "number of seedlings")

# Paired measures: computed from the historical AND future stacks together, and
# undefined for a single timeframe. They ignore input$time rather than
# inheriting historical/future/delta from it.
#
# Note these are NOT deltas of a per-timeframe measure, and the distinction
# matters: net richness change of zero is entirely compatible with complete
# species replacement. Both families have to exist.
#
# The pairing rule is a separate axis from the measure. Every formula here takes
# two composition vectors and does not care whether they are separated in time
# (same cell, two timeframes), in space (focal site vs every cell, one
# timeframe), or both (the compositional counterpart of the reverse/forward
# analog machinery). Only temporal pairing is wired up so far, but the measures
# are written so the others slot in without redefining them.
paired_measures <- c("turnover", "expected gains", "expected losses",
                     "colonization demand", "loss fraction",
                     "biggest winner", "biggest loser")

# Measures whose value is an entity identity rather than a quantity. They need
# a discrete palette, nearest-neighbour reprojection, and no scale transform,
# and they are meaningless as a difference -- so they never take the "change in
# X" relabelling the continuous per-timeframe measures get.
#
# Note "dominant entity" is NOT restricted to co-occurring entities: over a
# mutually exclusive classification the argmax is the most likely forest type,
# which is a perfectly good map.
# "biggest winner" and "biggest loser" are BOTH categorical and paired: their
# value is an entity identity, and they are defined only across two timeframes.
# Membership in paired_measures puts them under "change"; membership here gives
# them the discrete palette and nearest-neighbour reprojection.
categorical_measures <- c("dominant entity", "biggest winner", "biggest loser")

# Community-level reductions across the per-entity stack. Each is a function of
# the whole entity set rather than of any one entity, which is why they occupy
# the entity-picker slot rather than the `metric` dropdown: in All mode, "how
# are entities combined" IS the selection the picker makes.
#
# Two families, distinguished by whether they are computed per timeframe (and so
# inherit historical / future / delta from input$time) or are inherently paired
# and defined only as a change. Only the per-timeframe family is implemented so
# far; the paired measures (turnover, expected gains and losses) are next.
#
# Richness is the one linear reduction here, and linearity means it commutes
# with the interpolation -- it is also available far more cheaply as an
# All-trees variable. It is kept here because it belongs to this family
# conceptually, and because the two routes agreeing is a free check on the whole
# multi-entity path. The other three genuinely require per-entity surfaces:
# no plot-level aggregate can be interpolated to recover them.
comm_measures <- c("richness", "diversity", "evenness", "dominance",
                   "dominant entity", paired_measures)

# Mask a ratio layer where the community or entity it describes is barely
# present. Ratios stay numerically well behaved as their denominator approaches
# zero but stop meaning anything: evenness, dominance, the regeneration
# fractions -- and, shortly, Bray-Curtis turnover and colonization demand --
# take their most extreme values exactly where nothing is there to be even,
# dominant, regenerating or turning over.
#
# ts_surfaces() already masks the fractions at denom == 0; this is the same
# idea extended to denom ~ 0, which is the case that actually produces the
# misleading map.
#
# Masking rather than fading, deliberately: the extreme values from empty cells
# also drive seq_domain(), compressing the well-supported part of the map into a
# narrow band of the palette. Removing them fixes the legend as much as the map.
apply_support_floor <- function(x, support, floor) {
      if (!length(floor) || is.na(floor) || floor <= 0) return(x)
      terra::ifel(support >= floor, x, NA)
}
