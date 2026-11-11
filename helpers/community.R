# Community-composition mode: run limits, measures, support floors, and the
# per-entity reduction.

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

# Ceiling on a community-composition run, in gigabytes of per-entity surface
# held at once (both timeframes).
#
# The cost scales as cells x entities, and BOTH inputs are user-controlled: the
# grain slider goes to 1 km and the bounding box can span the data. At the
# default 10 km the stack is about 0.07 GB and lives in memory; at 1 km over the
# same box it is nearly 7 GB. On a shared server that is not a slow render, it
# is an out-of-memory failure that takes other sessions with it.
#
# Set generously enough to allow 2 km over a typical domain and to refuse 1 km
# over a wide one. Deployment-specific: raise it on a machine with headroom.
comm_max_gb <- 1.5

# Community measures that require entities able to CO-OCCUR at a plot.
#
# Forest types and groups are mutually exclusive plot classifications: each plot
# carries exactly one FORTYPCD, so multi_values() rows are one-hot. Richness is
# then 1 by construction, diversity is 1, dominance is 1 and evenness is 0 --
# every per-timeframe measure is degenerate.
#
# This is the plot-support constraint again rather than anything wrong with the
# classification: forest-type richness is meaningful over a cell, which contains
# many plots of potentially different types. It is degenerate only because the
# support is one plot.
#
# The paired measures will NOT be restricted this way. Bray-Curtis between two
# probability vectors over forest types measures how the type distribution
# shifts, which is well defined for a mutually exclusive classification.
cooccurrence_measures <- c("richness", "diversity", "evenness", "dominance")

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

# Label and step for the support floor, which are properties of the DENOMINATOR
# rather than of the layer being masked -- and the denominator changes with the
# variable, so a fixed label would be wrong most of the time.
#
# Single-entity ratios (the regeneration fractions) are conditional on the
# entity's occurrence probability, so the floor is on 0-1. Community measures
# are normalized by the summed per-entity response, whose units are whatever the
# variable measures.
support_floor_spec <- function(multi, fia_var, track) {
      # Defaults are nonzero because the unmasked layers are actively
      # misleading: which.max always names a winner, and a ratio is at its most
      # extreme where its denominator is smallest, so the empty ground is
      # exactly where these maps look most dramatic. The values below are
      # "barely anything there" in each unit -- deliberately low, since the
      # control exists to be tuned.
      if (!isTRUE(multi))
            return(list(label = "Minimum occurrence probability",
                        step = 0.05, default = 0.1))

      unit <- switch(fia_var,
                     "total basal area"      = list("Minimum ft\u00b2 BA / acre", 5, 5),
                     "number of trees"       = list("Minimum trees / acre", 5, 5),
                     "number of saplings"    = list("Minimum saplings / acre", 5, 5),
                     "number of seedlings"   = list("Minimum seedlings / acre", 10, 10),
                     "proportion basal area" = list("Minimum summed proportion", 0.1, 0.1),
                     # Occurrence variables: the denominator is expected richness.
                     # 0.5 rather than 1: a cell with 0.6 expected species is
                     # more likely than not to hold a tree, so cutting at 1 would
                     # discard genuinely occupied sparse woodland.
                     list(paste("Minimum expected",
                                if (identical(track, "Species")) "species" else "classes"),
                          0.5, 0.5))
      list(label = unit[[1]], step = unit[[2]], default = unit[[3]])
}

# Mask a ratio layer where the community or entity it describes is barely
# present. Ratios stay numerically well behaved as their denominator approaches
# zero but stop meaning anything: evenness, dominance, the regeneration
# fractions -- and, shortly, Bray-Curtis turnover and colonization demand --
# take their most extreme values exactly where nothing is there to be even,
# dominant, regenerating or turning over.
#
# combine_ratio() already masks the fractions at denom == 0; this is the same
# idea extended to denom ~ 0, which is the case that actually produces the
# misleading map.
#
# Masking rather than fading, deliberately: the extreme values from empty cells
# also drive seq_domain(), compressing the well-supported part of the map into a
# narrow band of the palette. Removing them fixes the legend as much as the map.
apply_support_floor <- function(x, support, floor) {
      if (!isTruthy(floor) || is.na(floor) || floor <= 0) return(x)
      terra::ifel(support >= floor, x, NA)
}

# Row-wise version of the community reductions in composition_layer(), for a
# plot x entity matrix. Used to paint the FIA markers with the OBSERVED
# counterpart of whatever the raster is modelling, on the same scale.
#
# Kept deliberately parallel to the raster path -- same formulas, different
# substrate -- since the two being comparable is what makes the markers a check
# on the surface rather than decoration.
#
# One interpretive wrinkle: for an occurrence variable the observed matrix is
# 0/1, so observed diversity collapses to observed richness and observed
# dominance to 1/richness. The MODELLED versions do not, because they act on
# continuous probabilities. The markers and the raster are both correct and are
# not measuring quite the same thing there.
comm_reduce_rows <- function(y, measure) {
      # Force the argument in this frame rather than inside as.matrix(): if a
      # caller passes a reactive whose req() fails, the resulting condition
      # escaping from S4 dispatch is reported as an evaluation error instead of
      # a silent cancellation.
      force(y)
      if (is.null(y) || !ncol(y)) return(NULL)
      v <- as.matrix(y)
      v[is.na(v)] <- 0
      s <- rowSums(v)

      # Categorical: which entity is largest at the plot. Returned as a column
      # index into `y`, the same index space the raster's which.max() uses, so
      # the two share a code vocabulary before lumping.
      if (measure %in% categorical_measures) {
            # Only the per-timeframe one has an observed counterpart: a plot is
            # surveyed at a single point in time, so there is no observed
            # "biggest winner" to paint the markers with.
            if (measure %in% paired_measures) return(NULL)
            out <- max.col(v, ties.method = "first")
            out[s == 0] <- NA_integer_
            return(out)
      }

      out <- if (identical(measure, "richness")) {
            s
      } else {
            # H = -sum(q log q) with q = v / s; see composition_layer().
            vlv <- rowSums(ifelse(v > 0, v * log(v), 0))
            H <- log(s) - vlv / s
            switch(measure,
                   "diversity" = exp(H),
                   "evenness"  = H / log(ncol(v)),
                   "dominance" = do.call(pmax, as.data.frame(v)) / s,
                   NULL)
      }
      if (is.null(out)) return(NULL)

      # A plot with nothing recorded has no composition to describe -- diversity,
      # evenness and dominance are all undefined when the community is empty.
      #
      # Richness is the exception: zero is a real, well-defined observation, and
      # it is exactly the value the model should be predicting low at. NA-ing
      # those plots dropped the entire bottom of the calibration curve and pulled
      # them out of the cross-validation pool, which is also why community
      # richness disagreed with the All-trees species-richness variable -- that
      # one counts zero correctly.
      if (!identical(measure, "richness")) out[s == 0] <- NA_real_
      out
}

# Which community measures are valid for a given track and variable. Invalid
# combinations are made unselectable rather than selectable-and-warned: these
# are permanent incompatibilities, not preconditions the user could satisfy by
# supplying something. (Contrast eval_scope, where an unavailable option means
# "no polygon selected yet" and saying so is more useful than hiding it.)
comm_measure_choices <- function(track, fia_var, time = "historical") {
      if (!isTruthy(fia_var) || !fia_var %in% multi_vars) return(character(0))

      static <- setdiff(comm_measures, paired_measures)
      # Mutually exclusive classifications make every per-timeframe measure
      # degenerate.
      if (track %in% c("Forest type", "Forest group"))
            static <- setdiff(static, cooccurrence_measures)
      # Richness counts entities, so it needs an occurrence variable: the same
      # sum over an abundance variable is community abundance.
      if (!fia_var %in% occurrence_vars) static <- setdiff(static, "richness")

      # Timeframe is the first cut, because "state or change" is a smaller and
      # more natural question than picking out of one mixed list. Under a single
      # timeframe only the static measures apply; under change they are joined by
      # the paired ones, and relabelled so the whole list reads as change
      # measures. Only the LABELS change -- the values stay the measure ids, so
      # nothing downstream has to know about this.
      if (!identical(time, "delta")) return(stats::setNames(static, static))

      # Categorical measures are dropped rather than relabelled: the difference
      # of two entity identities is not an entity identity. Their change
      # counterparts are separate measures (biggest winner, biggest loser).
      chg <- setdiff(static, categorical_measures)
      c(stats::setNames(chg, paste("change in", chg)),
        stats::setNames(paired_measures, paired_measures))
}

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
