# Map layers from a run's surfaces --------------------------------------------------
#
# One run serves every timeframe, measure and floor, so these take the surfaces
# and the display choices rather than settings. `time` is "historical",
# "future" or "delta" (future minus historical).

# Pick a timeframe from per-timeframe layers, or difference them
by_time <- function(time, bsl, fut) {
      switch(time,
             "historical" = bsl,
             "future"     = fut,
             "delta"      = fut - bsl,
             stop("Unknown timeframe '", time, "'"))
}

#' The single-entity variable's layer
#'
#' @param surfaces `ts_surfaces()`'s list.
#' @param support_floor For the regeneration fractions, the occurrence
#'   probability below which cells are masked (0 masks nothing). Applied per
#'   timeframe, so a change is blank wherever either end is unsupported.
ts_value_layer <- function(surfaces, fia_var, time = "historical", support_floor = 0) {
      # For the pooled-ratio variables the surfaces carry the occurrence layer
      # alongside the response; it is the denominator the fractions are
      # conditional on.
      is_frac <- grepl("fraction", fia_var)
      floor_at <- function(tf) {
            x <- tf$response
            if (is_frac && !is.null(tf$occ)) x <- apply_support_floor(x, tf$occ, support_floor)
            x
      }
      by_time(time, floor_at(surfaces$bsl), floor_at(surfaces$fut))
}

#' An evidence layer: what the plots behind each cell's estimate amount to
#'
#' @param layer "sum_weights" (support), "effective_weight" (closeness) or
#'   "spread".
#' @param time For change, the harmonic mean of the two timeframes, which
#'   leans toward the weaker: a change is only as well supported as its less
#'   supported end.
ts_evidence_layer <- function(surfaces, layer, time = "historical") {
      switch(time,
             "delta" = 2 / (1 / surfaces$fut[[layer]] + 1 / surfaces$bsl[[layer]]),
             "historical" = surfaces$bsl[[layer]],
             "future" = surfaces$fut[[layer]],
             stop("Unknown timeframe '", time, "'"))
}

#' A community measure's layer, from a community run's surfaces
#'
#' For the quantitative measures; the categorical ones (`categorical_measures`)
#' come from `ts_community_categories()`. Per-timeframe measures follow `time`;
#' the paired ones (`paired_measures`) compare the two timeframes and ignore it.
#'
#' @param surfaces `ts_multi_surfaces()`'s list.
#' @param support_floor Mask ratio measures where the community they describe
#'   sums to less than this, in the variable's units (0 masks nothing).
ts_community_layer <- function(surfaces, measure, time = "historical", support_floor = 0) {
      if (measure %in% categorical_measures)
            stop("'", measure, "' is categorical; use ts_community_categories()")
      if (!measure %in% comm_measures) stop("Unknown community measure '", measure, "'")

      reduce <- function(r) {
            ent <- setdiff(names(r), "sum_weights")
            v <- r[[ent]]

            # Everything below is layer-wise terra arithmetic rather than an
            # app() callback over cells, which matters at 143 layers.
            s <- terra::app(v, fun = "sum", na.rm = TRUE)

            # Shannon entropy of the composition, normalized to proportions
            # q = v / sum(v):  H = -sum(q log q) = log(s) - sum(v log v) / s.
            # ifel() guards log(0): a zero entity contributes nothing, and
            # 0 * -Inf would otherwise be NaN.
            shannon <- function() {
                  vlv <- terra::app(terra::ifel(v > 0, v * log(v), 0),
                                    fun = "sum", na.rm = TRUE)
                  log(s) - vlv / s
            }

            out <- switch(measure,
                          # Expected richness. Summing occurrence probabilities
                          # is deliberate: thresholding at p > 0.5 would discard
                          # the difference between 0.51 and 0.99, and Sum(p) is
                          # unbiased if p is calibrated.
                          "richness" = s,

                          # Hill N1 -- the effective number of entities, in the
                          # same units as richness, so the two are comparable.
                          "diversity" = exp(shannon()),

                          # Pielou's J, against the FULL entity set rather than
                          # a count present: with continuous probabilities
                          # almost every entity is nonzero somewhere. This makes
                          # evenness depend on the size of the entity set, which
                          # moves with theta_geog and the bbox.
                          "evenness" = shannon() / log(terra::nlyr(v)),

                          # Share of the community held by its largest member
                          "dominance" = terra::app(v, fun = "max", na.rm = TRUE) / s)

            # app(na.rm = TRUE) returns 0 where every layer is NA, i.e. outside
            # analog range -- where a single-entity surface is blank. Restore
            # nodata there. Every entity layer shares one NA mask (one call,
            # one set of weights, no NAs in the responses), so masking by any
            # one is enough.
            out <- terra::mask(out, v[[1]])

            # Support floor, per timeframe so a change is blank wherever
            # EITHER timeframe is unsupported. `s` is the natural denominator:
            # expected richness for an occurrence variable, total basal area or
            # stems for an abundance one. Richness is exempt: it IS the
            # denominator.
            if (!identical(measure, "richness"))
                  out <- apply_support_floor(out, s, support_floor)

            names(out) <- measure
            out
      }

      # Paired reduction: both stacks at once, entity by entity. The
      # per-entity differences cannot be recovered from the reduced surfaces
      # -- Sum(max(f-h, 0)) is not a function of Sum(f) and Sum(h). `a` is the
      # reference composition and `b` the comparison; for temporal pairing,
      # historical and future.
      reduce_paired <- function(a, b) {
            ent <- setdiff(names(a), "sum_weights")
            va <- a[[ent]]
            vb <- b[[ent]]

            d      <- vb - va
            gains  <- terra::app(terra::ifel(d > 0,  d, 0), fun = "sum", na.rm = TRUE)
            losses <- terra::app(terra::ifel(d < 0, -d, 0), fun = "sum", na.rm = TRUE)
            sa     <- terra::app(va, fun = "sum", na.rm = TRUE)
            sb     <- terra::app(vb, fun = "sum", na.rm = TRUE)

            out <- switch(measure,
                          # Bray-Curtis dissimilarity: sum|a-b| / sum(a+b).
                          # 0 is identical composition, 1 is no overlap.
                          "turnover" = (gains + losses) / (sa + sb),

                          # Absolute expected arrivals and departures, in the
                          # variable's own units
                          "expected gains"  = gains,
                          "expected losses" = losses,

                          # Gains as a share of the future community: the
                          # fraction of what will be there that is not there now
                          "colonization demand" = gains / sb,

                          # Losses as a share of the CURRENT community: how much
                          # of what is here now is expected to depart
                          "loss fraction" = losses / sa)

            out <- terra::mask(out, va[[1]])

            # Ratios only, each floored against its own denominator. The
            # absolute sums stay meaningful at low support.
            if (identical(measure, "turnover"))
                  out <- apply_support_floor(out, (sa + sb) / 2, support_floor)
            if (identical(measure, "colonization demand"))
                  out <- apply_support_floor(out, sb, support_floor)
            if (identical(measure, "loss fraction"))
                  out <- apply_support_floor(out, sa, support_floor)

            names(out) <- measure
            out
      }

      if (measure %in% paired_measures) return(reduce_paired(surfaces$bsl, surfaces$fut))
      by_time(time, reduce(surfaces$bsl), reduce(surfaces$fut))
}

#' A categorical community measure: which entity, per cell
#'
#' "dominant entity" is the largest entity in a timeframe ("historical" or
#' "future"); "biggest winner" and "biggest loser" the largest per-entity
#' gain and loss between them.
#'
#' @param top_n Entities kept by name, ranked by the number of cells they
#'   take; the rest are lumped as "other" (code 0).
#' @return list(raster, codes, labels): the raster's values are `codes`, which
#'   are column indices into the entity set (0 = other), and `labels` names
#'   them.
ts_community_categories <- function(surfaces, measure, time = "historical",
                                    support_floor = 0, top_n = 10) {
      if (!measure %in% categorical_measures)
            stop("'", measure, "' is not a categorical measure")
      ent <- setdiff(names(surfaces$bsl), "sum_weights")

      if (measure %in% paired_measures) {
            # Argmax and argmin of the per-entity DIFFERENCE, which cannot be
            # recovered from the reduced surfaces.
            va <- surfaces$bsl[[ent]]
            vb <- surfaces$fut[[ent]]
            dv <- vb - va
            idx <- if (identical(measure, "biggest winner")) terra::which.max(dv)
            else terra::which.min(dv)

            # Support over both eras together, so a cell that is empty now but
            # gains entities still qualifies. Requiring historical support
            # would mask exactly the arrivals this layer is for.
            support <- (terra::app(va, fun = "sum", na.rm = TRUE) +
                              terra::app(vb, fun = "sum", na.rm = TRUE)) / 2
      } else {
            r <- if (identical(time, "future")) surfaces$fut else surfaces$bsl
            v <- r[[ent]]
            idx <- terra::which.max(v)
            support <- terra::app(v, fun = "sum", na.rm = TRUE)
      }

      # Support floor FIRST. which.max always names a winner, including where
      # essentially nothing is expected; masking before the frequency count
      # also keeps those cells out of the top-N ranking.
      idx <- apply_support_floor(idx, support, support_floor)

      # Lump all but the most frequent N. Ranking by cell count answers "what
      # covers this map", which is what a discrete palette can show.
      f <- terra::freq(idx)
      f <- f[order(-f$count), , drop = FALSE]
      top <- utils::head(f$value, max(1, top_n))
      other <- setdiff(f$value, top)
      if (length(other)) idx <- terra::subst(idx, other, 0)
      names(idx) <- measure

      list(raster = idx,
           codes  = c(if (length(other)) 0, top),
           labels = c(if (length(other)) "other", ent[top]))
}
