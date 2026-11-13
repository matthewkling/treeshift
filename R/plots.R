# FIA plots and their responses ------------------------------------------------------

#' The analog pool's bounding box
#'
#' The model bounding box padded by the kernel's reach (3 x theta_geog km), so
#' plots just outside it still serve as analogs for cells near its edge. Km to
#' degrees with a conservative ratio: 1 degree of latitude is ~111 km, of
#' longitude ~71 km at 50 N, so /100 over-pads by 10-30% at CONUS latitudes,
#' which is harmless.
ts_pool_bbox <- function(bbox_x, bbox_y, theta_geog) {
      pad_deg <- theta_geog * 3 / 100
      list(x = c(bbox_x[1] - pad_deg, bbox_x[2] + pad_deg),
           y = c(bbox_y[1] - pad_deg, bbox_y[2] + pad_deg))
}

# FIA rows within the pool bounding box
fia_in_box <- function(data, pool_bbox) {
      dplyr::filter(data$fia,
                    dplyr::between(x, pool_bbox$x[1], pool_bbox$x[2]),
                    dplyr::between(y, pool_bbox$y[1], pool_bbox$y[2]))
}

#' One row per FIA plot in the pool, summarized around a focal species
#'
#' The row order defines the plot index that the pool and every response are
#' aligned to: plots sorted by x, then y.
#'
#' Summarized with rowsum() and tabulate() over a plot index rather than a
#' grouped summarize(): with tens of thousands of plots, evaluating each summary
#' once per plot took several seconds, and this runs on every model run.
#'
#' @param species The focal species. Its columns (`pres`, `batot`, ...) are
#'   used by the Species track only; the plot totals serve the other tracks.
ts_plots <- function(data, pool_bbox, species) {
      # Plot-level type->group lookup (one row per forest_type code)
      ft_to_fg <- dplyr::distinct(data$forest_types, forest_type, forest_group)

      # Rows grouped by plot, keeping their order within each plot
      f <- fia_in_box(data, pool_bbox)
      f <- f[order(f$x, f$y), ]
      first <- c(TRUE, diff(f$x) != 0 | diff(f$y) != 0)
      plot <- cumsum(first)
      n <- sum(first)

      focal <- f$species == species & !is.na(f$species)
      sums <- function(v) unname(rowsum(v, plot, reorder = FALSE)[, 1])
      any_of <- function(cond) as.integer(sums(as.integer(cond)) > 0)
      # Distinct species per plot among the rows meeting cond
      n_species <- function(cond) {
            keep <- cond & !is.na(f$species)
            p <- plot[keep]
            tabulate(p[!duplicated(data.frame(p, f$species[keep]))], nbins = n)
      }

      tibble::tibble(
            x = f$x[first], y = f$y[first],
            samp_method_cd = sums(f$samp_method_cd) / tabulate(plot, nbins = n),
            pres = any_of(focal), # focal species presence/absence (any life stage)
            pres_seed = any_of(focal & f$n_seedlings > 0), # focal species seedling presence
            pres_sap = any_of(focal & f$n_saplings > 0), # focal species sapling presence
            pres_tree = any_of(focal & f$n_trees > 0), # focal species adult (tree) presence
            batot = sums(f$basal_area * focal), # total focal species basal area
            baprop = batot / sums(f$basal_area), # proportion ba
            baplot = sums(f$basal_area), # total plot ba, all species
            nttot = sums(f$n_trees * focal), # focal species tree count
            ntplot = sums(f$n_trees), # total plot tree count, all species
            satot = sums(f$n_saplings * focal), # focal species sapling count
            saplot = sums(f$n_saplings), # total plot sapling count, all species
            nstot = sums(f$n_seedlings * focal), # focal species seedling count
            nsplot = sums(f$n_seedlings), # total plot seedling count, all species
            # Plot species richness, by stage. A count of distinct species is
            # an aggregation ACROSS species done at the plot, exactly like
            # baplot and ntplot -- which is what makes it an All-trees variable
            # rather than something only the multi-entity path can produce.
            # Because AIM is a weighted mean (linear) with entity-independent
            # weights, generalizing this equals summing the per-entity
            # occurrence surfaces. The two routes should agree to floating
            # point, which makes each a check on the other.
            nspp      = n_species(TRUE),
            nspp_seed = n_species(f$n_seedlings > 0),
            nspp_sap  = n_species(f$n_saplings > 0),
            nspp_tree = n_species(f$n_trees > 0),
            forest_type = f$forest_type[first]) |> # plot-level FORTYPCD
            dplyr::left_join(ft_to_fg, by = "forest_type") |>
            dplyr::mutate(baprop = ifelse(baplot == 0, 0, baprop),
                          # Per-plot regeneration fractions (focal species), on an
                          # occurrence basis to match the modeled surface (see
                          # ts_ratio_pieces()): 0 = adults only, 0.5 = stage and
                          # adults both present, 1 = stage only. NA where the
                          # species is absent as both that stage and adults.
                          seed_frac = ifelse((nstot + nttot) > 0,
                                             (nstot > 0) / ((nstot > 0) + (nttot > 0)), NA_real_),
                          sap_frac  = ifelse((satot + nttot) > 0,
                                             (satot > 0) / ((satot > 0) + (nttot > 0)), NA_real_))
}

#' The single response, aligned with `ts_plots()` rows
#'
#' @return A numeric vector, or NULL for a variable the track doesn't offer.
ts_response <- function(plots, track, fia_var, forest_type = NULL, forest_group = NULL) {
      po <- plots
      if (track == "Species") {
            switch(fia_var,
                   "presence probability"  = po$pres,
                   "seedling presence"     = po$pres_seed,
                   "sapling presence"      = po$pres_sap,
                   "tree presence"         = po$pres_tree,
                   "proportion basal area" = po$baprop,
                   "total basal area"      = po$batot,
                   "number of trees"       = po$nttot,
                   "number of saplings"    = po$satot,
                   "number of seedlings"   = po$nstot,
                   "seedling fraction"     = po$seed_frac,
                   "sapling fraction"      = po$sap_frac)
      } else if (track == "Forest type") {
            if (is.null(forest_type)) stop("No forest type given")
            as.integer(!is.na(po$forest_type) & po$forest_type == as.integer(forest_type))
      } else if (track == "Forest group") {
            if (is.null(forest_group)) stop("No forest group given")
            as.integer(!is.na(po$forest_group) & po$forest_group == as.integer(forest_group))
      } else if (track == "All trees") {
            switch(fia_var,
                   "presence probability" = as.integer((po$ntplot + po$saplot) > 0 | po$nsplot > 0),
                   "seedling presence"    = as.integer(po$nsplot > 0),
                   "sapling presence"     = as.integer(po$saplot > 0),
                   "tree presence"        = as.integer(po$ntplot > 0),
                   "total basal area"     = po$baplot,
                   "number of trees"      = po$ntplot,
                   "number of saplings"   = po$saplot,
                   "number of seedlings"  = po$nsplot,
                   "species richness"          = po$nspp,
                   "seedling species richness" = po$nspp_seed,
                   "sapling species richness"  = po$nspp_sap,
                   "adult species richness"    = po$nspp_tree,
                   # Occurrence-based, matching ts_ratio_pieces() and the
                   # Species-track fractions in ts_plots(): 0 = adults only,
                   # 0.5 = stage and adults both present, 1 = stage only.
                   "seedling fraction"    = ifelse((po$nsplot + po$ntplot) > 0,
                                                   (po$nsplot > 0) / ((po$nsplot > 0) + (po$ntplot > 0)), NA_real_),
                   "sapling fraction"     = ifelse((po$saplot + po$ntplot) > 0,
                                                   (po$saplot > 0) / ((po$saplot > 0) + (po$ntplot > 0)), NA_real_))
      }
}

#' Numerator and denominator of a pooled-ratio variable, per plot
#'
#' Pooled-ratio variables (proportion basal area, the regeneration fractions)
#' are modeled as a ratio of means, sum(w * num) / sum(w * denom), not a mean
#' of per-plot ratios. Occurrence is denom > 0.
#'
#' @return list(num, denom), or NULL for a variable that isn't a pooled ratio.
ts_ratio_pieces <- function(plots, track, fia_var) {
      po <- plots
      # Regeneration fractions are a fraction of OCCURRENCE weight, not of
      # stem density. Seedlings and saplings are tallied on the microplot
      # (each stem = 74.965/acre) while adults come off the subplot (each =
      # 6.018/acre), so a density ratio applies a ~12.5x amplification to
      # every regen stem before dividing: one seedling + one adult reads
      # 0.93 rather than 0.5, squashing the balanced core up against the
      # pure-colonization end. Using presence indicators puts both stages on
      # the same footing and gives a readable axis: 0 = adults only
      # (relict/trailing), ~0.5 = both present (core), 1 = regen only
      # (colonizing/leading).
      #
      # Proportion basal area is deliberately left on a density basis: it is
      # focal BA over stand BA, with each stem's own TPA applied to both
      # sides, so it never crosses the microplot/subplot boundary and isn't
      # distorted the way the regen fractions are.
      if (track == "Species") {
            switch(fia_var,
                   "proportion basal area" = list(num = po$batot, denom = po$baplot),
                   "seedling fraction"     = list(num = as.integer(po$nstot > 0),
                                                  denom = as.integer(po$nstot > 0) +
                                                        as.integer(po$nttot > 0)),
                   "sapling fraction"      = list(num = as.integer(po$satot > 0),
                                                  denom = as.integer(po$satot > 0) +
                                                        as.integer(po$nttot > 0)),
                   NULL)
      } else if (track == "All trees") {
            switch(fia_var,
                   "seedling fraction" = list(num = as.integer(po$nsplot > 0),
                                              denom = as.integer(po$nsplot > 0) +
                                                    as.integer(po$ntplot > 0)),
                   "sapling fraction"  = list(num = as.integer(po$saplot > 0),
                                              denom = as.integer(po$saplot > 0) +
                                                    as.integer(po$ntplot > 0)),
                   NULL)
      } else NULL
}

#' Per-entity responses for a community run: one column per entity
#'
#' One row per plot, row-aligned with `plots`, holding the response each entity
#' would have as the focal one. Alignment is the whole ballgame: a mismatch
#' silently attaches every entity's response to the wrong plot and gives
#' plausible-looking maps. It is guaranteed by construction (the join runs from
#' the plots' coordinates) and asserted on the coordinates themselves.
#'
#' @return A data frame with an `n_plots` attribute: each entity's number of
#'   occupied plots.
ts_multi_response <- function(data, plots, pool_bbox, track, fia_var) {
      po <- plots
      xy <- po[, c("x", "y")]

      if (track %in% c("Forest type", "Forest group")) {
            # Membership is plot-level -- one type per plot -- so the wide
            # form is an indicator matrix and fia_var does not apply.
            key <- if (track == "Forest type") po$forest_type else po$forest_group
            lvls <- sort(unique(key[!is.na(key)]))
            if (!length(lvls)) stop("No forest types among the plots")
            w <- as.data.frame(
                  vapply(lvls, function(k) as.integer(!is.na(key) & key == k),
                         integer(nrow(po))),
                  check.names = FALSE)

            # Columns are named for the classes rather than their FORTYPCD
            # integers: the names ARE the entity labels everywhere downstream.
            nm <- if (track == "Forest type") ts_forest_type_name(data, as.integer(lvls))
            else ts_forest_group_name(data, as.integer(lvls))
            nm[is.na(nm)] <- paste(track, lvls[is.na(nm)])
            names(w) <- make.unique(as.character(nm))

      } else {
            if (!fia_var %in% multi_vars)
                  stop("'", fia_var, "' has no per-entity form for a community run")

            # Same summaries ts_plots() computes, but grouped by entity as
            # well as plot.
            long <- fia_in_box(data, pool_bbox) |>
                  dplyr::filter(!is.na(species)) |>
                  dplyr::group_by(x, y, species) |>
                  dplyr::summarize(
                        # Presence is "the entity was recorded at this
                        # plot", matching ts_plots()'s any(focal) -- it is
                        # not conditioned on any stage count.
                        pres      = 1L,
                        pres_seed = as.integer(any(n_seedlings > 0)),
                        pres_sap  = as.integer(any(n_saplings > 0)),
                        pres_tree = as.integer(any(n_trees > 0)),
                        batot     = sum(basal_area),
                        nttot     = sum(n_trees),
                        satot     = sum(n_saplings),
                        nstot     = sum(n_seedlings),
                        .groups = "drop")

            val_col <- switch(fia_var,
                              "presence probability"  = "pres",
                              "seedling presence"     = "pres_seed",
                              "sapling presence"      = "pres_sap",
                              "tree presence"         = "pres_tree",
                              "total basal area"      = "batot",
                              "proportion basal area" = "batot",
                              "number of trees"       = "nttot",
                              "number of saplings"    = "satot",
                              "number of seedlings"   = "nstot")

            wide <- long |>
                  dplyr::select(x, y, species, value = dplyr::all_of(val_col)) |>
                  tidyr::pivot_wider(names_from = species, values_from = value,
                                     values_fill = 0)

            # Plots where no entity was recorded were dropped by the
            # !is.na(species) filter and come back as all-zero rows: absence
            # is data, and dropping those plots would bias every entity's
            # surface upward.
            j <- dplyr::left_join(xy, wide, by = dplyr::join_by(x, y))
            stopifnot(nrow(j) == nrow(xy),
                      identical(j$x, xy$x), identical(j$y, xy$y))

            w <- as.data.frame(j[, setdiff(names(j), c("x", "y")), drop = FALSE],
                               check.names = FALSE)
            w[is.na(w)] <- 0

            if (identical(fia_var, "proportion basal area")) {
                  # The denominator is the plot total across ALL entities, so
                  # it is shared. Plots with no basal area get 0 rather than
                  # NaN, matching ts_plots()'s baprop guard.
                  denom <- po$baplot
                  w[] <- lapply(w, function(v) ifelse(denom == 0, 0, v / denom))
            }
      }

      # Occupied-plot count per entity: the pool-level support signal
      attr(w, "n_plots") <- vapply(w, function(v) sum(v > 0, na.rm = TRUE), integer(1))
      w
}
