# Analog impact surfaces --------------------------------------------------------
#
# AIM surfaces are computed on plot and cell coordinates in CONUS Albers
# equal-area (EPSG:5070), in km so the geographic bandwidth keeps its units. The
# FFT method needs planar coordinates, and planar distance also makes the exact
# method about 1.5x faster than great-circle distance; over CONUS the two differ
# by about 1%. Results are written back onto the lon/lat modeling grid.
aim_crs <- "EPSG:5070"

aim_xy <- function(xy) {
      p <- terra::project(as.matrix(xy[, c("x", "y")]), "EPSG:4326", aim_crs) / 1000
      colnames(p) <- c("x", "y")
      p
}

aim_kernels <- function(theta_clim, theta_geog) {
      list(env  = analogs::kernel(max = theta_clim * 3, weight = "gaussian", theta = theta_clim),
           geog = analogs::kernel(max = theta_geog * 3, weight = "gaussian", theta = theta_geog))
}

# FFT runs are tiled geographically, so each tile's grid spans only the climates
# found there. ~300 km tiles were fastest in benchmarks.
aim_tile_km <- 300

# Which method a run uses: "exact" visits every cell-plot pair, "fft" (binned
# convolution) costs about the same whatever the grain, and its error is small
# (99.9th percentile ~0.007 on a 0-1 occurrence scale; identical NA cells). Exact
# wins at coarse grain or narrow bandwidths, FFT at fine grain and wide ones.
#
# The choice must be a deterministic function of the settings, so that the same
# settings always give the same map -- hence estimated work, not timings. Exact
# work is the number of cell-plot pairs within the kernel limits, from analog
# counts at a fixed sample of cells. FFT work is the number of grid nodes, from
# each tile's extent in geography and climate (nodes are a quarter-bandwidth
# apart, plus a 3-bandwidth halo), times the number of kernel sums. The
# coefficients were fitted on runs at 1-10 km over 6-20 degree domains with
# 1-80 outcome columns; only their ratio matters, so machine speed cancels.
# One method serves every timeframe of a run, so differences between timeframes
# never mix methods.
aim_method <- function(pts_list, index, theta_clim, theta_geog, n_y, n_sums) {
      env_vars <- setdiff(names(pts_list[[1]]), c("x", "y"))
      if (length(env_vars) > 3) return("exact")   # FFT supports up to 3

      est <- vapply(pts_list, function(pts) {
            n <- nrow(pts)
            probe <- pts[unique(round(seq(1, n, length.out = min(n, 2000)))), , drop = FALSE]
            count <- suppressWarnings(analogs::analog_search(
                  x = probe, pool = index, select = "all", stat = "count",
                  env = analogs::kernel(max = theta_clim * 3),
                  geog = analogs::kernel(max = theta_geog * 3)))$count
            pairs <- n * mean(count, na.rm = TRUE)

            spans <- pts |>
                  dplyr::group_by(tx = floor(x / aim_tile_km), ty = floor(y / aim_tile_km)) |>
                  dplyr::summarise(dplyr::across(dplyr::everything(), ~ max(.x) - min(.x)), .groups = "drop")
            nodes_along <- function(span, theta) (span + 6 * theta) / (0.25 * theta) + 1
            nodes <- sum(nodes_along(spans$x, theta_geog) * nodes_along(spans$y, theta_geog) *
                         Reduce(`*`, lapply(spans[env_vars], nodes_along, theta = theta_clim)))

            c(exact = 5.0e-6 * n + 4.5e-9 * pairs * (1 + 0.20 * (n_y - 1)),
              fft   = 1.8e-6 * n + 2.1e-9 * nodes * (1 + 0.42 * n_sums))
      }, numeric(2))
      if (sum(est["fft", ]) < sum(est["exact", ])) "fft" else "exact"
}

# AIM surfaces for several focal rasters (e.g. historical and future climate)
# against one pool, on each raster's own grid. `pool` is pool()'s list, of which
# this uses the AIM-coordinate plot table (`aim_pts`) and its index
# (`aim_index`). Returns a list of SpatRasters named like `focals`, with the
# layers analog_impact() would give for a raster.
aim_surfaces <- function(focals, pool, y, theta_clim, theta_geog, stat) {
      k <- aim_kernels(theta_clim, theta_geog)
      cells <- lapply(focals, function(r) as.data.frame(r, xy = TRUE, cells = TRUE, na.rm = TRUE))
      pts <- lapply(cells, function(fd) data.frame(aim_xy(fd), fd[, setdiff(names(fd), c("cell", "x", "y")), drop = FALSE]))

      run <- function(method) lapply(pts, function(p) {
            if (method == "exact") {
                  suppressWarnings(analogs::analog_impact(
                        x = p, pool = pool$aim_index, env = k$env, geog = k$geog,
                        stat = stat, y = y, coord_type = "projected"))
            } else {
                  span <- vapply(p[, c("x", "y")], function(v) diff(range(v)), numeric(1))
                  suppressWarnings(analogs::tiled_analog_search(
                        x = p, pool = pool$aim_pts, y = y,
                        n_tiles = max(1, round(prod(span) / aim_tile_km^2)),
                        fun = analogs::analog_impact, env = k$env, geog = k$geog,
                        stat = stat, method = "fft", coord_type = "projected",
                        progress = FALSE))
            }
      })

      # Kernel sums FFT computes: one per outcome column (two with
      # weighted_sd), the weights, and the squared weights if a stat needs them
      n_sums <- ncol(y) * (1 + ("weighted_sd" %in% stat)) + 1 +
            any(c("ess", "effective_weight") %in% stat)
      method <- aim_method(pts, pool$aim_index, theta_clim, theta_geog,
                           n_y = ncol(y), n_sums = n_sums)
      res <- if (method == "fft") tryCatch(run("fft"), error = function(e) {
            # A tile whose FFT grid exceeds the size limit (narrow bandwidths
            # over climatically varied ground) sends the whole run to exact,
            # rather than mixing methods between tiles. Any other error is real.
            if (!grepl("max_cells", conditionMessage(e))) stop(e)
            message("AIM: FFT grid too large for a tile; using the exact method")
            NULL
      })
      if (is.null(res)) res <- run("exact")

      Map(function(r, fd, out) {
            cols <- setdiff(names(out), c("index", "x", "y"))
            s <- terra::rast(r, nlyrs = length(cols), names = cols, vals = NA)
            s[fd$cell] <- as.matrix(out[, cols, drop = FALSE])
            s
      }, focals, cells, res)
}

# The analog pool ---------------------------------------------------------------

#' The outcome columns a pool carries for a single-entity run
#'
#' One `response` column, or for a pooled-ratio variable its occurrence,
#' numerator and denominator, which `ts_surfaces()` combines after
#' generalizing each. `samp_method_cd` rides along for plot styling and is not
#' an outcome.
ts_pool_values <- function(plots, response = NULL, ratio = NULL) {
      if (is.null(ratio)) {
            if (is.null(response)) stop("Give a response or ratio pieces")
            tibble::tibble(response = response,
                           samp_method_cd = plots$samp_method_cd)
      } else {
            tibble::tibble(occ = as.integer(ratio$denom > 0),
                           num = ratio$num,
                           denom = ratio$denom,
                           samp_method_cd = plots$samp_method_cd)
      }
}

#' The analog pool: the plots with their climates, indexed for search
#'
#' Plot climates come from the native-resolution rasters at each plot's own
#' location, so the pool doesn't depend on the grid's grain or crop.
#'
#' @param climate `ts_climate()`'s list.
#' @param values Outcome columns row-aligned with `plots` (`ts_pool_values()`),
#'   or NULL for a community run, whose responses are passed to
#'   `ts_multi_surfaces()` directly.
ts_pool <- function(plots, climate, values = NULL) {
      xy <- plots[, c("x", "y")]
      clim_hst <- terra::extract(climate$hst, xy, ID = FALSE, method = "bilinear")
      clim_fut <- terra::extract(climate$fut, xy, ID = FALSE, method = "bilinear")

      # The same plots, in the same rows, in AIM coordinates (see
      # aim_surfaces()), as a table for the FFT method and an index for the
      # exact one. index_hst/index_fut stay lon/lat for the focal-site analog
      # search.
      aim_pts <- data.frame(aim_xy(xy), clim_hst)

      list(index_hst = analogs::build_analog_index(cbind(xy, clim_hst), coord_type = "lonlat",
                                                   geog_res_adj = 4),
           index_fut = analogs::build_analog_index(cbind(xy, clim_fut), coord_type = "lonlat",
                                                   geog_res_adj = 4),
           aim_pts = aim_pts,
           aim_index = analogs::build_analog_index(aim_pts, coord_type = "projected",
                                                   geog_res_adj = 4),
           clim_hst = clim_hst,
           clim_fut = clim_fut,
           values = values)
}


# Surfaces ---------------------------------------------------------------------

# Stats that describe the WEIGHTS rather than an outcome column: support (sum
# of weights) and closeness (sum of squared weights over sum of weights).
# analog_impact() returns one layer for each however many outcome columns `y`
# carries.
weight_stats <- c("sum_weights", "effective_weight")

# analog_impact() names outcome layers `weighted_mean_<column>` when `y` has
# two or more columns but bare `weighted_mean` when it has one, so stripping
# the prefix is unreliable. Rename positionally instead: outcome layers are
# whatever isn't a weight stat, in y-column order.
name_outcomes <- function(r, y_cols, other = character(0)) {
      names(r)[names(r) == "weighted_sd"] <- "spread"
      is_outcome <- !names(r) %in% c(other, "spread")
      stopifnot(sum(is_outcome) == length(y_cols))
      names(r)[is_outcome] <- y_cols
      r
}

#' Single-entity surfaces for both timeframes
#'
#' @param grid `ts_climate_grid()`'s list.
#' @param pool `ts_pool()`'s list, with `values` from `ts_pool_values()`.
#' @return list(bsl, fut) of SpatRasters with layers `response`, the weight
#'   stats `sum_weights` and `effective_weight`, `spread` (the weighted SD of
#'   the response, single non-ratio variables only), and for a pooled ratio its
#'   `occ`, `num` and `denom`.
ts_surfaces <- function(grid, pool, fia_var, theta_clim, theta_geog) {
      y_cols <- setdiff(names(pool$values), "samp_method_cd")
      is_ratio <- "denom" %in% y_cols

      # Spread among the analog plots, for a single (non-ratio) response only:
      # a pooled ratio is not the mean of plot values, so their plot-to-plot
      # SD isn't its spread.
      spread_stat <- if (!is_ratio) "weighted_sd"

      surf <- aim_surfaces(list(bsl = grid$hst, fut = grid$fut), pool,
                           y = pool$values[, y_cols, drop = FALSE],
                           theta_clim = theta_clim, theta_geog = theta_geog,
                           stat = c("weighted_mean", spread_stat, weight_stats))
      surf <- lapply(surf, name_outcomes, y_cols = y_cols, other = weight_stats)

      # Pooled-ratio surface: ratio of the num and denom means, as the
      # standard `response` layer.
      #
      # num/denom is already conditional on occurrence -- absent plots
      # contribute 0 to both and drop out -- so it reads "the ratio among
      # plots where the focal entity occurs".
      #
      # The regeneration fractions stay conditional and are MASKED (NA) where
      # no occupied plot falls within kernel range. Multiplying by occurrence
      # instead would drive unsuitable ground to 0 -- but 0 already means
      # "occurs, adults only" (relict), so the product would conflate two very
      # different states. Masking keeps them distinct: blank = doesn't occur,
      # 0 = occurs but no regeneration. This also matches the per-plot values,
      # NA where the species is absent and 0 where it's adults-only.
      #
      # Proportion basal area keeps the occurrence weighting: there 0 means
      # "the focal species is 0% of a real stand" vs "there is no stand" --
      # both read as "the species isn't here", so collapsing them costs
      # little, and the unconditional surface is better behaved.
      if (is_ratio) {
            conditional <- grepl("fraction", fia_var)
            combine_ratio <- function(r) {
                  resp <- r[["num"]] / r[["denom"]]
                  if (!conditional) resp <- r[["occ"]] * resp
                  # denom == 0: no occupied plot in range. Undefined for the
                  # fractions; 0 for the occurrence-weighted proportion.
                  resp[!is.finite(resp[])] <- if (conditional) NA else 0
                  names(resp) <- "response"
                  c(r, resp)
            }
            surf <- lapply(surf, combine_ratio)
      }
      surf[c("fut", "bsl")]
}

#' Per-entity surfaces for a community run, both timeframes
#'
#' One analog_impact() run per timeframe with every entity as an outcome
#' column; kernel weights are computed once and reused across columns. Cost is
#' linear in cells x entities (see `ts_community_cost()`).
#'
#' @param y `ts_multi_response()`'s data frame.
#' @return list(bsl, fut) of SpatRasters: one layer per entity, named for it,
#'   plus `sum_weights`.
ts_multi_surfaces <- function(grid, pool, y, theta_clim, theta_geog) {
      if (!ncol(y)) stop("No entities to model")
      surf <- aim_surfaces(list(bsl = grid$hst, fut = grid$fut), pool, y = y,
                           theta_clim = theta_clim, theta_geog = theta_geog,
                           stat = c("weighted_mean", "sum_weights"))
      lapply(surf, name_outcomes, y_cols = names(y), other = "sum_weights")
}

#' Size and rough duration of a community run, before running it
#'
#' Memory is the per-entity surfaces for both timeframes. Timing coefficients
#' were measured: about 0.10 s fixed plus 0.012 s per entity at 31,624 cells,
#' per timeframe, scaling linearly in cells; it runs a little optimistic above
#' ~25 entities.
#'
#' @param max_gb The memory ceiling, for the smallest grain that fits it.
ts_community_cost <- function(n_cells, n_entities, clim_agg, max_gb) {
      gb <- n_cells * n_entities * 8 * 2 / 1024^3
      secs <- (n_cells / 31624) * (0.10 + 0.012 * n_entities) * 2

      # Smallest grain that would fit. Cells scale as the inverse square of
      # grain, so gb(g) = gb_now * (agg_now / g)^2.
      fits <- ceiling(clim_agg / sqrt(max_gb / gb))
      list(entities = n_entities, cells = n_cells, gb = gb, secs = secs,
           # NA when even the coarsest grain offered (10) is not enough --
           # then the bounding box is the thing to change.
           grain_that_fits = if (is.finite(fits) && fits <= 10) max(fits, 1) else NA_integer_)
}
