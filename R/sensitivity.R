# Focal-site sensitivity ----------------------------------------------------------
#
# How much a site's estimate changes when one choice at a time is varied from
# the current settings: the climate model, the emissions pathway, either
# bandwidth, and the sample of FIA plots. The panel makes no probability
# statement and doesn't combine the rows, which aren't on a common scale.

#' How much a site's estimate changes under alternative choices
#'
#' Every estimate is a weighted mean (or pooled ratio) over a local pool of
#' plots around the site, with weights computed here from the plots' geographic
#' and climatic distances to the site's exact climate -- the same Gaussian
#' kernels, cut off at three bandwidths, that `aim_kernels()` gives
#' `analogs::analog_impact()` (a test holds them equal). Computing the weights
#' once per climate lets every bandwidth reuse one distance search, and gives
#' the bootstrap its fixed weights. The estimate can differ slightly from the
#' mapped cell's, which comes from the cell's aggregated climate.
#'
#' @param settings `ts_settings()` for a single-entity run (a storyline is
#'   resolved first).
#' @param data `ts_data()`.
#' @param site c(x, y), lon/lat.
#' @param mult Bandwidth multipliers; the evaluation module's bracket. Must
#'   include 1.
#'
#' The climate-variable row swaps each selected variable for each unused one.
#' Keeping the number of variables fixed keeps the climate bandwidth roughly
#' comparable, since it is a distance across all of them; it isn't exact,
#' because correlated variables measure distance differently.
#' @param n_boot Bootstrap draws for the plot sample.
#' @return A tibble, one row per (row, alternative): `row` ("gcm", "ssp",
#'   "theta_geog", "theta_clim", "clim_vars" or "sample"), `alternative` (a
#'   model, a pathway, a bandwidth, or comma-separated climate variables, as
#'   character), `current` (whether it is the current setting), `hst`, `fut`
#'   and `delta` (fut - hst), `occ_hst` and `occ_fut` (the occurrence
#'   probability a pooled ratio is conditional on, else NA), the evidence
#'   behind each estimate as the map's layers measure it (`support_*`, the sum
#'   of weights, and `closeness_*`, the sum of squared weights over the sum of
#'   weights), `n_plots` (plots with nonzero weight at either timeframe), and
#'   for the
#'   "sample" row the 5%-95% bootstrap bounds `hst_lower`, `hst_upper`, ...,
#'   `delta_upper`. Estimates are NA where no plot is in range. Attribute
#'   "missing" names the alternative climate files that aren't available.
ts_site_sensitivity <- function(settings, data, site, mult = c(.25, .5, 1, 2, 4),
                                n_boot = 500) {
      s <- ts_resolve_scenario(settings, data)
      if (ts_multi(s)) stop("Sensitivity is for single-entity runs")
      if (!1 %in% mult) stop("`mult` must include 1, the current bandwidth")
      gcm_now <- ts_scenario_model(s)

      # The local pool: every plot the widest kernel could reach, independent
      # of the bounding box. Km to degrees, with longitude's shrinking with
      # latitude.
      reach <- 3 * max(mult) * s$theta_geog
      pb <- list(x = site[1] + c(-1, 1) * reach / (111 * cos(site[2] * pi / 180)),
                 y = site[2] + c(-1, 1) * reach / 111)
      plots <- ts_plots(data, pb, s$species)
      ratio <- ts_ratio_pieces(plots, s$track, s$fia_var)
      response <- if (is.null(ratio))
            ts_response(plots, s$track, s$fia_var, s$forest_type, s$forest_group)
      if (is.null(ratio) && is.null(response))
            stop("The ", s$track, " track has no variable '", s$fia_var, "'")
      y <- ts_pool_values(plots, response, ratio)
      y <- as.matrix(y[, setdiff(names(y), "samp_method_cd")])
      xy <- plots[, c("x", "y")]

      # Climate-variable sets: each selected variable swapped for each unused
      # one, so the number of variables stays the same
      swaps <- tidyr::expand_grid(out = s$clim_vars, into = setdiff(clim_vars(), s$clim_vars))
      all_vars <- c(s$clim_vars, unique(swaps$into))
      hst <- terra::rast(data$clim_historical)[[all_vars]]
      plot_clim <- terra::extract(hst, xy, ID = FALSE, method = "bilinear")
      pool <- ts_aim_pool(xy, plot_clim[s$clim_vars])

      # The site's climate: historical, and future under each model and
      # pathway for the selected period
      site_xy <- data.frame(x = site[1], y = site[2])
      site_climate <- function(r) data.frame(aim_xy(site_xy),
                                             terra::extract(r, site_xy, ID = FALSE, method = "bilinear"))
      futures <- tidyr::expand_grid(ssp = sort(unique(data$scenarios$ssp)),
                                    model = sort(unique(data$scenarios$model)))
      futures <- futures[futures$ssp == s$scenario_ssp | futures$model == gcm_now, ]
      futures$path <- ts_scenario_path(data, futures$ssp, s$scenario_period, futures$model)
      missing <- basename(futures$path[!file.exists(futures$path)])
      futures <- futures[file.exists(futures$path), ]
      now <- futures$ssp == s$scenario_ssp & futures$model == gcm_now
      if (!any(now))
            stop("No climate data for ", basename(ts_scenario_path(data, s$scenario_ssp,
                                                                   s$scenario_period, gcm_now)))

      # Each climate's plot distances, by default to the widest kernels' reach
      distances <- function(q, index = pool$aim_index, vars = s$clim_vars, m = max(mult)) {
            q <- q[c("x", "y", vars)]
            if (anyNA(q)) return(NULL)   # off the climate grid
            suppressWarnings(analogs::analog_search(
                  x = q, pool = index, select = "all", stat = NULL,
                  env = analogs::kernel(max = 3 * m * s$theta_clim),
                  geog = analogs::kernel(max = 3 * m * s$theta_geog), coord_type = "projected"))
      }
      site_hst <- site_climate(hst)
      site_fut_now <- site_climate(terra::rast(futures$path[now])[[all_vars]])
      d_hst <- distances(site_hst)
      d_fut <- lapply(futures$path, function(p) distances(site_climate(terra::rast(p)[[s$clim_vars]])))
      d_fut_now <- d_fut[[which(now)]]

      # Kernel weights per pool plot (0 out of range)
      weights <- function(d, theta_geog = s$theta_geog, theta_clim = s$theta_clim) {
            w <- numeric(nrow(y))
            if (is.null(d)) return(w)
            k <- exp(-0.5 * (d$geog_dist / theta_geog)^2) * exp(-0.5 * (d$env_dist / theta_clim)^2)
            k[d$geog_dist > 3 * theta_geog | d$env_dist > 3 * theta_clim] <- 0
            w[d$analog_index] <- k
            w
      }

      # Estimates from weighted column sums: `ws` is a matrix of weight
      # vectors (one column per estimate) over the pool plots `at`; returns
      # the estimate and, for a pooled ratio, its occurrence, one per column
      estimate <- function(ws, at = seq_len(nrow(y))) {
            sw <- colSums(ws)
            m <- crossprod(ws, y[at, , drop = FALSE]) / sw   # estimates x outcome columns
            m[sw == 0, ] <- NA
            if (is.null(ratio)) return(list(value = m[, "response"], occ = rep(NA_real_, ncol(ws))))
            list(value = ts_combine_ratio(m[, "occ"], m[, "num"], m[, "denom"], s$fia_var),
                 occ = m[, "occ"])
      }

      rows <- function(row, alternative, current, w_hst, w_fut) {
            h <- estimate(w_hst)
            f <- estimate(w_fut)
            support <- function(w) unname(colSums(w))
            closeness <- function(w) unname(colSums(w^2) / colSums(w))
            tibble::tibble(row = row, alternative = as.character(alternative), current = current,
                           hst = unname(h$value), fut = unname(f$value),
                           delta = unname(f$value - h$value),
                           occ_hst = unname(h$occ), occ_fut = unname(f$occ),
                           support_hst = support(w_hst), support_fut = support(w_fut),
                           closeness_hst = closeness(w_hst), closeness_fut = closeness(w_fut),
                           n_plots = unname(colSums(w_hst > 0 | w_fut > 0)))
      }
      w_hst <- weights(d_hst)
      w_fut <- weights(d_fut_now)
      same <- function(n) matrix(w_hst, nrow(y), n)
      by_scenario <- function(keep) {
            f <- futures[keep, ]
            list(f = f, w = vapply(d_fut[keep], weights, numeric(nrow(y))))
      }
      gcm <- by_scenario(futures$ssp == s$scenario_ssp)
      ssp <- by_scenario(futures$model == gcm_now)
      geog <- lapply(mult, function(m) list(h = weights(d_hst, theta_geog = s$theta_geog * m),
                                            f = weights(d_fut_now, theta_geog = s$theta_geog * m)))
      clim <- lapply(mult, function(m) list(h = weights(d_hst, theta_clim = s$theta_clim * m),
                                            f = weights(d_fut_now, theta_clim = s$theta_clim * m)))
      # The current set reuses the current weights, so its estimate is
      # identical to the other rows'
      vars <- lapply(seq_len(nrow(swaps)), function(i)
            replace(s$clim_vars, s$clim_vars == swaps$out[i], swaps$into[i]))
      sets <- lapply(vars, function(v) {
            index <- ts_aim_pool(xy, plot_clim[v])$aim_index
            list(h = weights(distances(site_hst, index, v, 1)),
                 f = weights(distances(site_fut_now, index, v, 1)))
      })
      sets <- c(list(list(h = w_hst, f = w_fut)), sets)
      vars <- c(list(s$clim_vars), vars)
      pick <- function(l, part) vapply(l, `[[`, numeric(nrow(y)), part)

      out <- dplyr::bind_rows(
            rows("gcm", gcm$f$model, gcm$f$model == gcm_now, same(nrow(gcm$f)), gcm$w),
            rows("ssp", ssp$f$ssp, ssp$f$ssp == s$scenario_ssp, same(nrow(ssp$f)), ssp$w),
            rows("theta_geog", s$theta_geog * mult, mult == 1, pick(geog, "h"), pick(geog, "f")),
            rows("theta_clim", s$theta_clim * mult, mult == 1, pick(clim, "h"), pick(clim, "f")),
            rows("clim_vars", vapply(vars, paste, "", collapse = ", "), seq_along(vars) == 1,
                 pick(sets, "h"), pick(sets, "f")),
            dplyr::bind_cols(rows("sample", "current", TRUE, as.matrix(w_hst), as.matrix(w_fut)),
                      ts_site_bootstrap(w_hst, w_fut, estimate, n_boot,
                                        seed = list(site, s[surfaces_setting_ids], n_boot))))
      structure(out, missing = missing)
}

#' The plot-sample row: a bootstrap over the plots in range
#'
#' Weights depend only on distances, so they stay fixed while the plots are
#' resampled with replacement from those weighted at either timeframe. Both
#' timeframes share each draw, so the interval on the change reflects their
#' shared plots. Seeded from `seed` (hashed), so the same site and settings
#' always give the same interval, without touching the session's random
#' numbers.
#'
#' @param estimate `ts_site_sensitivity()`'s estimator: a matrix of weight
#'   vectors and the pool plots they weight.
#' @return A one-row tibble of 5% and 95% bounds for hst, fut and delta.
ts_site_bootstrap <- function(w_hst, w_fut, estimate, n_boot, seed) {
      in_range <- which(w_hst > 0 | w_fut > 0)
      n <- length(in_range)
      bounds <- function(x) stats::quantile(x, c(.05, .95), na.rm = TRUE, names = FALSE)
      if (n == 0) {
            b <- list(hst = c(NA, NA), fut = c(NA, NA), delta = c(NA, NA))
      } else {
            seed <- strtoi(substr(digest::digest(seed, algo = "xxhash32"), 1, 7), 16L)
            draws <- withr::with_seed(seed, stats::rmultinom(n_boot, n, rep(1, n)))
            # Each draw's weights over the plots in range: a plot drawn k
            # times counts k times
            h <- estimate(draws * w_hst[in_range], in_range)$value
            f <- estimate(draws * w_fut[in_range], in_range)$value
            b <- list(hst = bounds(h), fut = bounds(f), delta = bounds(f - h))
      }
      tibble::tibble(hst_lower = b$hst[1], hst_upper = b$hst[2],
                     fut_lower = b$fut[1], fut_upper = b$fut[2],
                     delta_lower = b$delta[1], delta_upper = b$delta[2])
}
