# Focal-site sensitivity (R/sensitivity.R)

withr::with_dir(app_dir, source("scripts/load_core.R"))
data <- withr::with_dir(app_dir, ts_data("data"))

# Near Flagstaff, AZ, where ponderosa pine is common
site <- c(-111.7, 35.2)
pipo <- function(...) ts_settings(species = "Pinus ponderosa", ...)

# The site's estimate straight from analogs::analog_impact(), over the same
# local pool, at the given bandwidths: the reference the hand-computed kernel
# weights must reproduce
reference <- function(settings, theta_geog = settings$theta_geog, theta_clim = settings$theta_clim) {
      s <- ts_resolve_scenario(settings, data)
      reach <- 3 * 4 * s$theta_geog
      pb <- list(x = site[1] + c(-1, 1) * reach / (111 * cos(site[2] * pi / 180)),
                 y = site[2] + c(-1, 1) * reach / 111)
      plots <- ts_plots(data, pb, s$species)
      ratio <- ts_ratio_pieces(plots, s$track, s$fia_var)
      response <- if (is.null(ratio)) ts_response(plots, s$track, s$fia_var)
      climate <- ts_climate(data, s$clim_vars, s$scenario_ssp, s$scenario_period, s$scenario_model)
      pool <- ts_pool(plots, climate, ts_pool_values(plots, response, ratio))
      cols <- setdiff(names(pool$values), "samp_method_cd")
      k <- aim_kernels(theta_clim, theta_geog)
      xy <- data.frame(x = site[1], y = site[2])
      out <- lapply(list(hst = climate$hst, fut = climate$fut), function(r) {
            q <- data.frame(aim_xy(xy), terra::extract(r, xy, ID = FALSE, method = "bilinear"))
            a <- analogs::analog_impact(x = q, pool = pool$aim_index, y = pool$values[, cols, drop = FALSE],
                                        env = k$env, geog = k$geog,
                                        stat = c("weighted_mean", weight_stats), coord_type = "projected")
            a <- name_outcomes(a[, setdiff(names(a), c("index", "x", "y")), drop = FALSE], cols, weight_stats)
            c(value = if (is.null(ratio)) a$response else ts_combine_ratio(a$occ, a$num, a$denom, s$fia_var),
              support = a$sum_weights, closeness = a$effective_weight)
      })
      structure(c(hst = out$hst[["value"]], fut = out$fut[["value"]]),
                evidence = c(support_hst = out$hst[["support"]], support_fut = out$fut[["support"]],
                             closeness_hst = out$hst[["closeness"]], closeness_fut = out$fut[["closeness"]]))
}
value_of <- function(x) { attributes(x) <- list(names = names(x)); x }

test_that("the current settings give one estimate, in every row, as analog_impact() would", {
      d <- ts_site_sensitivity(pipo(), data, site)
      now <- d[d$current, ]
      expect_setequal(now$row, c("gcm", "ssp", "theta_geog", "theta_clim", "clim_vars", "sample"))
      for (col in c("hst", "fut", "delta", "support_hst", "support_fut"))
            expect_identical(length(unique(now[[col]])), 1L, info = col)
      ref <- reference(pipo())
      expect_equal(c(hst = now$hst[1], fut = now$fut[1]), value_of(ref), tolerance = 1e-10)
      # Evidence as the map's support and closeness layers measure it
      expect_equal(unlist(now[1, names(attr(ref, "evidence"))]), attr(ref, "evidence"), tolerance = 1e-10)

      # And away from the current settings
      s <- pipo()
      wide <- d[d$row == "theta_clim" & d$alternative == as.character(s$theta_clim * 2), ]
      expect_equal(c(hst = wide$hst, fut = wide$fut), value_of(reference(s, theta_clim = s$theta_clim * 2)),
                   tolerance = 1e-10)

      # A swapped climate variable is a run on that set of variables
      swap <- d[d$row == "clim_vars" & d$alternative == "aet, tmean", ]
      expect_identical(nrow(swap), 1L)
      expect_equal(c(hst = swap$hst, fut = swap$fut),
                   value_of(reference(pipo(clim_vars = c("aet", "tmean")))), tolerance = 1e-10)
      expect_identical(sum(d$row == "clim_vars"), 1L + 2L * 4L)   # current, and 2 x 4 swaps

      # Climate scenarios don't change the historical estimate
      scen <- d[d$row %in% c("gcm", "ssp"), ]
      expect_true(all(scen$hst == now$hst[1]))
      expect_identical(nrow(scen), 5L + 3L)
      expect_true(all(d$hst_lower[d$row == "sample"] <= now$hst[1],
                      d$hst_upper[d$row == "sample"] >= now$hst[1]))
})

test_that("a pooled ratio runs, combined as the map combines it", {
      s <- pipo(stage = "seedling", metric = "regeneration fraction")
      d <- ts_site_sensitivity(s, data, site)
      now <- d[d$row == "sample", ]
      expect_equal(c(hst = now$hst, fut = now$fut), value_of(reference(s)), tolerance = 1e-10)
      expect_false(is.na(now$occ_hst))
      expect_false(is.na(now$hst_lower))
})

test_that("the bootstrap is repeatable and leaves the session's random numbers alone", {
      # The same result whatever state the session's random numbers are in
      set.seed(1)
      a <- ts_site_sensitivity(ts_settings(), data, site)
      set.seed(2)
      b <- ts_site_sensitivity(ts_settings(), data, site)
      expect_identical(a, b)

      # The bootstrap draws from its own seed. (analogs::build_analog_index()
      # does advance the session's random numbers, as the map's pool does too.)
      w <- c(0, 1, 0.5, 0.2, 0)
      est <- function(ws, at) list(value = colSums(ws) / colSums(ws))
      set.seed(1)
      before <- .Random.seed
      ts_site_bootstrap(w, w, est, 50, seed = "x")
      expect_identical(.Random.seed, before)
})

test_that("the estimate doesn't depend on where the bounding box ends", {
      # A specific model: a storyline would pick its model for each box
      a <- ts_site_sensitivity(ts_settings(scenario_model = "mri-esm2-0"), data, site)
      # The site just inside the box's west edge, so the box's own analog pool
      # would miss plots to the west
      b <- ts_site_sensitivity(ts_settings(scenario_model = "mri-esm2-0",
                                           bbox_x = c(site[1] - 0.05, -105), bbox_y = c(33, site[2] + 0.05)),
                               data, site)
      # The seed hashes the settings, so only the point estimates must match
      cols <- c("row", "alternative", "current", "hst", "fut", "delta", "n_plots")
      expect_identical(a[cols], b[cols])
})

test_that("the ratio helper gives the same values on rasters as on vectors, blank without climate", {
      r <- terra::rast(nrows = 2, ncols = 3, nlyrs = 3, names = c("occ", "num", "denom"),
                       vals = c(0.5, 0.2, 0, NA, 1, 0.8,
                                0.3, 0.1, 0, NA, 0.4, 0,
                                0.6, 0.2, 0, NA, 0.5, 0.4))
      v <- terra::values(r, dataframe = TRUE)
      for (var in c("seedling fraction", "proportion basal area")) {
            from_raster <- terra::values(ts_combine_ratio(r$occ, r$num, r$denom, var), mat = FALSE)
            expect_identical(from_raster, ts_combine_ratio(v$occ, v$num, v$denom, var), info = var)
      }
      # No occupied plot in range (cell 3): blank for a fraction, 0 for a
      # proportion. No climate (cell 4): blank for both.
      expect_identical(ts_combine_ratio(v$occ, v$num, v$denom, "seedling fraction")[3:4], c(NA_real_, NA))
      expect_identical(ts_combine_ratio(v$occ, v$num, v$denom, "proportion basal area")[3:4], c(0, NA))
})
