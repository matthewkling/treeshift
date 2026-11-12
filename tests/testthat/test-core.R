# The modeling core (R/), run outside Shiny as a script would run it

withr::with_dir(app_dir, source("scripts/load_core.R"))
data <- withr::with_dir(app_dir, ts_data("data"))

# A small domain, so each run takes seconds
small <- function(...) ts_settings(bbox_x = c(-112, -108), bbox_y = c(38, 42), ...)

test_that("ts_settings() resolves the variable and rejects what it can't run", {
      expect_equal(ts_settings()$fia_var, "presence probability")
      expect_equal(ts_settings(stage = "seedling", metric = "regeneration fraction")$fia_var,
                   "seedling fraction")
      expect_error(ts_settings(theta = 1), "Unknown settings")
      expect_error(ts_settings(metric = "abundance"), "No variable")   # needs a stage
      expect_error(ts_settings(track = "Forest type"), "forest_type")
      expect_error(ts_settings(clim_vars = "elevation"), "Unknown climate")
      expect_error(ts_model(ts_settings(scenario_model = "nope"), data), "No climate data")
})

test_that("community richness equals the All-trees richness variable", {
      # AIM is a weighted mean with entity-independent weights, so summing
      # per-entity occurrence surfaces must equal generalizing the per-plot
      # species count. Two independent routes through the pipeline -- a check
      # on the per-entity response matrix's row alignment in particular.
      comm <- ts_model(small(entity_scope = "multi"), data)
      all <- ts_model(small(track = "All trees", metric = "species richness"), data)
      a <- terra::values(ts_layer(comm, measure = "richness"), mat = FALSE)
      b <- terra::values(ts_layer(all), mat = FALSE)
      expect_identical(is.na(a), is.na(b))
      expect_lt(max(abs(a - b), na.rm = TRUE), 1e-8)
      expect_gt(sum(!is.na(a)), 0)
})

test_that("a regeneration fraction is the pooled ratio, blank where the species is absent", {
      res <- ts_model(small(species = "Pinus ponderosa", stage = "seedling",
                            metric = "regeneration fraction"), data)
      s <- res$surfaces$bsl
      ratio <- terra::values(s$num / s$denom, mat = FALSE)
      ratio[!is.finite(ratio)] <- NA
      expect_identical(terra::values(s$response, mat = FALSE), ratio)
      expect_true(all(is.na(terra::values(s$response, mat = FALSE)[
            terra::values(s$denom, mat = FALSE) == 0])))
      expect_false("spread" %in% names(s))   # no spread for a pooled ratio
})

test_that("change is future minus historical, and evidence uses the harmonic mean", {
      res <- ts_model(small(), data)
      expect_identical(terra::values(ts_layer(res, "delta"), mat = FALSE),
                       terra::values(res$surfaces$fut$response - res$surfaces$bsl$response,
                                     mat = FALSE))
      sw <- function(tf) terra::values(res$surfaces[[tf]]$sum_weights, mat = FALSE)
      expect_equal(terra::values(ts_evidence_layer(res$surfaces, "sum_weights", "delta"),
                                 mat = FALSE),
                   2 / (1 / sw("fut") + 1 / sw("bsl")))
})

test_that("the per-entity response matrix is row-aligned with the plots", {
      pb <- ts_pool_bbox(c(-112, -108), c(38, 42), 30)
      plots <- ts_plots(data, pb, "Pinus ponderosa")
      y <- ts_multi_response(data, plots, pb, "Species", "presence probability")
      expect_equal(nrow(y), nrow(plots))
      expect_identical(y[["Pinus ponderosa"]],
                       as.numeric(ts_response(plots, "Species", "presence probability")))
})
