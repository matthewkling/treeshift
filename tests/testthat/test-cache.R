# The result cache (R/cache.R): what its keys depend on, and that a cached result
# is the computed one

withr::with_dir(app_dir, source("scripts/load_core.R"))

# Stand-in input files: the keys hash the files, so they needn't be real data
fake_data <- function() {
      dir <- withr::local_tempdir(.local_envir = parent.frame())
      dir.create(file.path(dir, "clim_future"))
      files <- c("fia.csv", "forest_type_codes.csv", "clim_scales.rds", "clim_historical.tif",
                 "wind_quantized.tif", "clim_future/ssp585_2041-2070_gfdl-esm4.tif",
                 "clim_future/ssp585_2041-2070_mri-esm2-0.tif")
      for (f in files) writeLines(f, file.path(dir, f))
      list(dir = dir, clim_historical = file.path(dir, "clim_historical.tif"),
           clim_future_dir = file.path(dir, "clim_future"))
}
key <- function(settings, data) withr::with_dir(app_dir, ts_surfaces_key(settings, data))

test_that("the surfaces key changes with every setting it is computed from", {
      d <- fake_data()
      base <- key(ts_settings(), d)
      changed <- list(species = "Pinus ponderosa", stage = "seedling", clim_vars = c("aet", "ppt"),
                      scenario_model = "mri-esm2-0", bbox_x = c(-116, -101), clim_agg = 20,
                      theta_geog = 40, theta_clim = 0.2, track = "All trees")
      for (nm in names(changed))
            expect_false(identical(key(do.call(ts_settings, changed[nm]), d), base), info = nm)
})

test_that("the surfaces key ignores what the run doesn't read", {
      d <- fake_data()
      # The species is read on the Species track only
      expect_identical(key(ts_settings(track = "All trees", species = "Pinus ponderosa"), d),
                       key(ts_settings(track = "All trees"), d))
      # Nor do display settings or how a number arrived (slider double vs integer)
      s <- ts_settings()
      s2 <- s; s2$clim_agg <- 10L; s2$entity_scope <- "single"
      expect_identical(key(s2, d), key(s, d))
})

test_that("the surfaces key changes when an input file's content changes", {
      d <- fake_data()
      base <- key(ts_settings(), d)
      writeLines("revised", file.path(d$dir, "fia.csv"))
      expect_false(identical(key(ts_settings(), d), base))
      file.remove(file.path(d$dir, "fia.csv"))
      expect_null(key(ts_settings(), d))   # a missing input disables caching
})

test_that("the connectivity key includes the focal site only when an electrode uses it", {
      d <- fake_data()
      cs <- function(sources = "grid cells")
            c(ts_settings(), list(cs_sources = sources, cs_targets = "grid cells",
                                  cs_conductance_layer = "suitability", cs_loss_layer = "none",
                                  cs_conductance_scale = 0, cs_loss_scale = -1))
      spec <- function(s, site) ts_connectivity_spec(s, site)
      expect_identical(spec(cs(), c(-105, 40)), spec(cs(), c(-106, 41)))
      expect_identical(spec(cs(), c(-105, 40)), spec(cs(), NULL))
      site_spec <- spec(cs("selected site"), c(-105, 40))
      expect_false(identical(site_spec, spec(cs("selected site"), c(-106, 41))))
      expect_false(identical(site_spec, spec(cs("selected site"), NULL)))
      expect_false(identical(spec(cs("selected site"), NULL), spec(cs("FIA plots"), NULL)))
})

test_that("a cached result is the computed one, and a bad entry is recomputed", {
      skip_on_cran()
      data <- withr::with_dir(app_dir, ts_data("data"))
      s <- ts_settings(bbox_x = c(-112, -108), bbox_y = c(38, 42))
      cache <- ts_cache(withr::local_tempdir())
      k <- key(s, data)
      runs <- 0
      compute <- function() { runs <<- runs + 1; ts_model(s, data)$surfaces }

      live <- ts_cached(cache, k, compute)
      cached <- ts_cached(cache, k, compute)
      expect_equal(runs, 1)
      for (t in c("bsl", "fut")) {
            expect_identical(names(cached[[t]]), names(live[[t]]))
            expect_identical(terra::values(cached[[t]]), terra::values(live[[t]]))
            expect_identical(as.vector(terra::ext(cached[[t]])), as.vector(terra::ext(live[[t]])))
            expect_identical(terra::crs(cached[[t]]), terra::crs(live[[t]]))
      }
      # Down to the file a user downloads
      f <- withr::local_tempfile(fileext = c(".tif", ".tif"))
      terra::writeRaster(live$bsl, f[1]); terra::writeRaster(cached$bsl, f[2])
      expect_identical(unname(tools::md5sum(f[1])), unname(tools::md5sum(f[2])))

      # A corrupted entry is treated as missing
      writeLines("not an rds file", file.path(cache$info()$dir, paste0(k, ".rds")))
      again <- ts_cached(cache, k, compute)
      expect_equal(runs, 2)
      expect_identical(terra::values(again$bsl), terra::values(live$bsl))

      # No cache, or no key, just computes
      ts_cached(NULL, k, compute)
      ts_cached(cache, NULL, compute)
      expect_equal(runs, 4)
})
