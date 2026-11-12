# The evidence map layers: analog support, closeness and spread
# (server/aim.R, server/map_layers.R).

load_app_globals()

test_that("aim_surfaces() gives support, closeness and spread that match hand-computed weights", {
      set.seed(1)
      n <- 300
      xy <- data.frame(x = runif(n, -111, -108), y = runif(n, 39, 42))
      clim <- data.frame(c1 = rnorm(n), c2 = rnorm(n))
      aim_pts <- data.frame(aim_xy(xy), clim)
      pool <- list(aim_pts = aim_pts,
                   aim_index = analogs::build_analog_index(aim_pts, coord_type = "projected",
                                                           geog_res_adj = 4))
      y <- data.frame(response = rgamma(n, 2))

      focal <- terra::rast(xmin = -110, xmax = -109, ymin = 40, ymax = 41, res = 0.25, nlyrs = 2,
                           names = c("c1", "c2"), crs = "EPSG:4326")
      terra::values(focal) <- matrix(rnorm(2 * terra::ncell(focal)), ncol = 2)

      theta_clim <- 0.8
      theta_geog <- 60
      surf <- aim_surfaces(list(bsl = focal), pool, y, theta_clim, theta_geog,
                           stat = c("weighted_mean", "weighted_sd", "sum_weights",
                                    "effective_weight"))$bsl
      expect_setequal(names(surf), c("weighted_mean", "weighted_sd", "sum_weights",
                                     "effective_weight"))

      fd <- as.data.frame(focal, xy = TRUE, cells = TRUE)
      for (i in c(1, 7, 16)) {
            fxy <- aim_xy(fd[i, c("x", "y")])
            dg <- sqrt((aim_pts$x - fxy[1, "x"])^2 + (aim_pts$y - fxy[1, "y"])^2)
            dc <- sqrt((clim$c1 - fd$c1[i])^2 + (clim$c2 - fd$c2[i])^2)
            w <- exp(-dg^2 / (2 * theta_geog^2)) * exp(-dc^2 / (2 * theta_clim^2))
            w[dg > 3 * theta_geog | dc > 3 * theta_clim] <- 0
            m <- sum(w * y$response) / sum(w)
            cell <- fd$cell[i]
            expect_equal(surf$sum_weights[cell][[1]], sum(w), tolerance = 1e-6)
            expect_equal(surf$effective_weight[cell][[1]], sum(w^2) / sum(w), tolerance = 1e-6)
            expect_equal(surf$weighted_sd[cell][[1]],
                         sqrt(sum(w * (y$response - m)^2) / sum(w)), tolerance = 1e-6)
      }
})


skip_on_cran()
app <- shinytest2::AppDriver$new(app_dir, width = 1400, height = 900,
                                 load_timeout = 600000, timeout = 600000)
withr::defer(app$stop(), teardown_env())
settle(app)

# The note under the raster menu. Outputs in a collapsed panel aren't
# rendered, so this opens the Map layers panel first.
stat_note <- function(app) {
      app$run_js("var b = Array.from(document.querySelectorAll('.accordion-button')).find(b => b.textContent.trim() === 'Map layers'); if (b.classList.contains('collapsed')) b.click();")
      settle(app)
      app$get_js("document.getElementById('stat_note').innerText")
}

test_that("spread for a presence variable is the SD of a 0/1 outcome", {
      set_in_order(app, bbox_x = c(-112, -104), bbox_y = c(37, 43), stat = "evidence", evidence_stat = "spread")
      apply_settings(app)
      expect_identical(stat_note(app), "")
      spread <- raster_values(app)
      expect_gt(sum(is.finite(spread)), 100)
      expect_gte(min(spread, na.rm = TRUE), 0)
      expect_lte(max(spread, na.rm = TRUE), 0.5 + 1e-9)
})

test_that("spread is withheld for a ratio variable, with a note saying why", {
      set_in_order(app, metric = "proportion basal area")
      apply_settings(app)
      expect_match(stat_note(app), "isn't shown for ratio variables", fixed = TRUE)
      set_in_order(app, metric = "occurrence")
      apply_settings(app)
})

test_that("the evidence layers draw in their valid ranges for a continuous variable", {
      set_in_order(app, species = "Pinus ponderosa", metric = "basal area")
      apply_settings(app)
      expect_identical(app$get_value(input = "fia_var"), "total basal area")
      expect_identical(stat_note(app), "")

      spread <- raster_values(app)
      expect_gt(sum(is.finite(spread)), 100)
      expect_gte(min(spread, na.rm = TRUE), 0)

      set_in_order(app, stat = "evidence", evidence_stat = "closeness")
      closeness <- raster_values(app)
      expect_gte(min(closeness, na.rm = TRUE), 0)
      expect_lte(max(closeness, na.rm = TRUE), 1)

      set_in_order(app, stat = "evidence", evidence_stat = "support")
      support <- raster_values(app)
      # drawn in fixed classes, labelled in plain units
      legend <- app$get_js("document.querySelector('#map .info.legend').innerText")
      expect_match(legend, "under 1", fixed = TRUE)
      expect_match(legend, "100 or more", fixed = TRUE)
      # support = effective number of plots x closeness, and that number is at least 1
      both <- is.finite(support) & is.finite(closeness)
      expect_true(all(support[both] >= closeness[both] - 1e-9))

      # Change shows the support of the weaker end; spread has no change version
      set_in_order(app, time = "delta")
      expect_gt(sum(is.finite(raster_values(app))), 100)
      set_in_order(app, stat = "evidence", evidence_stat = "spread")
      expect_match(stat_note(app), "not for change", fixed = TRUE)
})


test_that("transformed legends are labelled in the variable's units", {
      set_in_order(app, time = "historical", stat = "FIA variable")
      raw <- raster_values(app)          # downloads are untransformed
      top <- max(raw, na.rm = TRUE)
      legend_numbers <- function() {
            txt <- app$get_js("document.querySelector('#map .info.legend').innerText")
            as.numeric(gsub(",", "", regmatches(txt, gregexpr("-?[0-9][0-9,]*\\.?[0-9]*", txt))[[1]]))
      }
      for (tr in c("log10", "sqrt", "uniform")) {
            set_in_order(app, rast_transform = tr)
            n <- legend_numbers()
            n <- n[n != 2]                 # "ft2" in the units line
            expect_true(all(n >= 0 & n <= top * 1.01), info = tr)
            # Transformed units would top out near log10(top) or sqrt(top)
            expect_gt(max(n), sqrt(top) + 1)
      }
      set_in_order(app, rast_transform = "identity")
})
