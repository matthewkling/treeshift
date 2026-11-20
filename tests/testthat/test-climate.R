# Climate storylines (R/climate.R): which model each one picks

withr::with_dir(app_dir, source("scripts/load_core.R"))

pick <- function(data, storyline, bbox_x = c(-110, -106), bbox_y = c(38, 42))
      ts_storyline_model(data, storyline, "ssp370", "2041-2070", bbox_x, bbox_y)

test_that("storylines pick the least, median and most drying models", {
      d <- fake_cube(list(a = 3, b = 1, c = 5, d = 2, e = 4), withr::local_tempdir())
      expect_identical(vapply(storyline_choices, function(s) pick(d, s), ""),
                       c("Least drying" = "b", "Median drying" = "a", "Most drying" = "c"))
      expect_identical(names(ts_storyline_ranking(d, "ssp370", "2041-2070", c(-110, -106), c(38, 42))),
                       c("b", "d", "a", "e", "c"))
      expect_error(pick(d, "wettest"), "Unknown climate storyline")
})

test_that("tied models keep model order", {
      d <- fake_cube(list(a = 2, b = 1, c = 1, d = 2), withr::local_tempdir())
      expect_identical(pick(d, "least"), "b")
      expect_identical(pick(d, "median"), "c")   # rank ceiling(4 / 2) = 2
      expect_identical(pick(d, "most"), "d")
})

test_that("the ranking is over the bounding box, which may be smaller than a cell", {
      # Model a dries most in the northwest cell only, model b everywhere else
      d <- fake_cube(list(a = c(9, rep(0, 15)), b = 1), withr::local_tempdir())
      expect_identical(pick(d, "most"), "b")
      expect_identical(pick(d, "most", c(-109.8, -109.6), c(41.6, 41.8)), "a")
      expect_error(pick(d, "most", c(-100, -99), c(38, 39)), "No climate data|extents do not overlap")
})

test_that("a missing layer is an error, not a smaller ranking", {
      d <- fake_cube(list(a = 1, b = 2), withr::local_tempdir())
      d$scenarios <- rbind(d$scenarios, tibble::tibble(ssp = "ssp370", period = "2041-2070", model = "c"))
      expect_error(pick(d, "median"), "lacks ssp370_2041-2070_c_cwd")
})

test_that("resolving fills in a storyline's model and leaves a chosen model alone", {
      d <- fake_cube(list(a = 3, b = 1, c = 5), withr::local_tempdir())
      s <- ts_settings(bbox_x = c(-110, -106), bbox_y = c(38, 42))
      expect_null(s$scenario_model)
      expect_error(ts_scenario_model(s), "ts_resolve_scenario")
      expect_identical(ts_resolve_scenario(s, d)$scenario_model, "a")
      expect_identical(ts_resolve_scenario(ts_settings(scenario_storyline = "most",
                                                       bbox_x = c(-110, -106), bbox_y = c(38, 42)), d)$scenario_model, "c")
      chosen <- ts_settings(scenario_model = "b")
      expect_identical(ts_resolve_scenario(chosen, d), chosen)
})
