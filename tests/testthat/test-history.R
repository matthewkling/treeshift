# Browser history (server/url_state.R): each model update is a Back step, and
# Back and Forward restore that step's settings

skip_on_cran()
load_app_globals()

test_that("Back and Forward step through model updates", {
      a <- open_session(app_dir)
      on.exit(a$stop(), add = TRUE)
      steps <- function() a$get_js("history.length")
      start <- steps()

      set_in_order(a, species = "Pinus ponderosa"); apply_settings(a)
      set_in_order(a, theta_geog = 45); apply_settings(a)
      expect_identical(steps(), start + 2L)
      # Display settings replace the current step rather than adding one
      set_in_order(a, time = "future")
      expect_identical(steps(), start + 2L)

      expect_true(history_nav(a, "back"))
      expect_identical(page_input(a, "species"), "Pinus ponderosa")
      expect_identical(page_input(a, "theta_geog"), 30L)

      # Back to the link the session opened with
      expect_true(history_nav(a, "back"))
      expect_identical(page_input(a, "species"), "Pinus edulis")

      expect_true(history_nav(a, "forward"))
      expect_identical(page_input(a, "species"), "Pinus ponderosa")
      expect_identical(notices(a), "")
      expect_identical(a$get_js("document.getElementById('update_model').disabled"), TRUE)
})

test_that("moving between About sections doesn't reload", {
      a <- open_session(app_dir)
      on.exit(a$stop(), add = TRUE)
      a$run_js("document.querySelector(\".navbar .nav-link[data-value='about']\").click()")
      settle(a)
      a$run_js("document.querySelector(\"a[href='#about-data']\").click()")
      expect_false(history_nav(a, "back"))
      expect_identical(a$get_js("document.querySelector('.navbar .nav-link.active').textContent.trim()"), "About")
})
