# Model settings and the Update map button (helpers/model_settings.R,
# server/model_settings.R)

load_app_globals()

test_that("server code reads model settings only through their applied values", {
      # Only the controls' own menu logic and the button's code read the
      # controls directly. Anywhere else, a direct read would mix pending
      # settings into results computed from the applied ones.
      pattern <- paste0("input\\$(", paste(model_setting_ids, collapse = "|"), ")\\b")
      files <- setdiff(list.files(file.path(app_dir, "server"), "[.]R$", full.names = TRUE),
                       file.path(app_dir, "server", c("controls.R", "model_settings.R")))
      hits <- unlist(lapply(files, function(f) {
            code <- readLines(f, warn = FALSE)
            code <- code[!grepl("^\\s*#", code)]
            hit <- grep(pattern, code, value = TRUE, perl = TRUE)
            if (length(hit)) paste0(basename(f), ": ", trimws(hit))
      }))
      expect_length(hits, 0)
})

test_that("settings compare by value, and multi-selects as sets", {
      expect_true(same_setting(30L, 30))
      expect_true(same_setting(c("cwd", "aet"), c("aet", "cwd")))
      expect_false(same_setting(c(-116, -100), c(-116, -101)))
      expect_true(same_setting(NULL, character(0)))
      expect_false(same_setting("aet", c("aet", "cwd")))
})

test_that("only changes to settings in use are pending", {
      d <- url_defaults(url_ctx)
      now <- modifyList(d, list(forest_type = "221", theta_geog = 45))
      expect_identical(pending_settings(now, d), "theta_geog")
      now$track <- "Forest type"
      expect_setequal(pending_settings(now, d), c("track", "forest_type", "theta_geog"))
      # Connectivity settings count only while connectivity is on
      now <- modifyList(d, list(cs_loss_scale = -2))
      expect_length(pending_settings(now, d), 0)
      now$cs_enable <- TRUE
      expect_identical(pending_settings(now, d), "cs_loss_scale")
})

test_that("the pending line names each change once", {
      expect_identical(pending_phrase(c("metric", "fia_var", "theta_geog")),
                       "Variable and geographic bandwidth changed.")
})

test_that("changes wait for the Update button", {
      skip_on_cran()
      a <- open_session(app_dir)
      on.exit(a$stop(), add = TRUE)
      button <- function() a$get_js("(() => { const b = document.getElementById('update_model'); return b.textContent.trim() + (b.disabled ? ' [disabled]' : ''); })()")
      washed <- function() a$get_js("document.querySelector('.ts-map-stale') !== null")

      expect_identical(button(), "Map is up to date [disabled]")
      expect_false(washed())
      before <- list(query = address_bar(a), raster = raster_md5(a))

      set_in_order(a, species = "Pinus ponderosa", theta_geog = 45)
      expect_identical(button(), "Update map")
      expect_identical(a$get_js("document.querySelector('.ts-pending-note').textContent"),
                       "Species and geographic bandwidth changed.")
      expect_true(washed())
      # The map, the link, and the downloads still describe the applied settings
      expect_identical(address_bar(a), before$query)
      expect_identical(raster_md5(a), before$raster)

      apply_settings(a)
      expect_identical(button(), "Map is up to date [disabled]")
      expect_false(washed())
      expect_match(address_bar(a), "&sp=Pinus_ponderosa&.*&tg=45")
      expect_false(identical(raster_md5(a), before$raster))

      # Display settings still apply right away
      set_in_order(a, time = "future")
      expect_identical(button(), "Map is up to date [disabled]")
      expect_match(address_bar(a), "&t=fut&", fixed = TRUE)
})
