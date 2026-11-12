# Links in the running app: settings made in one session, carried by the
# address bar's link into a fresh session, give the same inputs and the same
# raster.

skip_on_cran()
load_app_globals()

app <- shinytest2::AppDriver$new(app_dir, width = 1400, height = 900,
                                 load_timeout = 600000, timeout = 600000)
withr::defer(app$stop(), teardown_env())
settle(app)
base <- sub("\\?.*$", "", app$get_url())

# Make settings in a fresh session, then open its link in another
check_link <- function(make_settings) {
      a <- open_session(base)
      on.exit(a$stop(), add = TRUE)
      make_settings(a)
      apply_settings(a)
      query <- address_bar(a)
      b <- open_session(paste0(base, query))
      on.exit(b$stop(), add = TRUE)

      expect_identical(address_bar(b), query)
      ids <- inputs_set_by(query)
      expect_identical(b$get_values(input = ids)$input, a$get_values(input = ids)$input)
      expect_identical(notices(b), "")
      expect_identical(raster_md5(b), raster_md5(a))
      query
}

test_that("a species link with a polygon and site reopens the same results", {
      q <- check_link(function(a) {
            set_in_order(a, species = "Pinus ponderosa", stage = "seedling", metric = "abundance",
                         time = "delta", theta_geog = 45, clim_vars = c("tmean", "ppt", "cwd"),
                         polygon_collection = "Counties", polygon_unit = "county:08069",
                         analog_direction = "Reverse analogs")
            click_map(a, -105.5, 40.6)
      })
      expect_match(q, "poly=county:08069", fixed = TRUE)
      expect_match(q, "site=-105.5,40.6", fixed = TRUE)
})

test_that("a forest type link with map layers reopens the same results", {
      check_link(function(a) {
            set_in_order(a, track = "Forest type", forest_type = "221", time = "future",
                         bbox_x = c(-112, -104), stat = "closeness", vector = "none")
      })
})

test_that("a community link reopens the same results", {
      check_link(function(a) {
            set_in_order(a, bbox_x = c(-110, -104), bbox_y = c(35, 40), entity_scope = "multi",
                         comm_measure = "dominant entity", comm_top_n = 6)
      })
})

test_that("a link with connectivity settings reopens them", {
      q <- check_link(function(a) {
            set_in_order(a, track = "All trees", stage = "adult", cs_enable = TRUE,
                         cs_sources = "FIA plots", cs_loss_layer = "mean suitability",
                         cs_loss_scale = -2)
      })
      expect_match(q, "cs=1&csrc=plots&ctgt=cells&ccl=suit&cll=suit&ccs=0&cls=-2", fixed = TRUE)
})

test_that("a bad link opens on defaults and says what it couldn't use", {
      b <- open_session(paste0(base, "?e=sp&sp=Pinus_notreal&m=abund&poly=wilderness:654&tg=500"))
      on.exit(b$stop(), add = TRUE)
      n <- notices(b)
      expect_match(n, "Species: 'Pinus_notreal' isn't a recognized value", fixed = TRUE)
      expect_match(n, "Focal polygon 'wilderness:654' isn't in this version's boundary library", fixed = TRUE)
      expect_match(n, "Geographic bandwidth: '500'", fixed = TRUE)
      expect_identical(address_bar(b), url_encode(url_defaults(url_ctx), url_ctx))
})
