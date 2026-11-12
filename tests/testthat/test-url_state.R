# The URL codec (helpers/url_state.R), without running the app

load_app_globals()
ctx <- url_context()
d <- url_defaults(ctx)

roundtrip <- function(s) {
      q <- url_encode(s, ctx)
      r <- url_decode(q, ctx)
      list(query = q, again = url_encode(r$state, ctx), problems = r$problems)
}

test_that("the default link is fully specified", {
      expect_identical(
            url_encode(d, ctx),
            "?v=1&e=sp&es=single&sp=Pinus_edulis&st=any&m=occ&t=hist&g=10&bx=-116,-100&by=29,43&cv=aet,cwd&ssp=ssp585&era=2041-2070&gcm=gfdl-esm4&tg=30&tc=0.1")
})

test_that("links round-trip without problems", {
      states <- list(
            species = modifyList(d, list(
                  species = "Pinus ponderosa", stage = "seedling", metric = "abundance",
                  fia_var = "number of seedlings", time = "delta", polygon_unit = "county:08069",
                  site = c(-105.5, 40.6), clim_vars = c("tmean", "ppt", "cwd"),
                  scenario_model = unname(ctx$models[1:2]), theta_geog = 45, theta_clim = 0.25,
                  bbox_x = c(-110.5, -101))),
            community = modifyList(d, list(
                  entity_scope = "multi", comm_measure = "dominant entity", comm_top_n = 7,
                  support_floor = 1.5, stat = "support")),
            forest_type = modifyList(d, list(
                  track = "Forest type", forest_type = ctx$forest_types[5], time = "future",
                  vector = "none", analog_direction = "Reverse analogs")),
            connectivity = modifyList(d, list(
                  track = "All trees", metric = "species richness", stage = "adult",
                  fia_var = resolve_var("All trees", "species richness", "adult"),
                  cs_enable = TRUE, cs_sources = "FIA plots", cs_loss_layer = "mean suitability",
                  cs_loss_scale = -2.5, stat = "connectivity", cs_stat = "voltage",
                  time = "delta", vector = "current field")))
      for (nm in names(states)) {
            r <- roundtrip(states[[nm]])
            expect_identical(r$again, r$query, label = nm)
            expect_length(r$problems, 0)
      }
})

test_that("settings that don't apply are left out", {
      q <- url_encode(modifyList(d, list(track = "Forest type")), ctx)
      expect_false(grepl("&sp=", q, fixed = TRUE))
      expect_false(grepl("&st=", q, fixed = TRUE))
      expect_true(grepl("&ft=", q, fixed = TRUE))
      # Connectivity settings appear only when it is on
      expect_false(grepl("&cs=|&csrc=", url_encode(d, ctx)))
})

test_that("a short link takes defaults for what it leaves out", {
      r <- url_decode("?sp=Pinus_ponderosa&st=adult&m=abund", ctx)
      expect_length(r$problems, 0)
      expect_identical(r$state$species, "Pinus ponderosa")
      expect_identical(r$state$fia_var, "number of trees")
      expect_identical(r$state$theta_geog, d$theta_geog)
})

test_that("bad values fall back to defaults and are reported", {
      r <- url_decode(paste0("?v=2&e=sp&sp=Pinus_notreal&st=any&m=abund&g=0.5&bx=-100,-116",
                             "&tg=500&poly=wilderness:654&site=10,10&cv=cwd,cwd&zz=1"), ctx)
      s <- r$state
      expect_identical(s[c("species", "metric", "clim_agg", "bbox_x", "theta_geog", "clim_vars")],
                       d[c("species", "metric", "clim_agg", "bbox_x", "theta_geog", "clim_vars")])
      expect_identical(s$polygon_unit, "")
      expect_null(s$site)
      expect_length(r$problems, 10)
      expect_true(any(grepl("newer version", r$problems)))
      expect_true(any(grepl("doesn't recognize: zz", r$problems)))
      expect_true(any(grepl("^Variable: 'abund' isn't available with the link's other settings", r$problems)))
      expect_true(any(grepl("wilderness:654", r$problems)))
})

test_that("dependent settings are checked against the settings they depend on", {
      # Categorical measures aren't offered as change
      r <- url_decode("?es=multi&cm=domid&t=chg", ctx)
      expect_identical(r$state$comm_measure, "richness")
      expect_length(r$problems, 1)
      # Connectivity is single-entity only
      r <- url_decode("?es=multi&cs=1", ctx)
      expect_false(r$state$cs_enable)
      expect_length(r$problems, 1)
      # The connectivity raster needs connectivity on
      r <- url_decode("?r=conn", ctx)
      expect_identical(r$state$stat, "FIA variable")
      # ...and is shown as change
      r <- url_decode("?cs=1&r=conn", ctx)
      expect_identical(r$state$time, "delta")
})

test_that("multi-selects are written in menu order", {
      a <- url_encode(modifyList(d, list(clim_vars = c("cwd", "tmean"))), ctx)
      b <- url_encode(modifyList(d, list(clim_vars = c("tmean", "cwd"))), ctx)
      expect_identical(a, b)
})
