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
            "?v=1&e=sp&es=single&sp=Pinus_edulis&st=any&m=occ&t=hist&g=10&bx=-116,-100&by=29,43&cv=aet,cwd&ssp=ssp370&era=2041-2070&sto=median&gcm=mpi-esm1-2-hr&tg=30&tc=0.1&crch=100&clu=50")
})

test_that("links round-trip without problems", {
      states <- list(
            species = modifyList(d, list(
                  species = "Pinus ponderosa", stage = "seedling", metric = "abundance",
                  fia_var = "number of seedlings", time = "delta", polygon_unit = "county:08069",
                  site = c(-105.5, 40.6), clim_vars = c("tmean", "ppt", "cwd"),
                  scenario_model = unname(ctx$models[2]), theta_geog = 45, theta_clim = 0.25,
                  bbox_x = c(-110.5, -101))),
            community = modifyList(d, list(
                  entity_scope = "multi", comm_measure = "dominant entity", comm_top_n = 7,
                  support_floor = 1.5, stat = "evidence", evidence_stat = "closeness")),
            forest_type = modifyList(d, list(
                  track = "Forest type", forest_type = ctx$forest_types[5], time = "future",
                  vector = "none", analog_direction = "Reverse analogs")),
            connectivity = modifyList(d, list(
                  track = "All trees", metric = "basal area", stage = "any",
                  fia_var = resolve_var("All trees", "basal area", "any"),
                  cs_reach = "Inf", cs_landuse = 75, stat = "connectivity", time = "future",
                  vector = "current field")),
            connected_habitat = modifyList(d, list(
                  stat = "connectivity", conn_stat = "colonization", time = "delta")),
            habitat_within_reach = modifyList(d, list(
                  stat = "connectivity", conn_stat = "reach", time = "future")),
            habitat_within_reach_change = modifyList(d, list(
                  stat = "connectivity", conn_stat = "reach_tracking", time = "delta")))
      # A specific model rather than a storyline
      states$species["scenario_storyline"] <- list(NULL)
      for (nm in names(states)) {
            r <- roundtrip(states[[nm]])
            expect_identical(r$again, r$query, label = nm)
            expect_length(r$problems, 0)
      }
})

test_that("a connectivity link without a measure opens on corridors", {
      r <- url_decode(sub("&r=conn&rcn=corr", "&r=conn", url_encode(modifyList(d, list(stat = "connectivity")), ctx)), ctx)
      expect_identical(r$state$conn_stat, "corridors")
      expect_length(r$problems, 0)
})

test_that("a connectivity measure the timeframe doesn't offer is reported", {
      q <- url_encode(modifyList(d, list(stat = "connectivity", conn_stat = "tracking", time = "delta")), ctx)
      r <- url_decode(sub("t=chg", "t=hist", q), ctx)
      expect_identical(r$state$conn_stat, "corridors")
      expect_length(r$problems, 1)
})

test_that("settings that don't apply are left out", {
      q <- url_encode(modifyList(d, list(track = "Forest type")), ctx)
      expect_false(grepl("&sp=", q, fixed = TRUE))
      expect_false(grepl("&st=", q, fixed = TRUE))
      expect_true(grepl("&ft=", q, fixed = TRUE))
      # The dispersal reach appears only where connectivity applies
      expect_true(grepl("&crch=100", url_encode(d, ctx), fixed = TRUE))
      expect_false(grepl("&crch=", url_encode(modifyList(d, list(entity_scope = "multi")), ctx), fixed = TRUE))
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

test_that("the land use weight must be one of the slider's steps", {
      for (bad in c("30", "125", "lots")) {
            r <- url_decode(paste0("?clu=", bad), ctx)
            expect_identical(r$state$cs_landuse, d$cs_landuse)
            expect_true(any(grepl(paste0("^Land use weight: '", bad, "'"), r$problems)))
      }
      expect_identical(url_decode("?clu=25", ctx)$state$cs_landuse, 25)
})

test_that("links from before the connectivity measure was retired still open cleanly", {
      r <- url_decode("?r=conn&rcm=choke", ctx)
      expect_identical(r$state$stat, "connectivity")
      expect_length(r$problems, 0)
})

test_that("dependent settings are checked against the settings they depend on", {
      # Categorical measures aren't offered as change
      r <- url_decode("?es=multi&cm=domid&t=chg", ctx)
      expect_identical(r$state$comm_measure, "richness")
      expect_length(r$problems, 1)
      # Connectivity is single-entity only, and needs a habitat variable
      r <- url_decode("?es=multi&r=conn", ctx)
      expect_identical(r$state$stat, "FIA variable")
      expect_length(r$problems, 1)
      r <- url_decode("?e=all&m=rich&r=conn", ctx)
      expect_identical(r$state$stat, "FIA variable")
      # ...and is shown for the timeframe the link gives
      r <- url_decode("?r=conn&t=fut&crch=inf", ctx)
      expect_identical(r$state[c("stat", "time", "cs_reach")],
                       list(stat = "connectivity", time = "future", cs_reach = "Inf"))
})

test_that("multi-selects are written in menu order", {
      a <- url_encode(modifyList(d, list(clim_vars = c("cwd", "tmean"))), ctx)
      b <- url_encode(modifyList(d, list(clim_vars = c("tmean", "cwd"))), ctx)
      expect_identical(a, b)
})

test_that("a storyline link keeps its model; a model alone is a specific model", {
      st <- function(q) url_decode(q, ctx)
      most <- ts_storyline_model(ctx$data, "most", d$scenario_ssp, d$scenario_period,
                                 d$bbox_x, d$bbox_y)
      least <- ts_storyline_model(ctx$data, "least", d$scenario_ssp, d$scenario_period,
                                  d$bbox_x, d$bbox_y)
      # The default link names the storyline and its model
      expect_match(url_encode(d, ctx), paste0("&sto=median&gcm=", d$scenario_model, "&"), fixed = TRUE)

      # A storyline alone resolves, for the link's bounding box
      r <- st("?sto=most")
      expect_identical(r$state[c("scenario_storyline", "scenario_model")],
                       list(scenario_storyline = "most", scenario_model = most))
      expect_length(r$problems, 0)
      far <- st("?sto=most&bx=-90,-80&by=30,40")
      expect_identical(far$state$scenario_model,
                       ts_storyline_model(ctx$data, "most", d$scenario_ssp, d$scenario_period,
                                          c(-90, -80), c(30, 40)))

      # A model alone is that model, with no storyline
      r <- st(paste0("?gcm=", most))
      expect_null(r$state$scenario_storyline)
      expect_true("scenario_storyline" %in% names(r$state))
      expect_identical(r$state$scenario_model, most)
      expect_length(r$problems, 0)
      expect_false(grepl("sto=", url_encode(r$state, ctx), fixed = TRUE))

      # A storyline with a model keeps the model, saying so if the storyline
      # would now pick another
      r <- st(paste0("?sto=most&gcm=", least))
      expect_identical(r$state$scenario_model, least)
      expect_identical(r$state$scenario_storyline, "most")
      expect_match(r$problems, "was kept", fixed = TRUE)
      expect_length(st(paste0("?sto=most&gcm=", most))$problems, 0)
})
