# The settings record (helpers/provenance.R)

load_app_globals()
ctx <- url_context()
d <- url_defaults(ctx)
versions <- list(app = "abc1234", climate = "data-v1", boundaries = "boundaries-v1.2", fia = "0123456789ab")
record <- function(s) provenance_text(s, ctx, "https://example.org/treeshift/?v=1", versions,
                                      created = as.POSIXct("2026-11-11 12:00", tz = "UTC"))

test_that("the record has the link, settings, and versions", {
      r <- record(d)
      expect_identical(r[2], "Created: 2026-11-11 12:00 UTC")
      expect_true("  https://example.org/treeshift/?v=1" %in% r)
      expect_true(any(grepl("^  Entity: +species: Pinus edulis$", r)))
      expect_true(any(grepl("^  Life stage: +all life stages$", r)))
      expect_true(any(grepl("^  TreeShift: +commit abc1234$", r)))
      expect_true(any(grepl("^  Boundary library: +boundaries-v1.2$", r)))
      expect_true(any(grepl("^  Focal polygon: +none$", r)))
      expect_true(any(grepl("^  Longitude range: +-116 to -100$", r)))
})

test_that("the record lists only settings that apply, with menu labels", {
      r <- record(modifyList(d, list(entity_scope = "multi", comm_measure = "dominant entity",
                                     comm_top_n = 6, polygon_unit = "county:08069")))
      expect_true(any(grepl("^  Entity: +species: all species$", r)))
      expect_true(any(grepl("^  Categories shown: +6$", r)))
      expect_true(any(grepl("Larimer County, Colorado \\(county:08069\\)", r)))
      expect_false(any(grepl("Connectivity", r)))

      r <- record(modifyList(d, list(track = "Forest type", forest_type = "221")))
      expect_true(any(grepl("^  Entity: +forest type: 221 .+Ponderosa pine$", r)))
      expect_false(any(grepl("Life stage", r)))
})

test_that("the app's commit is found with or without git", {
      expect_match(app_commit(app_dir), "^([0-9a-f]{7,}( \\(with uncommitted changes\\))?|unknown)$")
      # Without git, the commit is read from .git directly
      no_git <- withr::with_envvar(c(PATH = ""), app_commit(app_dir))
      expect_match(no_git, "^([0-9a-f]{7} \\(uncommitted changes not checked\\)|unknown)$")
      if (!identical(no_git, "unknown"))
            expect_identical(substr(no_git, 1, 7), substr(app_commit(app_dir), 1, 7))
})
