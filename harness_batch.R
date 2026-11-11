# ---- run_batch.R --------------------------------------------------------
# Worked examples against the real treeshift input schema. No app.R changes
# needed: ts_start() opens the Spatial grid and Export panels once so Shiny
# resumes the snap button and the download handlers.

source("treeshift_harness.R")

APP <- "."


# --- one-off --------------------------------------------------------------
# Arguments in any order; the harness sequences them. Buttons are flags.

if (FALSE) {
      app <- ts_start(APP)

      ts_render(app,
                track = "Species", species = "Pinus ponderosa",
                metric = "occurrence", stage = "any",
                focal_domain = TRUE,
                stat = "FIA variable", time = "historical",
                clim_vars = c("cwd", "aet"), theta_clim = 0.1,
                download = c("raster", "metadata"),
                outdir = "scratch", tag = "pipo_occ_hist")

      ts_state(app)     # confirm fia_var resolved as expected
      app$stop()
}


# --- constants ------------------------------------------------------------

SPECIES <- "Pinus ponderosa"

CONST <- list(
      track           = "Species",
      species         = SPECIES,
      stat            = "FIA variable",

      focal_domain    = TRUE,   # snap the grid to the focal species' range

      # Values are the raw lowercase filename tokens, not the display labels:
      # the app parses data/clim_future/<ssp>_<period>_<model>.tif and splits
      # the three axes into independent dropdowns. Not every combination
      # exists on disk -- a missing one trips shiny::validate() in
      # scenario_file() rather than erroring. Check what you actually have:
      #     list.files("data/clim_future")
      # These are the app's own defaults, so they are known to resolve.
      scenario_ssp    = "ssp585",
      scenario_period = "2041-2070",
      scenario_model  = "gfdl-esm4",

      clim_vars       = c("cwd", "aet"),   # of tmean/ppt/aet/cwd/tmincm/tmaxwm
      clim_agg        = 10                 # km; a species-range grid at 1 km is slow
)


# --- rasters: occurrence + basal area, historical + future ----------------

grid_rast <- expand.grid(
      metric     = c("occurrence", "basal area"),
      time       = c("historical", "future"),
      theta_clim = c(0.05, 0.1, 0.25),   # slider range is .05-.5
      stringsAsFactors = FALSE
)

manifest_rast <- ts_batch(
      grid      = grid_rast,
      app       = APP,
      outdir    = "batch_output/rasters",
      download  = c("raster", "metadata"),
      constants = CONST,
      tag_fun   = function(row) sprintf("pipo_%s_%s_tc%s",
                                        gsub(" ", "-", row$metric),
                                        row$time, row$theta_clim)
)

# Note: `stage` is omitted above. Occurrence defaults to the first stage
# choice ("seedling") via the metric observer, NOT to "any". Add
# stage = "any" to CONST if that is what you want -- but it is only a legal
# control for occurrence/abundance/regeneration fraction, so a mixed-metric
# grid needs it as a grid column with NA for basal area, or two batches.


# --- regeneration fraction point data -------------------------------------
# Regeneration fraction is a metric, not a vector layer; `vector` chooses
# which plot layer gets exported by download_points.

grid_pts <- expand.grid(
      stage      = c("seedling", "sapling"),
      theta_clim = c(0.05, 0.1, 0.25),
      stringsAsFactors = FALSE
)

manifest_pts <- ts_batch(
      grid      = grid_pts,
      app       = APP,
      outdir    = "batch_output/points",
      download  = c("points", "metadata"),
      constants = c(CONST, list(metric = "regeneration fraction",
                                vector = "occupied FIA plots")),
      tag_fun   = function(row) sprintf("pipo_regen-%s_tc%s",
                                        row$stage, row$theta_clim)
)


# --- check ----------------------------------------------------------------

subset(rbind(manifest_rast, manifest_pts), status != "ok")

cat(readLines(list.files("batch_output/rasters", "_metadata\\.txt$",
                         full.names = TRUE)[1]), sep = "\n")
