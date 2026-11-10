library(jsonlite)
library(shiny)
library(leaflet)
library(ggplot2)
library(patchwork)
library(terra)
library(analogs) # renv::install("matthewkling/analogs")
library(dplyr)
library(readr)
library(bslib)
library(shinybusy)
library(leaflet.extras2)
library(sf)
library(leafgl)
library(circuitscaper)


## global setup ##

# Lazy Julia/Circuitscape initialization.
#
# cs_setup() warms up the Julia runtime and JIT-compiles Circuitscape, which
# is slow (tens of seconds). Doing it here at global scope makes EVERY new R
# process pay that cost at startup before serving its first page, even those who
# never open the connectivity panel.
#
# Instead, defer cs_setup() until the first time connectivity is actually
# computed (see connectivity()), and run it at most once per R process via
# this guard. Most users never touch connectivity, so they never pay the cost;
# the first connectivity user in a given process pays it once.
#
# LONGER-TERM: the better fix is to eliminate the warmup cost rather than
# defer it, by building a custom Julia sysimage with PackageCompiler.jl that
# bakes Circuitscape's compiled code in, and pointing JuliaCall at it. That
# makes the per-process warmup near-instant, at which point eager global init
# would be fine again. Lazy init is the stopgap until then.
.cs_initialized <- FALSE
ensure_cs_setup <- function() {
      if (.cs_initialized) return(invisible(FALSE))
      cs_setup()
      .cs_initialized <<- TRUE
      invisible(TRUE)
}

# Ensure the climate rasters are present locally, fetching them from the
# GitHub release if missing. The bulky climate .tif files are not committed
# to git; they are distributed as piggyback release assets (see README).
#
# Layout reconciliation: pb_download() writes every release asset by its
# (flat) asset name into a single directory. The app, however, expects
#   data/clim_historical.tif
#   data/wind_quantized.tif
#   data/clim_future/<ssp>_<period>_<model>.tif   (many files)
# To bridge the flat release namespace and this nested layout, all assets are
# downloaded into a temp staging dir, then sorted: the named top-level rasters
# go to data/, and every other .tif is treated as a future-scenario raster and
# moved to data/clim_future/. New scenario rasters therefore need no code
# change here; only a new *non-scenario* raster would require updating
# `top_level_rasters` below.
ensure_climate_data <- function(repo = "matthewkling/treeshift",
                                tag  = "data-v1",
                                data_dir = "data") {

      future_dir <- file.path(data_dir, "clim_future")

      # Top-level (non-scenario) rasters, by asset name. Anything else with a
      # .tif extension is assumed to be a future-scenario raster.
      top_level_rasters <- c("clim_historical.tif", "wind_quantized.tif")

      # Presence test: all named top-level rasters exist, and clim_future/
      # exists and is non-empty. (The scenario set is variable, so "non-empty
      # directory" is the test rather than enumerating every file.)
      have_top_level <- all(file.exists(file.path(data_dir, top_level_rasters)))
      have_future    <- dir.exists(future_dir) &&
            length(list.files(future_dir, pattern = "\\.tif$")) > 0

      if (have_top_level && have_future) {
            return(invisible(FALSE))   # data already present; nothing to do
      }

      message("Climate rasters missing; downloading from ", repo,
              " release '", tag, "'. This may take several minutes...")

      if (!requireNamespace("piggyback", quietly = TRUE)) {
            stop("Package 'piggyback' is required to fetch climate data. ",
                 "Install it with install.packages('piggyback').")
      }

      # Stage the download in a temp dir so the sort into final locations is
      # the only step that touches data/. The staging dir is removed on exit
      # (success or error) so a failed run leaves no partial mess behind.
      staging <- tempfile("treeshift_climdata_")
      dir.create(staging)
      on.exit(unlink(staging, recursive = TRUE), add = TRUE)

      # pb_download() with no `file` argument pulls every asset on the release.
      piggyback::pb_download(repo = repo, tag = tag, dest = staging)

      dir.create(future_dir, recursive = TRUE, showWarnings = FALSE)

      staged <- list.files(staging, pattern = "\\.tif$", full.names = TRUE)
      if (length(staged) == 0) {
            stop("No .tif assets found on release '", tag, "' of ", repo, ". ",
                 "Check the release tag and asset names at ",
                 "https://github.com/", repo, "/releases/tag/", tag)
      }

      # Sort each staged raster to its destination: named top-level rasters to
      # data/, everything else to data/clim_future/.
      for (src in staged) {
            nm   <- basename(src)
            dest <- if (nm %in% top_level_rasters) {
                  file.path(data_dir, nm)
            } else {
                  file.path(future_dir, nm)
            }
            if (!file.rename(src, dest)) {
                  # file.rename() fails across filesystems (temp dir may be on
                  # a different mount); fall back to copy-then-remove.
                  ok <- file.copy(src, dest, overwrite = TRUE)
                  if (!ok) stop("Failed to place climate raster: ", nm)
                  file.remove(src)
            }
      }

      # Verify the sort produced what the app needs, so a partial or
      # mis-named download surfaces here rather than as a cryptic rast() error
      # further down in global setup.
      missing_top <- top_level_rasters[
            !file.exists(file.path(data_dir, top_level_rasters))]
      if (length(missing_top) > 0) {
            stop("Climate data download is missing expected top-level raster(s): ",
                 paste(missing_top, collapse = ", "), ". ",
                 "Check the asset names on the release.")
      }
      if (length(list.files(future_dir, pattern = "\\.tif$")) == 0) {
            stop("Climate data download produced no future-scenario rasters in ",
                 future_dir, ".")
      }

      message("Climate rasters ready.")
      invisible(TRUE)
}

# Fetch climate data if absent. Runs once at app process startup, before any
# rast("data/...") call below. On the production shiny-server deployment,
# prefer running this once as an explicit provisioning step (from an R
# console) so app startup stays fast; this in-app call then just confirms
# presence and returns immediately.
ensure_climate_data()

clim_vars <- function() c("tmean", "ppt", "aet", "cwd", "tmincm", "tmaxwm")

# Native (un-standardized) units for each climate variable, used to label the
# climate-model selector chart after deltas are back-converted from z-scores.
# The stored climate rasters are standardized; a delta back-converts to native
# units by multiplying the z-score delta by the variable's SD (the mean
# cancels in a difference). NOTE: ppt is log10-transformed before standardizing,
# so its delta is in log10-space, not mm -- labelled accordingly.
clim_units <- c(tmean  = "\u00b0C",
                tmincm = "\u00b0C",
                tmaxwm = "\u00b0C",
                aet    = "mm",
                cwd    = "mm",
                ppt    = "log10 mm")

# Native-unit label for a climate variable. Climate rasters are stored in
# standardized (z-score) units, so any plot showing native units must first
# back-convert via clim_scales (z * sd + mean). Falls back to "standardized"
# for a variable with no unit entry.
# Colors offered for the focal polygon outline. A curated set rather than a full
# picker: the polygon is a boundary drawn over a colored raster, so what matters
# is contrast against a busy background, and a dozen strong, well-separated hues
# cover that better than an unconstrained wheel. The selectize `create` option
# still accepts anything typed, so a specific brand color is one keystroke away.
overlay_colors <- c("darkred", "#e05c00", "#c2b200", "#2b6a4d", "#1f78b4",
                    "#4c32a8", "#6a3d9a", "#e7298a", "black", "white", "gray40")

# Swatch renderer for the overlay color selectors. Defined once because both use
# it, and a divergence between them would be invisible until someone compared
# the two dropdowns side by side.
overlay_color_render <- function() {
      swatch <- paste0(
            "'<span style=\"display:inline-block;width:12px;height:12px;",
            "margin-right:8px;vertical-align:middle;border:1px solid #999;",
            "background:' + escape(item.value) + ';\"></span>'")
      I(paste0("{
            option: function(item, escape) {
                  return '<div style=\"padding:4px 8px;\">' + ", swatch, " + escape(item.label) + '</div>';
            },
            item: function(item, escape) {
                  return '<div>' + ", swatch, " + escape(item.label) + '</div>';
            }
      }"))
}

# Blend a color toward white. Used for the buffer ring, which follows the focal
# polygon's color but needs to stay distinguishable from it -- a tint plus the
# dashed line reads as "the same thing, extended" rather than as a second
# unrelated boundary.
lighten_color <- function(col, amount = 0.45) {
      v <- grDevices::col2rgb(col)[, 1]
      grDevices::rgb(t(v + (255 - v) * amount), maxColorValue = 255)
}

# One selector for either overlay. `create = TRUE` accepts any R color name or
# hex typed straight in, so the curated list is a starting point rather than a
# limit.
overlay_color_input <- function(id, label, current) {
      selectizeInput(id, label,
                     choices = union(overlay_colors, current),
                     selected = current,
                     options = list(create = TRUE, render = overlay_color_render()))
}

# Variables the climate-relationships plot can put on its x axis. Deliberately
# wider than the set the model is fitted on: a relationship against a variable
# the model did not use is still a real property of the landscape, and often the
# more interesting one -- richness peaking at mid elevation says something even
# though elevation is not a predictor.
#
# Kept entirely separate from the modeling path. Widening clim_pool() instead
# would be tidier-looking and much riskier: pool()'s analog index is built from
# whatever climate columns it is handed, so a missed subset would silently build
# a six-dimensional index instead of two and change every prediction with
# nothing on screen to say so.
#
# Extra (non-climate) rasters are declared here with their file and units. Any
# single-layer raster on the same grid slots in; back_transform() already leaves
# variables with no clim_scales entry untouched, which is right for a raster
# stored in native units.
display_extra_vars <- list(
      elevation = list(file = "data/elevation.tif", units = "m")
)
display_extra_vars <- display_extra_vars[
      vapply(display_extra_vars, function(x) file.exists(x$file), logical(1))]

# Available x-axis variables, and which of them the model can be fitted on.
display_vars <- function() c(clim_vars(), names(display_extra_vars))

unit_lab <- function(v) {
      if (!is.null(display_extra_vars[[v]])) return(display_extra_vars[[v]]$units)
      u <- unname(clim_units[v])
      if (is.null(u) || is.na(u)) "standardized" else u
}
bbox_res <- 1
clim_bbox <- as.vector(ext(rast("data/clim_historical.tif")))
clim_bbox <- c(floor(clim_bbox[c(1, 3)] / bbox_res) * bbox_res, ceiling(clim_bbox[c(2, 4)] / bbox_res) * bbox_res)

# Future climate scenarios: discover available .tif files in data/clim_future/.
# Filenames follow the convention <ssp>_<period>_<model>.tif (e.g.,
# ssp585_2041-2070_gfdl-esm4.tif). Three independent dropdowns expose the
# pathway / era / model axes; the active raster filename is reconstructed from
# the current selections.
.ssp_label_map <- c(ssp126 = "SSP1-2.6", ssp245 = "SSP2-4.5",
                    ssp370 = "SSP3-7.0", ssp585 = "SSP5-8.5")
.scenario_files <- list.files("data/clim_future", pattern = "\\.tif$", full.names = FALSE)
.scenario_parts <- strsplit(tools::file_path_sans_ext(.scenario_files), "_", fixed = TRUE)
scenario_meta <- tibble::tibble(
      file   = .scenario_files,
      ssp    = vapply(.scenario_parts, `[`, character(1), 1),
      period = vapply(.scenario_parts, `[`, character(1), 2),
      model  = vapply(.scenario_parts, `[`, character(1), 3)
) %>%
      mutate(ssp_label = ifelse(ssp %in% names(.ssp_label_map),
                                .ssp_label_map[ssp], toupper(ssp)),
             model_label = toupper(model))

# Per-axis choice vectors (label → value). Sorted so SSPs ascend by forcing,
# periods ascend chronologically (string sort works for "YYYY-YYYY"), and
# models alphabetize.
.ssp_meta <- scenario_meta %>% distinct(ssp, ssp_label) %>% arrange(ssp)
.period_meta <- scenario_meta %>% distinct(period) %>% arrange(period)
.model_meta <- scenario_meta %>% distinct(model, model_label) %>% arrange(model)

ssp_choices    <- setNames(as.list(.ssp_meta$ssp), .ssp_meta$ssp_label)
period_choices <- setNames(as.list(.period_meta$period), .period_meta$period)
model_choices  <- setNames(as.list(.model_meta$model), .model_meta$model_label)

# Defaults: prefer the published baseline scenario; fall back along each axis
# independently if a piece is missing.
.pick_default <- function(preferred, axis_values) {
      if (preferred %in% axis_values) preferred else axis_values[1]
}
default_ssp    <- .pick_default("ssp585",    .ssp_meta$ssp)
default_period <- .pick_default("2041-2070", .period_meta$period)
default_model  <- .pick_default("gfdl-esm4", .model_meta$model)

fia_all <- read_csv("data/fia.csv", show_col_types = FALSE)
spp <- sort(unique(fia_all$species))
forest_types_in_data <- sort(unique(na.omit(fia_all$forest_type)))
# Initial map bounds: bounding box of all FIA plots
fia_xy_bounds <- list(
      xmin = min(fia_all$x, na.rm = TRUE),
      xmax = max(fia_all$x, na.rm = TRUE),
      ymin = min(fia_all$y, na.rm = TRUE),
      ymax = max(fia_all$y, na.rm = TRUE)
)

# Forest type / group lookup table (FIA FORTYPCD)
forest_type_codes <- read_csv("data/forest_type_codes.csv", show_col_types = FALSE) %>%
      bind_rows(tibble::tibble(
            forest_type = 0L,
            forest_type_name = "Non-forest",
            forest_group = 0L,
            forest_group_name = "Non-forest",
            forest_group_abbr = "nonfor"
      ))

# Include the Non-forest sentinel in the in-data list so it appears in the
# dropdown alongside real codes. fia_plots_species_all() recodes NA -> 0 to
# match.
forest_types_in_data <- sort(unique(c(0L, forest_types_in_data)))

# Only expose codes that actually appear in the FIA data
.ft_lookup <- forest_type_codes %>%
      filter(forest_type %in% forest_types_in_data) %>%
      arrange(forest_type_name)
.fg_lookup <- forest_type_codes %>%
      distinct(forest_group, forest_group_name) %>%
      # restrict groups to those that have at least one type present in the data
      filter(forest_group %in% (forest_type_codes %>%
                                      filter(forest_type %in% forest_types_in_data) %>%
                                      pull(forest_group) %>% unique())) %>%
      arrange(forest_group_name)

# Named choice vectors: labels shown to the user, integer codes stored as values.
forest_type_choices <- setNames(
      as.list(.ft_lookup$forest_type),
      paste0(.ft_lookup$forest_type, " — ", .ft_lookup$forest_type_name)
)
forest_group_choices <- setNames(
      as.list(.fg_lookup$forest_group),
      paste0(.fg_lookup$forest_group, " — ", .fg_lookup$forest_group_name)
)

# Lookup helpers (vectorized; return NA for unknown codes)
ft_name <- function(code) {
      forest_type_codes$forest_type_name[match(code, forest_type_codes$forest_type)]
}
fg_name <- function(code) {
      # use distinct group rows to avoid duplicate matches
      g <- distinct(forest_type_codes, forest_group, forest_group_name)
      g$forest_group_name[match(code, g$forest_group)]
}

cs_stat_names <- c("current flow", "current on target", "current loss", "current pinch points",
                   "current direction", "voltage", "conductance", "loss")
cs_stats <- c("current_flow", "current_destinations", "current_loss", "current_power",
              "current_direction", "voltage", "conductance", "loss")

# Left-sidebar menu choices, in connectivity-off and connectivity-on variants.
# Connectivity is a niche module, so its entries are kept out of the core
# `stat` and `vector` menus unless the user flips "Enable connectivity" in the
# right-sidebar Connectivity panel. The `*_base` sets are the UI defaults; the
# observeEvent(input$cs_enable) handler in the server swaps between the two.
# Listed explicitly rather than appended so "none" stays last in the vector menu.
stat_choices_base <- c("FIA variable", "effective sample size" = "ESS",
                       "focal site similarity")
stat_choices_cs   <- c("FIA variable", "effective sample size" = "ESS",
                       "focal site similarity", "connectivity")

vector_choices_base <- c("occupied FIA plots", "all FIA plots", "focal site analogs",
                         "none")
vector_choices_cs   <- c("occupied FIA plots", "all FIA plots", "focal site analogs",
                         "current animation", "current field", "none")

## Variables available for a track AND life stage. Stage is now the
## upstream control (see the UI comment), so this is where the
## stage-specific gaps in the FIA data get expressed:
##
##   basal area              only under "any" -- FIA's basal_area is not
##                           broken out by size class
##   proportion basal area   same
##   abundance               needs a specific stage; there is no combined
##                           stem count across seedlings, saplings and trees
##   regeneration fraction   seedling or sapling only; it is a ratio
##                           AGAINST adults, so "adult" and "any" are
##                           meaningless
metric_choices <- function(track, stage = "any") {
      if (track %in% c("Forest type", "Forest group")) return("occurrence")
      any_stage <- identical(stage, "any")

      c("occurrence",
        if (!any_stage) "abundance",
        if (any_stage) "basal area",
        if (identical(track, "All trees")) "species richness",
        if (any_stage && identical(track, "Species")) "proportion basal area",
        if (stage %in% c("seedling", "sapling")) "regeneration fraction")
}

# There is no stage_choices() counterpart: since the dependency was inverted,
# stage is the upstream control and its options are the same four everywhere.
# The observer that used to repopulate them from `metric` is gone -- it was the
# old backwards direction, and it also overwrote the labeled choice vector
# ("all life stages") with bare values and reset the user's selection.

# Canonical fia_var for a (track, metric, stage) combination. Returns NULL for
# combinations metric_choices() would not have offered -- "abundance" has no
# "any" stage, for instance -- so callers must guard rather than assume a value.
resolve_var <- function(track, metric, stage) {
      if (track %in% c("Forest type", "Forest group")) return("presence probability")
      switch(metric,
             "occurrence" = switch(stage,
                                   "seedling" = "seedling presence",
                                   "sapling"  = "sapling presence",
                                   "adult"    = "tree presence",
                                   "any"      = "presence probability"),
             "abundance" = switch(stage,
                                  "seedling" = "number of seedlings",
                                  "sapling"  = "number of saplings",
                                  "adult"    = "number of trees"),
             "basal area"            = "total basal area",
             "species richness" = switch(stage,
                                         "seedling" = "seedling species richness",
                                         "sapling"  = "sapling species richness",
                                         "adult"    = "adult species richness",
                                         "any"      = "species richness"),
             "proportion basal area" = "proportion basal area",
             "regeneration fraction" = switch(stage,
                                              "seedling" = "seedling fraction",
                                              "sapling"  = "sapling fraction"))
}

# Climate variables selected at startup. Named once because three controls
# depend on it: clim_vars itself, and cr_v1, whose choices an observer then
# keeps in sync with it. Restating the default in each UI declaration is
# how the Variable dropdown ended up offering "abundance" at a life stage where
# it could not resolve.
clim_vars_default <- c("cwd", "aet")

# Occurrence-type variables. Richness is a sum of occurrence probabilities and
# is only defined over these: summing "total basal area" across entities gives
# community basal area, which is a perfectly good surface but is not richness.
occurrence_vars <- c("presence probability", "seedling presence",
                     "sapling presence", "tree presence")

# Variables whose value depends on the SUPPORT it was estimated over, and so
# need the plot-scale qualifier in the legend.
#
# Densities are not here: a kernel-weighted mean of per-acre basal area is the
# expected per-acre basal area whatever area you evaluate it over, so support
# affects their variance but not their expectation.
#
# Occurrence is different in kind, because presence is a threshold on an
# underlying density -- P(>=1 stem in 0.07 ha) and P(>=1 stem in a 10 km cell)
# are different numbers with no unit conversion between them.
#
# Two non-obvious members. "proportion basal area" is a ratio of means, which
# would be scale-free on its own, but combine_ratio() multiplies it by the
# occurrence layer, and that multiplier carries the support dependence. The
# regeneration fractions are ratios of occurrence-based counts, and the support
# does not cancel: both components saturate toward 1 as support grows, dragging
# the ratio toward 0.5.
support_dependent_vars <- c("presence probability", "seedling presence",
                            "sapling presence", "tree presence",
                            "proportion basal area",
                            "seedling fraction", "sapling fraction",
                            "species richness", "seedling species richness",
                            "sapling species richness", "adult species richness")

# Paired measures: computed from the historical AND future stacks together, and
# undefined for a single timeframe. They ignore input$time rather than
# inheriting historical/future/delta from it.
#
# Note these are NOT deltas of a per-timeframe measure, and the distinction
# matters: net richness change of zero is entirely compatible with complete
# species replacement. Both families have to exist.
#
# The pairing rule is a separate axis from the measure. Every formula here takes
# two composition vectors and does not care whether they are separated in time
# (same cell, two timeframes), in space (focal site vs every cell, one
# timeframe), or both (the compositional counterpart of the reverse/forward
# analog machinery). Only temporal pairing is wired up so far, but the measures
# are written so the others slot in without redefining them.
paired_measures <- c("turnover", "expected gains", "expected losses",
                     "colonization demand", "loss fraction",
                     "biggest winner", "biggest loser")

# Measures whose value is an entity identity rather than a quantity. They need
# a discrete palette, nearest-neighbour reprojection, and no scale transform,
# and they are meaningless as a difference -- so they never take the "change in
# X" relabelling the continuous per-timeframe measures get.
#
# Note "dominant entity" is NOT restricted to co-occurring entities: over a
# mutually exclusive classification the argmax is the most likely forest type,
# which is a perfectly good map.
# "biggest winner" and "biggest loser" are BOTH categorical and paired: their
# value is an entity identity, and they are defined only across two timeframes.
# Membership in paired_measures puts them under "change"; membership here gives
# them the discrete palette and nearest-neighbour reprojection.
categorical_measures <- c("dominant entity", "biggest winner", "biggest loser")

# Palette choices for the raster legend. The two sets are disjoint by design:
# a sequential or diverging ramp applied to entity identities implies an
# ordering that does not exist, and a qualitative set applied to a continuous
# surface is unreadable. One control swaps between them.
#
# The categorical options come from grDevices::palette.colors() rather than
# hcl.colors(): the hcl qualitative palettes sample a 1D path through HCL space,
# so past about eight categories adjacent colors differ only slightly and the
# result is a discretized gradient. scales::hue_pal() has the same problem --
# a hue sweep at fixed chroma and luminance.
#
# Polychrome 36 and Alphabet are built by maximizing perceptual distance across
# hue, chroma AND lightness together, which is what actually survives n = 10-25.
# Tableau 10 and Set 3 look better at small n but run out; Okabe-Ito is the
# colorblind-safe choice at eight. All are base R since 4.0, no new dependency.
rast_palettes_continuous  <- c("default", "viridis", "magma", "plasma", "inferno",
                               "RdYlBu", "Spectral", "RdBu", "BrBG")
rast_palettes_categorical <- c("default", "Polychrome", "Alphabet", "Tableau",
                               "Set 3", "Okabe-Ito")

# Analog impact surfaces --------------------------------------------------------
#
# AIM surfaces are computed on plot and cell coordinates in CONUS Albers
# equal-area (EPSG:5070), in km so the geographic bandwidth keeps its units. The
# FFT method needs planar coordinates, and planar distance also makes the exact
# method about 1.5x faster than great-circle distance; over CONUS the two differ
# by about 1%. Results are written back onto the lon/lat modeling grid.
aim_crs <- "EPSG:5070"

aim_xy <- function(xy) {
      p <- terra::project(as.matrix(xy[, c("x", "y")]), "EPSG:4326", aim_crs) / 1000
      colnames(p) <- c("x", "y")
      p
}

aim_kernels <- function(theta_clim, theta_geog) {
      list(env  = analogs::kernel(max = theta_clim * 3, weight = "gaussian", theta = theta_clim),
           geog = analogs::kernel(max = theta_geog * 3, weight = "gaussian", theta = theta_geog))
}

# FFT runs are tiled geographically, so each tile's grid spans only the climates
# found there. ~300 km tiles were fastest in benchmarks.
aim_tile_km <- 300

# Which method a run uses: "exact" visits every cell-plot pair, "fft" (binned
# convolution) costs about the same whatever the grain, and its error is small
# (99.9th percentile ~0.007 on a 0-1 occurrence scale; identical NA cells). Exact
# wins at coarse grain or narrow bandwidths, FFT at fine grain and wide ones.
#
# The choice must be a deterministic function of the settings, so that the same
# settings always give the same map -- hence estimated work, not timings. Exact
# work is the number of cell-plot pairs within the kernel limits, from analog
# counts at a fixed sample of cells. FFT work is the number of grid nodes, from
# each tile's extent in geography and climate (nodes are a quarter-bandwidth
# apart, plus a 3-bandwidth halo), times the number of kernel sums. The
# coefficients were fitted on runs at 1-10 km over 6-20 degree domains with
# 1-80 outcome columns; only their ratio matters, so machine speed cancels.
# One method serves every timeframe of a run, so differences between timeframes
# never mix methods.
aim_method <- function(pts_list, index, theta_clim, theta_geog, n_y, n_sums) {
      env_vars <- setdiff(names(pts_list[[1]]), c("x", "y"))
      if (length(env_vars) > 3) return("exact")   # FFT supports up to 3

      est <- vapply(pts_list, function(pts) {
            n <- nrow(pts)
            probe <- pts[unique(round(seq(1, n, length.out = min(n, 2000)))), , drop = FALSE]
            count <- suppressWarnings(analogs::analog_search(
                  x = probe, pool = index, select = "all", stat = "count",
                  env = analogs::kernel(max = theta_clim * 3),
                  geog = analogs::kernel(max = theta_geog * 3)))$count
            pairs <- n * mean(count, na.rm = TRUE)

            spans <- pts %>%
                  group_by(tx = floor(x / aim_tile_km), ty = floor(y / aim_tile_km)) %>%
                  summarise(across(everything(), ~ max(.x) - min(.x)), .groups = "drop")
            nodes_along <- function(span, theta) (span + 6 * theta) / (0.25 * theta) + 1
            nodes <- sum(nodes_along(spans$x, theta_geog) * nodes_along(spans$y, theta_geog) *
                         Reduce(`*`, lapply(spans[env_vars], nodes_along, theta = theta_clim)))

            c(exact = 5.0e-6 * n + 4.5e-9 * pairs * (1 + 0.20 * (n_y - 1)),
              fft   = 1.8e-6 * n + 2.1e-9 * nodes * (1 + 0.42 * n_sums))
      }, numeric(2))
      if (sum(est["fft", ]) < sum(est["exact", ])) "fft" else "exact"
}

# AIM surfaces for several focal rasters (e.g. historical and future climate)
# against one pool, on each raster's own grid. `pool` is pool()'s list, of which
# this uses the AIM-coordinate plot table (`aim_pts`) and its index
# (`aim_index`). Returns a list of SpatRasters named like `focals`, with the
# layers analog_impact() would give for a raster.
aim_surfaces <- function(focals, pool, y, theta_clim, theta_geog, stat) {
      k <- aim_kernels(theta_clim, theta_geog)
      cells <- lapply(focals, function(r) as.data.frame(r, xy = TRUE, cells = TRUE, na.rm = TRUE))
      pts <- lapply(cells, function(fd) data.frame(aim_xy(fd), fd[, setdiff(names(fd), c("cell", "x", "y")), drop = FALSE]))

      run <- function(method) lapply(pts, function(p) {
            if (method == "exact") {
                  suppressWarnings(analogs::analog_impact(
                        x = p, pool = pool$aim_index, env = k$env, geog = k$geog,
                        stat = stat, y = y, coord_type = "projected"))
            } else {
                  span <- vapply(p[, c("x", "y")], function(v) diff(range(v)), numeric(1))
                  suppressWarnings(analogs::tiled_analog_search(
                        x = p, pool = pool$aim_pts, y = y,
                        n_tiles = max(1, round(prod(span) / aim_tile_km^2)),
                        fun = analogs::analog_impact, env = k$env, geog = k$geog,
                        stat = stat, method = "fft", coord_type = "projected",
                        progress = FALSE))
            }
      })

      # Kernel sums FFT computes: one per outcome column, the weights, and the
      # squared weights if ess is requested
      method <- aim_method(pts, pool$aim_index, theta_clim, theta_geog,
                           n_y = ncol(y), n_sums = ncol(y) + 1 + ("ess" %in% stat))
      res <- if (method == "fft") tryCatch(run("fft"), error = function(e) {
            # A tile whose FFT grid exceeds the size limit (narrow bandwidths
            # over climatically varied ground) sends the whole run to exact,
            # rather than mixing methods between tiles. Any other error is real.
            if (!grepl("max_cells", conditionMessage(e))) stop(e)
            message("AIM: FFT grid too large for a tile; using the exact method")
            NULL
      })
      if (is.null(res)) res <- run("exact")

      Map(function(r, fd, out) {
            cols <- setdiff(names(out), c("index", "x", "y"))
            s <- terra::rast(r, nlyrs = length(cols), names = cols, vals = NA)
            s[fd$cell] <- as.matrix(out[, cols, drop = FALSE])
            s
      }, focals, cells, res)
}

# Ceiling on a community-composition run, in gigabytes of per-entity surface
# held at once (both timeframes).
#
# The cost scales as cells x entities, and BOTH inputs are user-controlled: the
# grain slider goes to 1 km and the bounding box can span the data. At the
# default 10 km the stack is about 0.07 GB and lives in memory; at 1 km over the
# same box it is nearly 7 GB. On a shared server that is not a slow render, it
# is an out-of-memory failure that takes other sessions with it.
#
# Set generously enough to allow 2 km over a typical domain and to refuse 1 km
# over a wide one. Deployment-specific: raise it on a machine with headroom.
comm_max_gb <- 1.5

# Community measures that require entities able to CO-OCCUR at a plot.
#
# Forest types and groups are mutually exclusive plot classifications: each plot
# carries exactly one FORTYPCD, so multi_values() rows are one-hot. Richness is
# then 1 by construction, diversity is 1, dominance is 1 and evenness is 0 --
# every per-timeframe measure is degenerate.
#
# This is the plot-support constraint again rather than anything wrong with the
# classification: forest-type richness is meaningful over a cell, which contains
# many plots of potentially different types. It is degenerate only because the
# support is one plot.
#
# The paired measures will NOT be restricted this way. Bray-Curtis between two
# probability vectors over forest types measures how the type distribution
# shifts, which is well defined for a mutually exclusive classification.
cooccurrence_measures <- c("richness", "diversity", "evenness", "dominance")

# Community-level reductions across the per-entity stack. Each is a function of
# the whole entity set rather than of any one entity, which is why they occupy
# the entity-picker slot rather than the `metric` dropdown: in All mode, "how
# are entities combined" IS the selection the picker makes.
#
# Two families, distinguished by whether they are computed per timeframe (and so
# inherit historical / future / delta from input$time) or are inherently paired
# and defined only as a change. Only the per-timeframe family is implemented so
# far; the paired measures (turnover, expected gains and losses) are next.
#
# Richness is the one linear reduction here, and linearity means it commutes
# with the interpolation -- it is also available far more cheaply as an
# All-trees variable. It is kept here because it belongs to this family
# conceptually, and because the two routes agreeing is a free check on the whole
# multi-entity path. The other three genuinely require per-entity surfaces:
# no plot-level aggregate can be interpolated to recover them.
comm_measures <- c("richness", "diversity", "evenness", "dominance",
                   "dominant entity", paired_measures)

# Label and step for the support floor, which are properties of the DENOMINATOR
# rather than of the layer being masked -- and the denominator changes with the
# variable, so a fixed label would be wrong most of the time.
#
# Single-entity ratios (the regeneration fractions) are conditional on the
# entity's occurrence probability, so the floor is on 0-1. Community measures
# are normalized by the summed per-entity response, whose units are whatever the
# variable measures.
support_floor_spec <- function(multi, fia_var, track) {
      # Defaults are nonzero because the unmasked layers are actively
      # misleading: which.max always names a winner, and a ratio is at its most
      # extreme where its denominator is smallest, so the empty ground is
      # exactly where these maps look most dramatic. The values below are
      # "barely anything there" in each unit -- deliberately low, since the
      # control exists to be tuned.
      if (!isTRUE(multi))
            return(list(label = "Minimum occurrence probability",
                        step = 0.05, default = 0.1))

      unit <- switch(fia_var,
                     "total basal area"      = list("Minimum ft\u00b2 BA / acre", 5, 5),
                     "number of trees"       = list("Minimum trees / acre", 5, 5),
                     "number of saplings"    = list("Minimum saplings / acre", 5, 5),
                     "number of seedlings"   = list("Minimum seedlings / acre", 10, 10),
                     "proportion basal area" = list("Minimum summed proportion", 0.1, 0.1),
                     # Occurrence variables: the denominator is expected richness.
                     # 0.5 rather than 1: a cell with 0.6 expected species is
                     # more likely than not to hold a tree, so cutting at 1 would
                     # discard genuinely occupied sparse woodland.
                     list(paste("Minimum expected",
                                if (identical(track, "Species")) "species" else "classes"),
                          0.5, 0.5))
      list(label = unit[[1]], step = unit[[2]], default = unit[[3]])
}

# Mask a ratio layer where the community or entity it describes is barely
# present. Ratios stay numerically well behaved as their denominator approaches
# zero but stop meaning anything: evenness, dominance, the regeneration
# fractions -- and, shortly, Bray-Curtis turnover and colonization demand --
# take their most extreme values exactly where nothing is there to be even,
# dominant, regenerating or turning over.
#
# combine_ratio() already masks the fractions at denom == 0; this is the same
# idea extended to denom ~ 0, which is the case that actually produces the
# misleading map.
#
# Masking rather than fading, deliberately: the extreme values from empty cells
# also drive seq_domain(), compressing the well-supported part of the map into a
# narrow band of the palette. Removing them fixes the legend as much as the map.
apply_support_floor <- function(x, support, floor) {
      if (!isTruthy(floor) || is.na(floor) || floor <= 0) return(x)
      terra::ifel(support >= floor, x, NA)
}

# Row-wise version of the community reductions in composition_layer(), for a
# plot x entity matrix. Used to paint the FIA markers with the OBSERVED
# counterpart of whatever the raster is modelling, on the same scale.
#
# Kept deliberately parallel to the raster path -- same formulas, different
# substrate -- since the two being comparable is what makes the markers a check
# on the surface rather than decoration.
#
# One interpretive wrinkle: for an occurrence variable the observed matrix is
# 0/1, so observed diversity collapses to observed richness and observed
# dominance to 1/richness. The MODELLED versions do not, because they act on
# continuous probabilities. The markers and the raster are both correct and are
# not measuring quite the same thing there.
comm_reduce_rows <- function(y, measure) {
      # Force the argument in this frame rather than inside as.matrix(): if a
      # caller passes a reactive whose req() fails, the resulting condition
      # escaping from S4 dispatch is reported as an evaluation error instead of
      # a silent cancellation.
      force(y)
      if (is.null(y) || !ncol(y)) return(NULL)
      v <- as.matrix(y)
      v[is.na(v)] <- 0
      s <- rowSums(v)

      # Categorical: which entity is largest at the plot. Returned as a column
      # index into `y`, the same index space the raster's which.max() uses, so
      # the two share a code vocabulary before lumping.
      if (measure %in% categorical_measures) {
            # Only the per-timeframe one has an observed counterpart: a plot is
            # surveyed at a single point in time, so there is no observed
            # "biggest winner" to paint the markers with.
            if (measure %in% paired_measures) return(NULL)
            out <- max.col(v, ties.method = "first")
            out[s == 0] <- NA_integer_
            return(out)
      }

      out <- if (identical(measure, "richness")) {
            s
      } else {
            # H = -sum(q log q) with q = v / s; see composition_layer().
            vlv <- rowSums(ifelse(v > 0, v * log(v), 0))
            H <- log(s) - vlv / s
            switch(measure,
                   "diversity" = exp(H),
                   "evenness"  = H / log(ncol(v)),
                   "dominance" = do.call(pmax, as.data.frame(v)) / s,
                   NULL)
      }
      if (is.null(out)) return(NULL)

      # A plot with nothing recorded has no composition to describe -- diversity,
      # evenness and dominance are all undefined when the community is empty.
      #
      # Richness is the exception: zero is a real, well-defined observation, and
      # it is exactly the value the model should be predicting low at. NA-ing
      # those plots dropped the entire bottom of the calibration curve and pulled
      # them out of the cross-validation pool, which is also why community
      # richness disagreed with the All-trees species-richness variable -- that
      # one counts zero correctly.
      if (!identical(measure, "richness")) out[s == 0] <- NA_real_
      out
}

# Which community measures are valid for a given track and variable. Invalid
# combinations are made unselectable rather than selectable-and-warned: these
# are permanent incompatibilities, not preconditions the user could satisfy by
# supplying something. (Contrast eval_scope, where an unavailable option means
# "no polygon uploaded yet" and saying so is more useful than hiding it.)
comm_measure_choices <- function(track, fia_var, time = "historical") {
      if (!isTruthy(fia_var) || !fia_var %in% multi_vars) return(character(0))

      static <- setdiff(comm_measures, paired_measures)
      # Mutually exclusive classifications make every per-timeframe measure
      # degenerate.
      if (track %in% c("Forest type", "Forest group"))
            static <- setdiff(static, cooccurrence_measures)
      # Richness counts entities, so it needs an occurrence variable: the same
      # sum over an abundance variable is community abundance.
      if (!fia_var %in% occurrence_vars) static <- setdiff(static, "richness")

      # Timeframe is the first cut, because "state or change" is a smaller and
      # more natural question than picking out of one mixed list. Under a single
      # timeframe only the static measures apply; under change they are joined by
      # the paired ones, and relabelled so the whole list reads as change
      # measures. Only the LABELS change -- the values stay the measure ids, so
      # nothing downstream has to know about this.
      if (!identical(time, "delta")) return(stats::setNames(static, static))

      # Categorical measures are dropped rather than relabelled: the difference
      # of two entity identities is not an entity identity. Their change
      # counterparts are separate measures (biggest winner, biggest loser).
      chg <- setdiff(static, categorical_measures)
      c(stats::setNames(chg, paste("change in", chg)),
        stats::setNames(paired_measures, paired_measures))
}

# FIA variables that generalize to a per-entity column in the multi-entity
# response matrix. The ratio variables (seedling fraction, sapling fraction)
# carry the occurrence-gated numerator/denominator structure that
# ratio_pieces() hands to aim() as three coupled columns; generalizing those
# multiplies the column count by three and is deferred.
#
# Note "proportion basal area" IS here: its denominator is the plot total,
# shared across entities, so it is a plain per-entity column divided by a
# common vector rather than a pooled ratio.
multi_vars <- c("presence probability", "seedling presence", "sapling presence",
                "tree presence", "total basal area", "proportion basal area",
                "number of trees", "number of saplings", "number of seedlings")

# Spatial scopes the model-evaluation module can score over. These are the
# EVALUATION set, never the analog pool: predictions are always built from
# every available analog, and the scope only decides which held-out plots the
# metrics are computed from. "Best bandwidth here" and "best bandwidth overall"
# are different answers, and this is what lets the module give the first one.
#
# Deliberately not the same three scopes as model_delta_region: a focal site is
# a point, and a point is not an evaluation set. The meaningful site-scoped
# quantity is the plots within the geographic cutoff of the click -- the same
# neighborhood the model draws its analogs from.
# The areas a module can summarize over, named once so a code means the same
# thing everywhere. This previously diverged: "site" meant the analog
# neighborhood in the evaluation and climate-relationship modules and the single
# clicked cell in the composition and climate-model modules. Nothing shared them
# yet, so nothing broke -- but the first helper used by both would have silently
# done the wrong one.
#
# "site" is the neighborhood; "cell" is the one clicked cell.
area_names <- c("Model bounding box"      = "domain",
                "Focal polygon"           = "polygon",
                "Buffered polygon"        = "buffered",
                "Focal site neighborhood" = "site",
                "Focal site"              = "cell")

# Modules take the subset that is meaningful for them. Evaluation and climate
# relationships exclude the single cell: one plot is not an evaluation set, and
# one cell is not a relationship.
eval_scope_names <- area_names[area_names != "cell"]

# Below this many held-out plots, metrics are not reported at all. A small unit
# can easily contain single-digit plot counts, and an AUC over eight plots is
# worse than no answer because it looks like one.
eval_min_n <- 10

# Display labels for the metrics cv_performance() returns, in the order they
# should be read. Both metric sets are listed; only the applicable ones appear
# for a given response (binary -> n/auc/tss/tss_threshold/brier,
# continuous -> n/rmse/mae/bias/r2).
eval_metric_labels <- c(n             = "n",
                        auc           = "AUC",
                        tss           = "TSS",
                        tss_threshold = "TSS thresh.",
                        brier         = "Brier",
                        rmse          = "RMSE",
                        mae           = "MAE",
                        bias          = "bias",
                        r2            = "R\u00b2")

# Fields the model-evaluation results depend on, and how to name them to the
# user. Note clim_agg is absent: eval_var_set() extracts climate from the raw
# rasters at plot coordinates, never from the aggregated modeling grid, so the
# grain slider does not affect evaluation at all.
eval_stamp_labels <- c(track      = "focal entity type",
                       focal      = "focal entity",
                       measure    = "community measure",
                       fia_var    = "FIA variable",
                       bbox       = "model bounding box",
                       clim_vars  = "climate variables",
                       theta_clim = "climate bandwidth",
                       theta_geog = "geographic bandwidth",
                       polygon    = "focal polygon",
                       site       = "focal site")

theme_set(theme_bw(base_size = 16))

# House style for facet strips: white text on a black background. Defined
# once here and added to every faceted plot via `+ strip_style`. It is a
# partial theme() (only the two strip properties), so each plot keeps its
# own legend / axis theme settings and simply layers this on top.
strip_style <- theme(strip.text = element_text(color = "white"),
                     strip.background = element_rect(fill = "black", color = "black"))


# HELPER FUNCTIONS ------------------------------------

## STALENESS -------------------------------------------
##
## Several outputs in this app are computed on an explicit button rather than
## reactively, which means they can outlive the settings that produced them.
## The mechanism is the same in each case: snapshot the inputs a result depends
## on at compute time ("stamp"), compare against the current inputs on every
## render, and mark the output rather than destroying it.
##
## Marking, not clearing, is deliberate. These results can take minutes to
## produce and stay readable after the settings move; blanking them throws away
## something the user may still want. Fading plus an explicit warning makes the
## result unmistakably provisional while leaving it legible.

# Which stamped fields differ between a stored stamp and the current one.
# `exempt` names fields that were deliberately varied by the run itself and so
# must not count as staleness -- a bandwidth sweep whose whole purpose is
# comparing bandwidths would otherwise mark itself stale the moment the user
# acts on its recommendation.
stamp_diff <- function(now, old, exempt = character()) {
      if (is.null(old)) return(character())
      keys <- setdiff(intersect(names(now), names(old)), exempt)
      keys[!vapply(keys, function(k) identical(now[[k]], old[[k]]), logical(1))]
}

# Human-readable phrase for a set of changed fields.
stale_phrase <- function(fields, labels) {
      nm <- unname(labels[fields])
      nm <- nm[!is.na(nm)]
      if (!length(nm)) return("settings")
      if (length(nm) == 1) return(nm)
      paste(paste(utils::head(nm, -1), collapse = ", "), "and", utils::tail(nm, 1))
}

# Inline warning banner. Returns NULL when nothing has changed, so it can be
# dropped straight into a renderUI.
stale_banner <- function(fields, labels,
                         action = "Re-run to update these results.") {
      if (!length(fields)) return(NULL)
      div(class = "ts-notice",
          sprintf("These results were produced before the %s changed. %s",
                  stale_phrase(fields, labels), action))
}

# Fade a ggplot and caption it when stale. Implemented as an in-plot overlay
# rather than CSS opacity so it needs no JS round-trip and survives into a
# downloaded image. Added last, so it sits above the plot's own layers.
gg_stale <- function(p, stale, label = "Stale \u2014 re-run to update") {
      if (!isTRUE(stale)) return(p)
      p +
            annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                     fill = "white", alpha = 0.5) +
            labs(caption = label) +
            # Same color as the .ts-notice rule, so an in-plot caption and a
            # sidebar notice read as the same kind of statement.
            theme(plot.caption = element_text(color = "#9a3412", face = "bold",
                                              hjust = 0, size = 11))
}

circleLegend <- function(title, labels, colors, radii = 6) {
      # If only one radius is provided, repeat it for all labels
      if(length(radii) == 1) radii <- rep(radii, length(labels))

      items <- Map(function(lbl, col, rad) {
            # Use the 'rad' variable to define SVG dimensions and circle size
            sprintf(
                  "<div style='display:flex; align-items:center; margin:4px 0;'>
         <svg width='%d' height='%d' style='flex:0 0 auto;'>
           <circle cx='%d' cy='%d' r='%d' fill='%s'></circle>
         </svg>
         <span style='margin-left:8px;'>%s</span>
       </div>",
                  12, 12, # Keep container size consistent for alignment
                  6, 6,   # Center the circle in the 12x12 box
                  rad, col, lbl
            )
      }, labels, colors, radii)

      HTML(sprintf(
            "<div class='leaflet-control'
          style='background: rgba(255,255,255,0.75);
                 padding: 10px 12px;
                 margin: 2px;
                 border-radius: 5px;
                 box-shadow: 0 0 15px rgba(0,0,0,0.2);
                 color: #444;
                 line-height: 18px;'>
       <div style='font-weight:600; font-size:14px; margin-bottom:8px;'>%s</div>
       %s
     </div>",
            title,
            paste(items, collapse = "")
      ))
}

geo_circle <- function(pts, radius, thresh = .1) {
      radius <- radius * 1000 # m to km
      thresh <- thresh * 1000

      angles <- seq(from = 0, to = 360, by = 2)
      vertices <- geosphere::destPoint(p = pts, b = angles, d = radius)

      if(pts[1,1] > thresh) vertices[,1] <- ifelse(vertices[,1] > 0, vertices[,1], vertices[,1] + 360)
      if(pts[1,1] < -thresh) vertices[,1] <- ifelse(vertices[,1] < 0, vertices[,1], vertices[,1] - 360)

      vertices
}

clim_circle <- function(pt, radius, means, sds) {
      angles <- seq(from = 0, to = 360, by = 2) / 180 * pi
      as.data.frame(cbind(v1 = pt[1] * sds[1] + means[1] + cos(angles) * radius * sds[1],
                          v2 = pt[2] * sds[2] + means[2] + sin(angles) * radius * sds[2]))
}


gradient_bearing <- function(target) {
      sobel_x <- matrix(c(-1, 0, 1,
                          -2, 0, 2,
                          -1, 0, 1), nrow=3, byrow=TRUE)

      sobel_y <- matrix(c( 1,  2,  1,
                           0,  0,  0,
                           -1, -2, -1), nrow=3, byrow=TRUE)

      gx <- focal(target, w = sobel_x)
      gy <- focal(target, w = sobel_y)

      (atan2(gy, gx) * 180 / pi + 360) %% 360
}


aligned_flux <- function(bearing, wind) {
      crs(bearing) <- crs(wind)

      # 1. Define the bin centers of your windscape layers
      # Assuming 8 layers: N, NE, E, SE, S, SW, W, NW
      wind <- wind[[c("N", "NE", "E", "SE", "S", "SW", "W", "NW")]]
      bin_centers <- seq(0, 315, by = 45)

      # 2. Normalize target bearing to [0, 360)
      bearing <- bearing %% 360

      # 3. Calculate weights and indices for the two layers
      # These are rasters
      weight_upper <- (bearing %% 45) / 45
      weight_lower <- 1 - weight_upper
      idx_lower <- floor(bearing / 45) + 1
      idx_upper <- (idx_lower %% 8) + 1 # Wraps 8 back to 1

      # 4. Use app() with a vectorized indexer
      # We combine the wind layers and the indices into one stack
      combined_stack <- c(wind, idx_lower, idx_upper, weight_lower, weight_upper)

      # The function now treats 'x' as a matrix where each column is a layer
      out <- app(combined_stack, function(x) {
            # x[, 1:8] are the wind layers
            # x[, 9] is idx_lower
            # x[, 10] is idx_upper
            # x[, 11] is weight_lower
            # x[, 12] is weight_upper

            # Vectorized indexing to pick the right wind values for every cell
            # We use a matrix-style index [row_index, col_index]
            rows <- 1:nrow(x)
            val_low <- x[cbind(rows, x[, 9])]
            val_up  <- x[cbind(rows, x[, 10])]

            # Weighted average
            return(val_low * x[, 11] + val_up * x[, 12])
      })

      return(out)
}

format_velocity_data <- function(curr_mag, flow_dir) {

      u <- curr_mag * cos(flow_dir * pi / 180)
      v <- curr_mag * sin(flow_dir * pi / 180)

      ext <- terra::ext(u)
      dims <- dim(u)

      # Extract values row-major, top to bottom
      u_vals <- as.vector(t(terra::as.matrix(u, wide = TRUE)))
      v_vals <- as.vector(t(terra::as.matrix(v, wide = TRUE)))

      # Replace NA with 0
      u_vals[is.na(u_vals)] <- 0
      v_vals[is.na(v_vals)] <- 0

      list(
            list(
                  header = list(
                        parameterCategory = 2,
                        parameterNumber = 2,
                        lo1 = ext$xmin,
                        la1 = ext$ymax,
                        lo2 = ext$xmax,
                        la2 = ext$ymin,
                        dx = terra::res(u)[1],
                        dy = terra::res(u)[2],
                        nx = dims[2],
                        ny = dims[1]
                  ),
                  data = u_vals
            ),
            list(
                  header = list(
                        parameterCategory = 2,
                        parameterNumber = 3,
                        lo1 = ext$xmin,
                        la1 = ext$ymax,
                        lo2 = ext$xmax,
                        la2 = ext$ymin,
                        dx = terra::res(u)[1],
                        dy = terra::res(u)[2],
                        nx = dims[2],
                        ny = dims[1]
                  ),
                  data = v_vals
            )
      ) |>
            jsonlite::toJSON(auto_unbox = TRUE, digits = 6)
}

make_spokes <- function(mag, dir, subsample = 10, arrow_scale = 0.5) {
      mag_sub <- terra::aggregate(mag, fact = subsample, fun = "mean", na.rm = TRUE)
      dir_sub <- terra::aggregate(dir, fact = subsample, fun = "mean", na.rm = TRUE)

      coords <- terra::xyFromCell(mag_sub, 1:ncell(mag_sub))
      m <- values(mag_sub, mat = FALSE)
      d <- values(dir_sub, mat = FALSE) / 180 * pi # to radians

      valid <- !is.na(m) & !is.na(d)

      spokes <- data.frame(
            x0 = coords[valid, 1],
            y0 = coords[valid, 2],
            x1 = coords[valid, 1] + cos(d[valid]) * m[valid] * arrow_scale,
            y1 = coords[valid, 2] + sin(d[valid]) * m[valid] * arrow_scale,
            mag = m[valid]
      )

      lines <- lapply(seq_len(nrow(spokes)), function(i) {
            st_linestring(matrix(c(spokes$x0[i], spokes$y0[i],
                                   spokes$x1[i], spokes$y1[i]),
                                 ncol = 2, byrow = TRUE))
      })

      points <- lapply(seq_len(nrow(spokes)), function(i) {
            st_point(c(spokes$x0[i], spokes$y0[i]))
      })

      list(
            lines = st_sf(mag = spokes$mag, geometry = st_sfc(lines, crs = 4326)),
            points = st_sf(mag = spokes$mag, geometry = st_sfc(points, crs = 4326))
      )
}

aggregate <- function(x, fact, ...){
      if(fact == 1) return(x)
      terra::aggregate(x, fact = fact, ...)
}

tooltip <- function(x, ...) bslib::tooltip(x, ..., options = list(delay = list(show = 2000, hide = 100)))


# UI ---------------------------------------------------
ui <- page_sidebar(

      tags$style(HTML("
  /* Only tiles belonging to the provider tile layer */
  #map .leaflet-layer.basemap-gray img.leaflet-tile {
    filter: grayscale(100%) brightness(100%) opacity(35%);
    -webkit-filter: grayscale(100%) brightness(100%) opacity(35%);
  }
")),

      tags$style(HTML("
  /* Notices come from two mechanisms that cannot share markup: inline divs in
     the sidebar, and validate() messages, which Shiny renders as plain text in
     place of an output. Styling them together is what makes them read as one
     kind of thing.

     Pitched between the two extremes it replaced: colored enough to register as
     something to read rather than blending into body text, but without the
     solid fill of an alert box. Most of these report an unmet precondition -- no polygon
     uploaded, no focal site set -- which is a normal state on the way to using
     the tool, and a full-width warning flag overstates it. */
  .ts-notice,
  .shiny-output-error-validation {
    display: block;
    color: #9a3412;
    background-color: transparent;
    border: none;
    border-left: 3px solid #c2410c;
    border-radius: 0;
    padding: 0.35rem 0 0.35rem 0.7rem;
    margin: 0.5rem 0;
    font-size: 0.85rem;
    line-height: 1.4;
    font-style: normal;
    text-align: left;
  }
")),

      tags$style(HTML("
  /* Raster legend: title beside the color bar rather than stacked above it.
     A multi-line title over a narrow gradient leaves a lot of dead width, and
     the title is the part that wraps. Targets .info.legend, which addLegend()
     produces; circleLegend() builds its own .leaflet-control container and is
     unaffected. */
  /* Gradient legends only: a discrete legend's swatch rows are themselves the
     container's children, so laying them out in a row strings every category
     across the screen. Those opt out via .legend-discrete. */
  #map .info.legend:not(.legend-discrete) {
    display: flex;
    flex-direction: row;
    align-items: flex-start;
    gap: 10px;
  }
  #map .info.legend > div:first-child {
    max-width: 11em;
    min-width: 6em;
    line-height: 1.2;
    text-align: left;
  }
")),

      tags$style(HTML("
  /* Remove padding on the main content wrapper(s) */
  .bslib-page-sidebar main,
  .bslib-page-sidebar .main {
    padding: 0 !important;
  }
")),

      tags$head(
            tags$style(HTML("
              /* Make the tooltip wider */
              .tooltip-inner {
                  max-width: 350px !important;
                  text-align: left;
                  text-transform: none !important;
                  line-height: 1.2 !important;
                  padding: 8px 12px;
                  font-size: 0.85rem;
              }

              /* Subtle styling for definition list descriptions */
              .tooltip-inner dd {
                  color: #aaa;
                  margin-left: 0;
                  margin-bottom: 0.5em;
              }

              /* Optional: make terms stand out more */
              .tooltip-inner dt {
                  color: #fff;
                  font-weight: 600;
              }
      "))
      ),

      theme = bs_theme(version = 5, bootswatch = "flatly"),

      fillable = TRUE,   # lets the page body fill the browser height

      title = "TreeShift [beta]",


      ## LEFT SIDEBAR ------------------------------------

      sidebar = sidebar(
            open = "always",
            width = 320,

            accordion(
                  id = "sidebar_acc",
                  open = "Info",
                  multiple = FALSE,

                  accordion_panel(
                        "Info",
                        helpText("This tool interpolates FIA variables across a spatial grid based on geographic proximity
                          and climatic similarity to FIA plots, optionally accounting for predicted climate change.
                          Use it to map forest properties like occurrence, abundance, and basal area of individual species,
                          forest types, and the tree community as a whole, and to assess potential future shifts in
                          the distributions of these variables."),
                        br(),
                        helpText("Here on the left, you'll find core controls for setting model parameters and
                                 choosing display variables. The righthand sidebar houses modules with their own
                                 dedicated controls and plots. Throughout, hover your cursor over controls for explanatory pop-ups."),
                        br(),
                        helpText("NOTE: This tool is a prototype and is still under development.
                          Please contact matthew.kling@umontana.edu with questions, suggestions, or bug reports.")
                  ),

                  accordion_panel(
                        "FIA variable",
                        ## Labels lowercased to match the other parameter menus;
                        ## values are unchanged, since switch() statements
                        ## throughout the server match on them.
                        selectInput("track", "Entity",
                                    choices = c("species"      = "Species",
                                                "forest type"  = "Forest type",
                                                "forest group" = "Forest group",
                                                "all trees"    = "All trees"),
                                    selected = "Species") %>%
                              tooltip(HTML("Choose the focal entity for the model.<br><br>
                                          <dl style='margin:0;'>
                                              <dt>Species</dt> <dd>Model the presence or abundance of a single focal tree species.</dd>
                                              <dt>Forest type</dt> <dd>Model the probability that a site is classified as a focal FIA forest type (a fine-grained classification).</dd>
                                              <dt>Forest group</dt> <dd>Model the probability that a site is classified within a focal FIA forest type group (a coarser classification).</dd>
                                              <dt>All trees</dt> <dd>Model community-wide forest properties (any-species presence or total basal area).</dd>
                                          </dl>")),
                        ## How MANY entities, orthogonal to which quantity is
                        ## measured on each. Unlabelled, so it reads as a
                        ## qualifier on the entity dropdown above rather than as
                        ## a separate axis. Option text is track-specific ("single
                        ## species" vs "single forest type") and set by an
                        ## observer; the underlying values stay single/multi.
                        ##
                        ## Hidden for All trees, which already aggregates the
                        ## community into a single response BEFORE interpolation
                        ## -- the opposite operation to interpolating each entity
                        ## and reducing afterward.
                        conditionalPanel(
                              condition = "input.track != 'All trees'",
                              selectInput("entity_scope", NULL,
                                          choices = c("single species" = "single",
                                                      "community composition" = "multi"),
                                          selected = "single") %>%
                                    tooltip(HTML("How many entities (species or forest types) to model at once.<br><br>
                                                <dl style='margin:0;'>
                                                    <dt>Single</dt> <dd>Model one focal species, forest type, or group.</dd>
                                                    <dt>Community composition</dt> <dd>Model every entity in the region at once, reducing results to a community-level summary.</dd>
                                                </dl>
                                                The life stage and variable settings apply either way."))
                        ),
                        conditionalPanel(
                              condition = "input.track == 'Species' && input.entity_scope != 'multi'",
                              selectizeInput("species", NULL, choices = spp, selected = "Pinus edulis") %>%
                                    tooltip("Select a focal species by scrolling or typing")
                        ),
                        conditionalPanel(
                              condition = "input.track == 'Forest type' && input.entity_scope != 'multi'",
                              selectizeInput("forest_type", NULL,
                                             choices = forest_type_choices,
                                             selected = forest_type_choices[[1]]) %>%
                                    tooltip("Select a focal FIA forest type by scrolling or typing. Labels show the FORTYPCD followed by the type name.")
                        ),
                        conditionalPanel(
                              condition = "input.track == 'Forest group' && input.entity_scope != 'multi'",
                              selectizeInput("forest_group", NULL,
                                             choices = forest_group_choices,
                                             selected = forest_group_choices[[1]]) %>%
                                    tooltip("Select a focal FIA forest type group by scrolling or typing.")
                        ),
                        ## Life stage closes the entity group -- kind, scope,
                        ## which one, which individuals -- all four unlabelled so
                        ## they read as a single cluster under one heading. The
                        ## scope row poses a question the picker row answers
                        ## ("single species" / "which species"), so stage trails
                        ## rather than splitting that pair.
                        ##
                        ## That placement inverts the old dependency: stage used
                        ## to be filtered by metric, and now metric is filtered by
                        ## stage (metric_choices takes it). A control whose
                        ## options are set by a control below it reads badly, and
                        ## the inverted direction is the more truthful one --
                        ## "seedling basal area" is unavailable because FIA does
                        ## not break basal area out by size class, not because
                        ## basal area lacks a stage.
                        ##
                        ## Forest type and group have no stage: their entity is a
                        ## plot classification, and multi_values() builds an
                        ## indicator matrix that ignores the variable entirely.
                        conditionalPanel(
                              condition = "input.track == 'Species' || input.track == 'All trees'",
                              selectInput("stage", NULL,
                                          c("seedlings only"  = "seedling",
                                            "saplings only"   = "sapling",
                                            "adults only"     = "adult",
                                            "all life stages" = "any"),
                                          selected = "any") %>%
                                    tooltip(HTML("Which tree size class to model. Variables that FIA does not break out by size class &mdash; basal area, and its proportion &mdash; are only offered under 'all life stages'.<br><br>
                                                <dl style='margin:0;'>
                                                    <dt>Seedlings</dt> <dd>&lt; 1 in. DBH/DRC.</dd>
                                                    <dt>Saplings</dt> <dd>1.0&ndash;4.9 in. DBH/DRC.</dd>
                                                    <dt>Adults</dt> <dd>&ge; 5.0 in. DBH/DRC.</dd>
                                                    <dt>All life stages</dt> <dd>Combined.</dd>
                                                </dl>"))
                        ),
                        ## Single-entity pickers and the community measure share
                        ## one slot: in All mode the reduction IS the selection.
                        ## Hidden rather than removed, so input$species keeps its
                        ## value while All mode is active.
                        ## Number of categories drawn for the categorical
                        ## measures. 143 entities is far past what any palette or
                        ## legend can carry, so all but the most frequent N are
                        ## lumped into a single "other" class. The right N differs
                        ## a lot between a polygon-scoped domain and a
                        ## region-wide one, so it is a control rather than a
                        ## constant.
                        conditionalPanel(
                              condition = paste0(
                                    "input.track != 'All trees' && input.entity_scope == 'multi' && ",
                                    "['", paste(categorical_measures, collapse = "','"),
                                    "'].indexOf(input.comm_measure) >= 0"),
                              numericInput("comm_top_n", "Categories shown", value = 10, min = 2, max = 25, step = 1) %>%
                                    tooltip("How many entities get their own color. The rest are pooled into 'other', ranked by how much of the mapped area each one covers.")
                        ),

                        ## Timeframe sits above the variable because it is the
                        ## first cut on the community measures below: "state or
                        ## change" narrows the list before the user reads it,
                        ## rather than the list being long and mixed and the
                        ## timeframe control quietly not applying to half of it.
                        ##
                        ## "change" covers both a difference of a static measure
                        ## (change in richness) and the paired measures, which are
                        ## change quantities in their own right and not deltas of
                        ## anything -- net richness change of zero is entirely
                        ## compatible with complete species replacement.
                        selectInput("time", "Timeframe",
                                    c("historical", "future", "change" = "delta")) %>%
                              tooltip("Choose whether to display historical (2000-2025) or future (scenario selected under 'Climate data') model predictions, or the change between them.
                                      Only relevant for `FIA variable` and `ESS` raster variables."),
                        ## Declared from metric_choices() rather than a literal
                        ## list. The syncing observer uses ignoreInit, so a
                        ## hard-coded list is the startup state -- and it held
                        ## "abundance", which is not valid at the default life
                        ## stage of "any". Selecting it made resolve_var() return
                        ## NULL and the app silently do nothing.
                        selectInput("metric", "Variable",
                                    choices = metric_choices("Species", "any")) %>%
                              tooltip(HTML("Choose the quantity to model.<br><br>
                                          <dl style='margin:0;'>
                                              <dt>Occurrence</dt> <dd>Probability the focal entity is present within a single FIA plot footprint (~0.07 ha) randomly located within the grid cell.</dd>
                                              <dt>Abundance</dt> <dd>Population density of the selected life stage, in stems per acre.</dd>
                                              <dt>Basal area</dt> <dd>Live tree and sapling basal area, in ft&sup2; per acre.</dd>
                                              <dt>Species richness</dt> <dd>Number of tree species present within a single FIA plot footprint (~0.07 ha).</dd>
                                              <dt>Proportion basal area</dt> <dd>The focal species' share of total stand basal area, multiplied by the probability that stand basal area is nonzero.</dd>
                                              <dt>Regeneration fraction</dt> <dd>Probability that the selected regeneration stage is present relative to adults, conditional on the species occurring on the plot. Calculated as p(regen) / (p(regen) + p(adult)). Adults only = 0; regen only = 1; adults and regen together = intermediate value (not precisely 0.5 due to differing sample areas).</dd>
                                          </dl>")),
                        ## The community reduction pairs with the variable, not
                        ## with the entity: "occurrence + richness" and "basal
                        ## area + dominance" jointly define the mapped quantity,
                        ## the same way the variable alone does in single-entity
                        ## mode. It sits below the variable because
                        ## comm_measure_choices() is filtered by it.
                        conditionalPanel(
                              condition = "input.track != 'All trees' && input.entity_scope == 'multi'",
                              selectInput("comm_measure", NULL,
                                          choices = comm_measure_choices("Species", "presence probability"),
                                          selected = "richness") %>%
                                    tooltip(HTML("How the per-entity (species or forest types) surfaces are combined into one map. Only combinations that are well defined for the current entity type and variable are listed.<br><br>
                                                <dl style='margin:0;'>
                                                    <dt>Richness</dt> <dd>Expected number of entities present, summed across per-entity occurrence probabilities. No threshold: a cell with 0.5 probability for each of two species has expected richness 1.</dd>
                                                    <dt>Diversity</dt> <dd>Effective number of entities (Hill N1 = exp of Shannon entropy), in the same units as richness. Low where one entity dominates, however many rare ones are present.</dd>
                                                    <dt>Evenness</dt> <dd>Pielou&#39;s J (Shannon entropy relative to its maximum given the species pool). 0 where a single species dominates, with higher values for more balanced communities.</dd>
                                                    <dt>Dominance</dt> <dd>Share of the community held by its single largest member (0&ndash;1).</dd>
                                                    <dt>Dominant entity</dt> <dd>Which entity is has the highest value at each cell. Drawn as discrete classes; all but the most frequent few are pooled into &quot;other&quot;. Most intuitive over an abundance variable &mdash; over occurrence it shows the most widely probable entity.</dd>
                                                </dl>
                                                The measures below compare the historical and future communities directly (they are not deltas of the measures above).<br>
                                                <dl style='margin:0;'>
                                                    <dt>Turnover</dt> <dd>Bray&ndash;Curtis dissimilarity between the historical and future communities. 0 is identical composition, 1 is no overlap.</dd>
                                                    <dt>Expected gains / losses</dt> <dd>How much arrives and how much departs, in the variable&#39;s own units, counted separately. Asymmetric on purpose &mdash; what arrives is a different question from what leaves.</dd>
                                                    <dt>Colonization demand</dt> <dd>Gains as a share of the future community: the fraction of what will be suitable there that is not there now.</dd>
                                                    <dt>Loss fraction</dt> <dd>The opposite: losses as a share of the community that is there now. Read alongside turnover, the two fractions show which side is driving it &mdash; high turnover with a low loss fraction is accumulation rather than replacement.</dd>
                                                    <dt>Biggest winner / loser</dt> <dd>Which single entity gains or loses the most at each cell, drawn as discrete classes. Says which entity is moving, where the measures above say how much. Use the Community composition panel to see the full ranking behind any one cell.</dd>
                                                </dl>")),
                              uiOutput("comm_note"),
                              uiOutput("comm_cost_note")
                        ),
                        ## Support floor for ratio-valued layers. Shown only when
                        ## the selected layer is a ratio -- the community
                        ## measures other than richness, and the regeneration
                        ## fractions. Richness IS the denominator, so masking it
                        ## by itself would be circular.
                        conditionalPanel(
                              condition = "(input.track != 'All trees' && input.entity_scope == 'multi' && input.comm_measure != 'richness' && input.comm_measure != 'expected gains' && input.comm_measure != 'expected losses') || input.metric == 'regeneration fraction'",
                              numericInput("support_floor", "Minimum expected species", value = 0.5, min = 0, step = 0.5) %>%
                                    tooltip(HTML("Blank out cells where there is too little predicted mass present for the ratio to mean anything.<br><br>
                                                Ratios stay well behaved as their denominator approaches zero but stop being informative: evenness and the regeneration fractions take their most extreme values where the community or species is barely present.<br><br>
                                                The label shows which quantity the threshold applies to &mdash; expected richness for occurrence variables, the summed response for abundance ones, and the species' own occurrence probability for the regeneration fractions. Set 0 for no masking.<br><br>
                                                This also rescales the color palette, since the extreme values from near-empty cells otherwise compress the rest of the map into a narrow band."))
                        ),

                        ## Canonical variable resolved from (track, metric, stage);
                        ## hidden because all downstream logic switches on input$fia_var.
                        tags$div(style = "display: none;",
                                 selectInput("fia_var", NULL,
                                             c("presence probability", "seedling presence",
                                               "sapling presence", "tree presence",
                                               "proportion basal area", "total basal area",
                                               "number of trees", "number of saplings",
                                               "number of seedlings",
                                               "species richness", "seedling species richness",
                                               "sapling species richness", "adult species richness",
                                               "seedling fraction", "sapling fraction")))
                  ),

                  ## Three nested spatial extents, one shared resolution.
                  ##
                  ## The model bounding box, the focal polygon, and the focal
                  ## site are all "where", differing only in dimension: an area
                  ## computed over, an area summarized within, a point queried at.
                  ## Splitting them across two panels (or, before that, across
                  ## two sidebars) made near-identical things look unrelated.
                  ##
                  ## Grain sits above all three rather than with the domain: it
                  ## is the resolution everything here is evaluated at. The
                  ## polygon density plot and the climate-models summarize-over
                  ## sample the display raster at this grain, and the focal-site
                  ## similarity surface is on the same grid.
                  ##
                  ## Note this panel deliberately mixes controls that invalidate
                  ## aim() (bbox, grain) with ones that do not (polygon, site).
                  ## The old left-panel rule kept those apart, but that rule had
                  ## no home for shared-but-not-model controls. If the deferred
                  ## update button lands, it should key off an explicit list of
                  ## model-parameter inputs rather than off panel membership --
                  ## otherwise every future placement decision silently becomes
                  ## an invalidation decision.
                  accordion_panel(
                        "Study area",
                        helpText("Where the model runs, and the places you compare against it: the model bounding box,
                                 an optional focal polygon, and an optional focal site. All three are evaluated at the
                                 grain size set below."), hr(),

                        sliderInput("clim_agg", "Raster grain size", 1, 10, value = 10, step = 1, post = " km") %>%
                              tooltip("Cell size of the model grid, and of every summary drawn from it. Larger values run faster but are less precise."),

                        hr(),

                        tags$strong("Model bounding box"),
                        helpText("The rectangulare region the model is computed over. A smaller box runs faster."),
                        sliderInput("bbox_x", "Longitude range",
                                    min = clim_bbox["xmin"], max = clim_bbox["xmax"],
                                    value = c(-116, -100)),
                        sliderInput("bbox_y", "Latitude range",
                                    min = clim_bbox["ymin"], max = clim_bbox["ymax"],
                                    value = c(29, 43)),
                        layout_column_wrap(
                              width = 1/2,
                              uiOutput("snap_focal_btn"),
                              actionButton("window_domain", "Snap to current view")
                        ),

                        hr(),

                        tags$strong("Focal polygon"),
                        helpText("A spatial polygon of your area of interest. Optional."),
                        uiOutput("polygon_status"),
                        fileInput("polygon_upload", NULL,
                                  multiple = TRUE,
                                  accept = c(".shp", ".dbf", ".shx", ".prj", ".kml", ".kmz",
                                             ".gpkg", ".geojson"),
                                  placeholder = "Shapefile, KML, GPKG, or GeoJSON") %>%
                              tooltip("Upload a polygon layer. For shapefiles, select all component files (.shp, .dbf, .shx, .prj)."),
                        actionButton("polygon_clear", "Clear") %>%
                              tooltip("Remove the polygon from the map and plots."),

                        hr(),

                        tags$strong("Focal site"),
                        helpText("A single location to query. Click the map to set one. Optional."),
                        uiOutput("site_status"),
                        ## The "highlight on map" button that used to sit here has
                        ## moved to the "How predictions are made" panel. This
                        ## group is for SETTING the site; showing which plots
                        ## inform it is that panel's subject.
                        actionButton("site_clear", "Clear") %>%
                              tooltip("Remove the focal site from the map and plots.")
                  ),

                  accordion_panel(
                        "Climate data",
                        helpText("Select climate variables and future scenarios."), hr(),
                        selectInput("clim_vars", "Climate variables", clim_vars(), selected = clim_vars_default, multiple = TRUE) %>%
                              tooltip(HTML("Choose one or more variables for measuring multivariate climate similarity:<br><br>
                                          <dl style='margin:0;'>
                                              <dt>CWD</dt> <dd>Climatic water deficit</dd>
                                              <dt>AET</dt> <dd>Actual evapotranspiration</dd>
                                              <dt>PPT</dt> <dd>Total annual precipitation (log-transformed)</dd>
                                              <dt>TMEAN</dt> <dd>Mean annual temperature</dd>
                                              <dt>TMINCM</dt> <dd>Avg min temp of coldest month</dd>
                                              <dt>TMAXWM</dt> <dd>Avg max temp of warmest month</dd>
                                          </dl>")),
                        selectInput("scenario_ssp", "Emissions pathway",
                                    choices = ssp_choices,
                                    selected = default_ssp) %>%
                              tooltip("Future emissions scenario (Shared Socioeconomic Pathway). Higher numbers correspond to higher emissions."),
                        selectInput("scenario_period", "Future time period",
                                    choices = period_choices,
                                    selected = default_period) %>%
                              tooltip("Future time period."),
                        selectizeInput("scenario_model", "Climate model(s)",
                                       choices = model_choices,
                                       selected = default_model,
                                       multiple = TRUE,
                                       options = list(plugins = list("remove_button"))) %>%
                              tooltip("Select one or more global climate models (GCMs). When multiple are selected, climate layers are averaged across the model ensemble.
                                      See the 'Climate models' panel in the right sidebar for help comparing GCMs and choosing among them.")
                  ),

                  accordion_panel(
                        "Interpolation",
                        helpText("Predictions for each grid cell are built from nearby FIA plots with similar climate. Plots are weighted
                                 more heavily the closer they are, in both senses, and these two settings control how quickly that
                                 falloff happens. They are standrd deviations of Gaussian weighting kernels.
                                 Plots more than three times either setting away are not used at all."), hr(),
                        sliderInput("theta_geog", "Geographic bandwidth (km)", 5, 100, value = 30) %>%
                              tooltip("Controls how plot weight decays with geographic distance from a grid cell.
                                      This distance also defines two of the areas used in some of the righthand modules: the 'focal site neighborhood' is everything within three bandwidths of the site you clicked, and the 'buffered polygon' is your polygon widened by the same amount."),
                        sliderInput("theta_clim", "Climate bandwidth (z)", .05, .5, value = .1) %>%
                              tooltip("Measured in standard deviations of climate across all FIA plots in the contiguous US, so 0.1 means plots within a tenth of a standard deviation count almost fully.")
                  ),

                  accordion_panel(
                        "Map layers",
                        helpText("Choose which model results are drawn on the map."), hr(),
                        selectInput("stat", "Raster variable", stat_choices_base) %>%
                              tooltip(HTML("<dl style='margin:0;'>
                                              <dt>FIA variable</dt> <dd>Predicted values of the variable selected above.</dd>
                                              <dt>Effective sample size</dt> <dd>How many FIA plots each cell's prediction effectively rests on. Nearby, climatically similar plots count close to fully and distant ones count for little, so this is usually well below the raw plot count. Low values mean a prediction is built on thin evidence.</dd>
                                              <dt>Focal site similarity</dt> <dd>How climatically similar each cell is to the site you clicked on the map. Which climates are compared depends on the analog direction set in the 'Focal site analogs' panel at right.</dd>
                                              <dt>Connectivity</dt> <dd>Circuit-based connectivity outputs between historical and future surfaces (available when the 'Lanscape connectivity' module is enabled); choose the specific connectivity variable below</dd>
                                          </dl>")),
                        conditionalPanel(
                              condition = "input.stat == 'connectivity'",
                              selectInput("cs_stat", "Connectivity variable", cs_stat_names) %>%
                                    tooltip(HTML("<dl style='margin:0;'>
                                              <dt>Current flow</dt> <dd>Total current passing through each cell; high values indicate important range shift corridors</dd>
                                              <dt>Current on target</dt> <dd>Current delivered to destination grounds; shows where dispersers successfully reach future-suitable habitat</dd>
                                              <dt>Current loss</dt> <dd>Current lost to non-target grounds; represents dispersers removed via mortality or settlement before reaching destinations</dd>
                                              <dt>Current pinch points</dt> <dd>Power dissipation (I²R); highlights bottlenecks where connectivity is both concentrated and tenuous — priorities for restoration or protection</dd>
                                              <dt>Current direction</dt> <dd>Bearing of net current flow through each cell</dd>
                                              <dt>Voltage</dt> <dd>How connected each cell is to source populations; high values indicate greater accessibility for dispersers</dd>
                                              <dt>Conductance</dt> <dd>Model input: ease of movement across the landscape (inverse of resistance)</dd>
                                              <dt>Loss</dt> <dd>Model input: rate at which dispersers are removed from the migrating pool at each cell</dd>
                                          </dl>"))
                        ),
                        selectInput("vector", "Vector variable", vector_choices_base) %>%
                              tooltip(HTML("Select which point or line layer is displayed:<br><br>
                                          <dl style='margin:0;'>
                                              <dt>Occupied FIA plots</dt> <dd>Plots where the selected species was recorded</dd>
                                              <dt>All FIA plots</dt> <dd>All plots including non-forested plots not surveyed</dd>
                                              <dt>Focal site analogs</dt> <dd>Plots within climatic and geographic range of the selected site. Which climates are matched &mdash; contemporary, reverse or forward &mdash; is set in the 'Focal site analogs' panel in the right sidebar.</dd>
                                              <dt>Current animation</dt> <dd>Animation with particle movement showing direction and magnitude of mean current flow (available when the connectivity module is enabled)</dd>
                                              <dt>Current field</dt> <dd>Vector field showing showing direction and magnitude of mean current flow (available when the connectivity module is enabled)</dd>
                                          </dl>"))
                  ),

                  accordion_panel(
                        "Cartography",
                        helpText("Adjust the styling of the map."), hr(),
                        sliderInput("vector_scale", "Vector point/line size", -1, 1, 0, step = .1) %>%
                              tooltip("Use this slider to adjust the size of FIA plot markers or the length of current lines"),
                        sliderInput("rast_opacity", "Raster opacity", 0, 1, .8, step = .01),
                        selectInput("rast_transform", "Scale transformation",
                                    choices = c("identity", "log10", "sqrt", "uniform")) %>%
                              tooltip("Choose transformation used to visualize raster values. Legend will display transformed values."),
                        ## Choices swap between continuous ramps and qualitative
                        ## sets depending on the layer, rather than showing both
                        ## lists at once and relying on the user to pick the
                        ## applicable one.
                        selectInput("rast_palette", "Color palette",
                                    choices = rast_palettes_continuous,
                                    selected = "default") %>%
                              tooltip("Override the default palette for the current raster variable. When the raster shows the historical or future FIA variable, the FIA plot markers and focal-site scatterplots adopt the same palette. For layers showing entity identities the options become qualitative palettes."),
                        input_switch("show_site", "Show focal site", TRUE),
                        uiOutput("site_color_ui"),
                        input_switch("show_polygon", "Show focal polygon", TRUE),

                        ## Hidden until a polygon exists -- an outline color for
                        ## nothing is just clutter. Rendered server-side rather
                        ## than with a conditionalPanel because the condition is
                        ## "a file has been uploaded", which is not visible to JS.
                        uiOutput("polygon_color_ui")
                  ),

                  accordion_panel(
                        "Export",
                        helpText("Download the raster data displayed on the map as a .tif, the underlying FIA plot data as a .csv,
                                  or the current model settings as a .txt for reproducibility."), hr(),
                        textInput("download_filename", "Filename (without extension)", "fia_aim_result"),
                        downloadButton("download_raster", "Download raster"),
                        br(),
                        downloadButton("download_points", "Download FIA points") %>%
                              tooltip(HTML("Plot-level FIA data for the points currently drawn on the map.<br><br>
                                           <dl style='margin:0;'>
                                               <dt>x, y</dt> <dd>Plot location (publicly released, fuzzed ~1.6 km)</dd>
                                               <dt>entity</dt> <dd>Selected species, forest type, or forest group</dd>
                                               <dt>variable</dt> <dd>Selected FIA variable</dd>
                                               <dt>value</dt> <dd>Observed historical value of that variable</dd>
                                               <dt>occupied</dt> <dd>Whether the entity is present at the plot</dd>
                                               <dt>sampled</dt> <dd>FALSE at plots that were never surveyed</dd>
                                           </dl>
                                           The analog layer adds climate/geographic distance from the focal site and
                                           the analog's kernel weight.")),
                        br(),
                        downloadButton("download_metadata", "Download metadata")
                  )
            )
      ),



      card(
            class = "h-100 w-100 rounded-0 overflow-hidden border-0",
            layout_sidebar(
                  fillable = TRUE,


                  ## RIGHT SIDEBAR ---------------------------------------

                  sidebar = sidebar(
                        # title = "Focal site analogs",
                        width = "40%",
                        position = "right",
                        open = TRUE,

                        accordion(
                              open = FALSE,
                              multiple = FALSE,

                              ## First panel: the map-level community measures
                              ## reduce ACROSS entities, which is exactly the
                              ## information this panel restores. A turnover or
                              ## richness surface says how much changes; this
                              ## says which entities changed, and is the only
                              ## view in the app where entity identity survives
                              ## the reduction.
                              accordion_panel(
                                    "Community composition",
                                    helpText("When you select the 'community composition' FIA variable at left,
                                             the map summarizes across species or forest types into one summary number per cell.
                                             This module lets you unpack the individual species behind those summaries. The two settings
                                             below apply to all three plots."), hr(),

                                    ## Both controls apply to every plot below, so
                                    ## they sit together above the rule rather
                                    ## than beside whichever plot they happen to
                                    ## precede. "Entities shown" is the count
                                    ## NAMED rather than pooled into "other", and
                                    ## it feeds all three.
                                    ## The full set: unlike the other modules, a
                                    ## single cell is meaningful here -- reading
                                    ## the composition at one place is the point
                                    ## of drilling in.
                                    selectInput("comm_scope", "Summarize over", area_names) %>%
                                          tooltip("Which part of the map these plots describe. All but the last average over many cells; 'Focal site' reads just the cell you clicked, so it shows what is at one place rather than a regional average."),
                                    sliderInput("comm_rank_n", "Entities shown", 5, 40, value = 15, step = 1) %>%
                                          tooltip("How many species get named individually. The rest are combined into a single 'other' group rather than left out, so the totals stay right."),
                                    hr(),

                                    ## Leads the panel: it is the only view here
                                    ## that shows the community as a whole. Band
                                    ## thickness is one entity and total column
                                    ## height is the community sum, so a
                                    ## community that reshuffles without changing
                                    ## size looks different from one that shrinks.
                                    helpText("Current and predicted future forest composition. Each band is one species; the
                                             height of the whole column is the community."),
                                    plotOutput("comm_composition") %>%
                                          tooltip("Each band runs from where a species is today, on the left, to where it is expected to be in the future, on the right. Labels sit on whichever side the band is thicker: on the left for species that decline, on the right for those that grow."),
                                    hr(),

                                    ## Second: the same numbers spread out, so
                                    ## rare-but-gaining stays distinct from
                                    ## abundant-but-declining. Absolute and
                                    ## proportional change rank differently and
                                    ## neither is wrong; the scatter shows both
                                    ## at once.
                                    helpText("A different view of the same data, showing whether gains are predicted in species that are already common or in ones
                                             that are currently rare. Rare and gaining sits at the top left; common and declining sits
                                             at the bottom right. Colored points are the species named above, gray are the rest."),
                                    plotOutput("comm_scatter"), hr(),

                                    ## Last: change without levels, for the cases
                                    ## where the ordering itself is the answer.
                                    ## Its "Rank by" control stays here because,
                                    ## unlike the two above, it applies to this
                                    ## plot alone.
                                    helpText("Species with the largest projected gains and losses."),
                                    selectInput("comm_rank_by", "Rank by",
                                                c("absolute change" = "absolute",
                                                  "proportional change" = "proportional")) %>%
                                          tooltip("Absolute change ranks by how much is gained or lost outright, so common species tend to lead. Proportional change measures that against how much is there now, so a rare species that doubles outranks a common one that shifts a little. Neither ordering is more correct than the other."),
                                    plotOutput("comm_ranking")
                              ),

                              ## Named for its contents rather than its input:
                              ## both plots are views of site_analogs(), the
                              ## analog SET of the focal site. The site itself is
                              ## now set under "Study area" in the left sidebar.
                              ##
                              ## The composition plot that used to live here has
                              ## been replaced by a stack-based one in the
                              ## Community composition panel: it was the only view
                              ## here where entity identity survived, and it
                              ## belonged with the other cross-entity displays.
                              ## Renamed from "Study area": that name now belongs
                              ## to the left panel that owns the polygon. This
                              ## panel is named for what it renders, matching the
                              ## rule that right-sidebar panels are local to one
                              ## output.
                              accordion_panel(
                                    "Value distributions",
                                    helpText("This module summarizes how the predicted values on the map differ across geographic areas.
                                             It shows one curve (or bar stack, for categorical variables) for the whole bounding box,
                                             and one for each area you have set: your focal polygon, that polygon
                                             plus a surrounding buffer, and the neighborhood around your focal site."), hr(),

                                    plotOutput("polygon_density"),
                                    selectInput("density_transform", "X-axis scale",
                                                choices = c("identity", "uniform", "sqrt", "log10"),
                                                selected = "identity") %>%
                                          tooltip("Optionally transform the x-axis for skewed variables. 'uniform' rank-transforms the values across the bounding box, spreading skewed distributions evenly.")
                              ),

                              accordion_panel(
                                    "Climate relationships",
                                    helpText("This plot shows how the mapped variable changes along a climate gradient. The curve is a smoothed
                                             relationship fitted to either the FIA plot observations or the modeled grid cells; the
                                             density histogram above it shows how much data sits at each point along the gradient.
                                             For categorical variables, this becomes one distribution per species instead, showing where along the
                                             gradient each one comes out on top."), hr(),

                                    ## Two columns of two, paired by what they
                                    ## control: the left column defines the x
                                    ## axis, the right column defines which data
                                    ## the curve is fitted to.
                                    layout_column_wrap(
                                          width = 1/2,
                                          selectInput("cr_v1", "X-axis variable", display_vars()) %>%
                                                tooltip("The variable the response is plotted against. It need not be one the model was fitted on; the plot says so when it isn't."),
                                          selectInput("cr_data", "Dataset", c("modeled grid cells", "FIA plots")) %>%
                                                tooltip("Whether to fit the curve to the model's predictions across the grid, or to the raw FIA plot observations."),
                                          selectInput("cr_clim", "Climate timeframe", c("historical", "future")) %>%
                                                tooltip("Which climate the x-axis is drawn from. Use 'future' to pair a future prediction with the future climate that produced it. FIA plot observations are historical, so pairing them with future climate relates an observed response to a climate the plot has not experienced."),
                                          ## Same scope vocabulary as the evaluation and composition
                                          ## modules, from the same named vector, so an area means the
                                          ## same thing wherever it is offered.
                                          selectInput("cr_scope", "Area", eval_scope_names) %>%
                                                tooltip("Which part of the map the relationship is fitted to. Narrowing it shows how the variable behaves in your area specifically, which can differ from the regional pattern.")
                                    ),


                                    plotOutput("climate_response", height = "600px")
                              ),

                              accordion_panel(
                                    "Model evaluation",
                                    helpText("This module uses 'leave-one-out cross validation' to quantify the accuracy of the model's present-day predictions.
                                             Each FIA plot is held out in turn, predicted from the others, and compared against what was actually measured there.
                                             The first plot below shows accuracy for the currently-displayed results, while the second section shows how accuracy
                                             varies over a range of alternative model parameters you could select.
                                             You can choose the spatial area ('evaluation scope') that accuracy is measured over."),
                                    br(),
                                    helpText("Note that these accuracy estimates and optimal parameter recommendations DO NOT apply to future predictions."),
                                    br(),

                                    selectInput("eval_scope", "Evaluation scope",
                                                choices = eval_scope_names,
                                                selected = "domain") %>%
                                          tooltip(HTML("Which held-out plots the accuracy metrics are computed from. The model itself is unchanged.<br><br>
                                                       <dl style='margin:0;'>
                                                           <dt>Model bounding box</dt> <dd>Every plot in the analog pool.</dd>
                                                           <dt>Focal polygon</dt> <dd>Plots inside the uploaded polygon. Often few; check the reported n.</dd>
                                                           <dt>Buffered polygon</dt> <dd>Plots within 3 geographic bandwidths of the polygon &mdash; the neighborhood predictions there are actually drawn from. Usually the better choice for a small unit.</dd>
                                                           <dt>Focal site neighborhood</dt> <dd>Plots within the geographic cutoff of the clicked site.</dd>
                                                       </dl>
                                                       Scopes whose input is not yet set under 'Study area' in the left sidebar will say so when selected.")),

                                    hr(),

                                    helpText("This plot compares model predictions against what was actually measured at each plot, using your current
                                             settings. Each red point is a group of plots with similar predicted values, plotted against what
                                             was really found there. Points on the dashed 1:1 line mean the predictions are unbiased, while
                                             deviations from the 1:1 line indicate bias; the shaded band is the IQR of observed values. Marks along the
                                             bottom show where the predictions actually fall. When you narrow the spatial scope, a gray line
                                             shows the same curve for plots in the model domain that are outside your scope."), br(),
                                    plotOutput("eval_scatter"),
                                    uiOutput("eval_scatter_stats"), hr(),

                                    helpText("This tool helps you explore how different parameter values affect predictive accuracy.
                                             For bandwidth parameters, it tests five values bracketing the value selected under 'Analog bandwidths'.
                                             For climate variables, it tests all combinations of the variables selected under 'Climate data'.
                                             The scatterplot above always reflects your current settings; the heatmap below shows the last completed
                                             run, and is faded and flagged if those settings have since changed."), br(),

                                    layout_column_wrap(
                                          width = 1/2,
                                          checkboxGroupInput("eval_params", NULL,
                                                             choiceNames = c("Test geographic bandwidths", "Test climate bandwidths", "Test climate variables"),
                                                             choiceValues = c("geog", "clim", "vars")),
                                          ## Choices swap with the response type; declared with the
                                          ## binary set because that matches the default variable.
                                          selectInput("eval_stat", "Evaluation statistic", c("AUC", "TSS", "Brier")),
                                          sliderInput("eval_n_vars", "Max n climate vars", 1, 6, value = c(1, 2), step = 1),
                                          actionButton("optimal_params", "Apply optimal parameters"),
                                          actionButton("run_eval", "Run evaluations")
                                    ),

                                    ## Directly above the heatmap, which is the
                                    ## only thing it describes. Sitting under the
                                    ## scope selector, it read as a caution about
                                    ## the scatterplot too -- and that plot always
                                    ## reflects current settings, so it is never
                                    ## the stale one.
                                    uiOutput("eval_stale"),
                                    plotOutput("eval_heatmap", height = "600px")
                              ),

                              accordion_panel(
                                    "Climate models",
                                    helpText("This module helps you choose among, and contextualize, climate models.
                                             It can take a couple minutes to run.
                                             For every available emissions scenario and GCM, it plots the projected climate change
                                             (future minus historical) as a trajectory across time periods, so
                                             you can see which models are hot/dry, warm/wet, and so on, and how
                                             the spread fans out over eras and emissions scenarios. Models in your
                                             current ensemble selection (set in the left sidebar) are drawn at full
                                             opacity, others faded. The x and y variables can be changed
                                             without re-clicking 'Update'."), hr(),
                                    layout_column_wrap(
                                          width = 1/2,
                                          selectInput("model_delta_region", "Summarize over",
                                                      choices = area_names,
                                                      selected = "domain") %>%
                                                tooltip("Region over which deltas are summarized. 'Model bounding box'
                                                        averages over the current bounding box; 'Focal site' samples
                                                        at the clicked site; 'Focal polygon' averages over the uploaded
                                                        polygon."),
                                          div(style = "padding-top: 1.9rem;",
                                              actionButton("model_delta_update", "Update chart") %>%
                                                    tooltip("Recompute the chart for the current axis variables and
                                                            summary region. The chart does not update automatically,
                                                            since computing deltas across the full model ensemble is
                                                            slow."))
                                    ),
                                    plotOutput("model_delta_plot"),
                                    layout_column_wrap(
                                          width = 1/2,
                                          selectInput("model_delta_x", "X variable",
                                                      choices = clim_vars()) %>%
                                                tooltip("Climate variable for the x-axis (delta = future - historical)."),
                                          selectInput("model_delta_y", "Y variable",
                                                      choices = clim_vars()) %>%
                                                tooltip("Climate variable for the y-axis (delta = future - historical).")
                                    )
                              ),

                              accordion_panel(
                                    "Landscape connectivity",
                                    helpText("This module lets you model circuit-based landscape connectivity between historical (source) and future (target) FIA surfaces.
                                             This approach can identify key movement corridors and pinch points for large-scale climate connectivity.
                                             Click 'enable' to activate options.
                                             To view connectivity results, select 'connectivity' as the raster variable under 'Map layers' at left.
                                             Julia compiler startup may take a couple minutes the first time you run a connectivity result, but should be fast thereafter."), hr(),
                                    uiOutput("cs_note"),
                                    conditionalPanel(
                                          condition = "input.track == 'All trees' || input.entity_scope != 'multi'",
                                          input_switch("cs_enable", "Enable connectivity", FALSE) %>%
                                                tooltip("Adds the connectivity raster and current animation/field vector layers
                                                        to the 'Map layers' menus in the left sidebar.")
                                    ),
                                    conditionalPanel(
                                          condition = "input.cs_enable && (input.track == 'All trees' || input.entity_scope != 'multi')",
                                          selectInput("cs_sources", "Origin electrodes (sources)", c("grid cells", "FIA plots", "selected site")) %>%
                                                tooltip("Sites where current is injected. If `grid cells` or `FIA plots` is selected, current is injected
                                            in proportion to modeled or observed (respectively) values of the selected FIA variable.
                                            If `selected site` is chosen, current is injected only at the focal site."),
                                          selectInput("cs_targets", "Destination electrodes (target grounds)", c("grid cells", "FIA plots", "selected site")) %>%
                                                tooltip("Grounding sites that 'pull' current and allow it to exit the circuit.
                                            If `grid cells` or `FIA plots` is selected, ground conductance is proportional to modeled values
                                            of the selected FIA variable. If `selected site` is chosen, the focal site."),
                                          selectInput("cs_conductance_layer", "Horizontal conductance", c("suitability", "wind", "suitability * wind", "uniform")) %>%
                                                tooltip("The inverse of landscape resistance, horizontal conductance determines how easily current can flow across a site.
                                            `Suitability` models conductance proportional to the mean of historical and future FIA raster values.
                                            `Wind` models conductance proportional to the frequency and speed with which wind blows toward targets."),
                                          selectInput("cs_loss_layer", "Loss (non-target grounds)", c("none", "uniform", "suitability" = "mean suitability")) %>%
                                                tooltip("Choose whether and how to model current 'leakage' representing dispersers that are removed from the pool
                                            (via mortality or permanent settlement) before reaching target grounds."),
                                          sliderInput("cs_conductance_scale", "log10 conductance scale", -2, 2, value = 0, step = .1) %>%
                                                tooltip("Relative strength of horizontal conductance, compared to source, target, and loss components."),
                                          sliderInput("cs_loss_scale", "log10 loss scale", -3, 2, value = -1, step = .1) %>%
                                                tooltip("Relative strength of loss conductance, compared to source, target, and horizontal components.
                                                        Ignored if loss is set to 'none.'")
                                    )
                              ),

                              ## Last, and deliberately so. Every other panel here
                              ## reports a result; this one shows the mechanism
                              ## behind them -- which plots inform a prediction,
                              ## and how the two bandwidths weight them. That is
                              ## worth having for anyone asking where the numbers
                              ## come from, but it is not what a manager opens
                              ## first.
                              ##
                              ## Named for its subject rather than its purpose,
                              ## even though the purpose is what the intro text
                              ## explains. The panel is not purely illustrative:
                              ## analog_direction set here also drives the 'focal
                              ## site similarity' raster and the 'focal site
                              ## analogs' vector layer. Matching those option
                              ## names exactly means cross-references from the
                              ## left sidebar are literal rather than approximate.
                              ##
                              ## The analog SET does not depend on the entity --
                              ## the same plots are near and climatically similar
                              ## whichever species is selected -- so only the
                              ## point coloring is entity-specific, and it follows
                              ## display_values() like the map markers do.
                              accordion_panel(
                                    "Focal site analogs",
                                    helpText(
                                          "This module helps you explore where a model prediction at a given location actually comes from:
                                          the 'analog' FIA plots that are geographically close and climatically similar.
                                          When you click the map to set a focal site, these figures show the FIA plots (drawn as points) inform the prediction.
                                          'Contemporary' analogs compare the focal cell's and plots' current climates, while 'reverse' and 'forward' analogs compare climates across time to make future projections.

                                          Reference lines mark the two bandwidths (dashed) and the hard cutoffs at three times each (solid).
                                          Weights, represented by shading in the righthand figure, fall off with both distances.

                                          Use the button below to see the same analog plots on the map.

                                          The analog type set here also governs two map layers: the 'focal site similarity'
                                          raster and the 'focal site analogs' points. Nothing here changes the model
                                          itself."), hr(),
                                    layout_column_wrap(
                                          width = 1/2,
                                          selectInput("analog_direction", NULL, c("Contemporary analogs", "Reverse analogs", "Forward analogs")) %>%
                                                tooltip(HTML("Which climates are being matched:
                                                <ul>
                                                    <li>Contemporary: plots whose present climate resembles the site's present climate. This is what the historical prediction is built from.</li>
                                                    <li>Reverse (inbound): plots whose present climate resembles the site's future climate. This is the basis for future model predictions in the app.</li>
                                                    <li>Forward (outbound): plots whose future climate resembles the site's present climate. Where the site's current conditions are heading; interesting, but not currently used in this app's core modeling machinery.</li>
                                                </ul>")),
                                          ## Kept, and relabeled. The old name
                                          ## described the mechanism rather than
                                          ## what turning it off would show. In a
                                          ## panel about how predictions are made,
                                          ## seeing which nearby plots the climate
                                          ## kernel rules out is half the point --
                                          ## otherwise the cutoff is a line with
                                          ## nothing outside it.
                                          input_switch("filter_clim", "Hide plots beyond the climate cutoff", TRUE) %>%
                                                tooltip("Turn this off to also show the plots that are close by but too climatically different to count. Those are what the climate bandwidth excludes. Plots beyond three geographic bandwidths are never shown either way.")
                                    ),

                                    actionButton("highlight_site", "Show these plots on the map") %>%
                                          tooltip("Switch the map's vector layer to this site's analog plots and zoom to its geographic cutoff, so the same plots can be seen in geographic space."),
                                    br(), br(),

                                    layout_column_wrap(
                                          width = 1/2,
                                          plotOutput("clim_plot"),
                                          plotOutput("dist_plot")
                                    )
                              )

                        )

                  ),

                  ## MAP ---------------------------------------
                  leafletOutput("map", height = "100%")
            )
      ),

      add_busy_spinner(spin = "cube-grid", position = "bottom-right")
)



# SERVER ---------------------------------------------------
server <- function(input, output, session) {

      point_radius <- reactive({10^input$vector_scale * 1000})

      green <- "darkgreen"

      # Active focal identifier and human-readable label
      ## All-entity ("community composition") mode. Guarded on track as well as
      ## the control, because input$entity_scope keeps its value while hidden --
      ## switching to All trees with the scope left on "multi" would otherwise
      ## try to run a community reduction over a track that is already an
      ## aggregate.
      multi_mode <- reactive({
            isTruthy(input$entity_scope) && input$entity_scope == "multi" &&
                  isTruthy(input$track) && input$track != "All trees"
      })

      focal_id <- reactive({
            if (multi_mode())
                  return(paste0("all_", tolower(gsub(" ", "_", input$track))))
            switch(input$track,
                   "Species"      = input$species,
                   "Forest type"  = input$forest_type,
                   "Forest group" = input$forest_group,
                   "All trees"    = "all_trees")
      })

      focal_label <- reactive({
            if (multi_mode()) return(switch(input$track,
                                            "Species"      = "All species",
                                            "Forest type"  = "All forest types",
                                            "Forest group" = "All forest groups"))
            switch(input$track,
                   "Species"      = input$species,
                   "Forest type"  = {
                         nm <- ft_name(as.integer(input$forest_type))
                         if (is.na(nm)) paste("Forest type", input$forest_type) else nm
                   },
                   "Forest group" = {
                         nm <- fg_name(as.integer(input$forest_group))
                         if (is.na(nm)) paste("Forest group", input$forest_group) else nm
                   },
                   "All trees"    = "All trees")
      })

      # TRUE if the active response is a 0/1 indicator (any presence variable:
      # "presence probability", "seedling presence", and the forest type/group
      # classification probabilities); FALSE if continuous (basal area / counts).
      response_is_binary <- reactive({
            # A community reduction is never binary even though it is built from
            # presence columns: richness is a count.
            if (multi_mode()) return(FALSE)
            isTruthy(input$fia_var) && grepl("presence", input$fia_var)
      })

      # Sequential color domain c(0, max) for markers/points, robust to the case
      # where the values are all NA or all 0 (e.g. a fraction variable for a
      # species absent from the current view) -- an empty/degenerate range makes
      # colorNumeric error, so fall back to c(0, 1).
      seq_domain <- function(x) {
            mx <- suppressWarnings(max(x, na.rm = TRUE))
            c(0, if (is.finite(mx) && mx > 0) mx else 1)
      }

      # Lower/upper bounds for the currently-displayed raster. Used by the
      # density plot to keep KDE mass from leaking past natural limits (e.g.,
      # presence probability is bounded [0, 1]; total BA is non-negative).
      # Returns c(-Inf, Inf) when no informative bound applies.
      display_bounds <- reactive({
            if (input$stat == "ESS") {
                  c(0, Inf)
            } else if (input$stat == "focal site similarity") {
                  c(0, Inf)
            } else if (input$stat == "connectivity") {
                  if (input$cs_stat == "current direction") {
                        c(0, 360)
                  } else {
                        # other connectivity stats: non-negative magnitudes
                        c(0, Inf)
                  }
            } else if (input$stat == "FIA variable") {
                  if (input$time == "delta") {
                        c(-Inf, Inf)
                  } else if (response_is_binary()) {
                        c(0, 1)
                  } else if (grepl("fraction", input$fia_var)) {
                        c(0, 1)  # regeneration fractions are bounded [0, 1]
                  } else {
                        c(0, Inf)
                  }
            } else {
                  c(-Inf, Inf)
            }
      })

      # TRUE when the FIA-plot markers (map and focal-site scatterplots) should
      # adopt the raster's palette. The markers always show historical FIA
      # response on a [0, max] sequential scale; they're commensurable with the
      # raster only when the raster is also showing FIA response on a
      # sequential scale -- i.e. stat == "FIA variable" and time is historical
      # or future (delta flips the raster to a diverging zero-centered scale).
      # In all other cases the markers keep their own canonical gray90 -> green
      # ramp. See palette discussion: this preserves the raster/marker
      # alignment that holds by default without ever giving the markers an
      # inappropriate scale.
      marker_follows_raster_palette <- reactive({
            # Not for categorical layers: the map markers take the raster's
            # colorFactor directly, and the focal-site scatterplots are
            # continuous, so handing them a qualitative palette name would
            # reach scale_fill_distiller() with a palette it does not know.
            input$stat == "FIA variable" && input$time != "delta" && !categorical_layer()
      })

      ## The panel plots are always categorical, but rast_palette holds a
      ## continuous ramp whenever the map is showing a continuous measure -- so
      ## without this they would fall back to the default palette the moment the
      ## user switched the map off dominant entity. Remembering the last
      ## qualitative choice keeps the panel on the palette the user picked.
      last_cat_palette <- reactiveVal("default")
      observeEvent(input$rast_palette, {
            if (input$rast_palette %in% rast_palettes_categorical)
                  last_cat_palette(input$rast_palette)
      })

      ## Swap the palette choices to match the layer type. observeEvent so the
      ## handler's read of input$rast_palette is isolated and cannot re-trigger
      ## its own write. "default" is in both sets, so it is the natural landing
      ## point when a selection does not survive the swap.
      observeEvent(categorical_layer(), {
            ch <- if (categorical_layer()) rast_palettes_categorical
            else rast_palettes_continuous
            cur <- if (isTruthy(input$rast_palette)) input$rast_palette else "default"
            updateSelectInput(session, "rast_palette", choices = ch,
                              selected = if (cur %in% ch) cur else "default")
      }, ignoreInit = TRUE)

      # ggplot fill-scale equivalent of the leaflet choose_palette() logic.
      # Resolve one of rast_palettes_categorical to n colors. Falls back to the
      # default hue ramp for any name it does not recognise, which is what makes
      # the choice-set swap safe: a stale continuous name arriving here during a
      # flush produces a usable palette rather than an error.
      # Drop neutrals: near-black, near-white and low-saturation grays. Gray is
      # reserved for the pooled "other" band, so a gray in the palette reads as
      # "other"; near-white leaves an entity with no visible swatch at all; and
      # black is heavy enough to dominate a map it only holds a slice of.
      # Okabe-Ito leads with black and Polychrome 36's second entry is near-white,
      # so both need this.
      drop_neutral_colors <- function(cols) {
            h <- grDevices::rgb2hsv(grDevices::col2rgb(cols))
            # Saturation alone catches grays AND near-whites (near-white is
            # bright but unsaturated); the value floor catches near-black. An
            # upper value bound would be wrong -- a pure hue like magenta or
            # green sits at v = 1 and is vivid, not washed out.
            keep <- h["s", ] > 0.25 & h["v", ] > 0.35
            if (!any(keep)) return(cols)
            cols[keep]
      }

      # When more entities are shown than the palette holds, recycled colors are
      # shaded rather than repeated exactly. Two entities the same color is
      # unavoidable past the palette length, but a darker and a lighter version
      # of one hue can still be told apart, where two identical swatches cannot.
      recycle_shaded <- function(cols, n) {
            if (n <= length(cols)) return(cols[seq_len(n)])
            k <- length(cols)
            vapply(seq_len(n), function(i) {
                  base <- cols[((i - 1) %% k) + 1]
                  cycle <- (i - 1) %/% k
                  if (cycle == 0) return(base)
                  h <- grDevices::rgb2hsv(grDevices::col2rgb(base))
                  # Alternate darker / lighter on successive passes.
                  darker <- (cycle %% 2) == 1
                  v <- if (darker) max(0.20, h["v", 1] * 0.62) else min(0.98, h["v", 1] * 1.30)
                  sat <- if (darker) h["s", 1] else h["s", 1] * 0.55
                  grDevices::hsv(h["h", 1], min(1, sat), v)
            }, character(1))
      }

      cat_palette_colors <- function(n, pal) {
            if (!pal %in% rast_palettes_categorical) pal <- "default"
            if (n < 1) return(character(0))

            pick <- function(name, max_n) {
                  cols <- drop_neutral_colors(
                        unname(grDevices::palette.colors(max_n, palette = name)))
                  recycle_shaded(cols, n)
            }

            switch(pal,
                   "Polychrome" = pick("Polychrome 36", 36),
                   "Alphabet"   = pick("Alphabet", 26),
                   "Tableau"    = pick("Tableau 10", 10),
                   "Set 3"      = pick("Set 3", 12),
                   "Okabe-Ito"  = pick("Okabe-Ito", 8),
                   # Default adapts: Tableau 10 is the nicer set at small n, but
                   # it runs out, and past that perceptual separation matters
                   # more than individual colors looking good.
                   if (n <= 9) pick("Tableau 10", 10) else pick("Polychrome 36", 36))
      }

      # rast_palette values are heterogeneous: "default" -> the canonical
      # gray90 -> green gradient; "viridis"/"magma"/"plasma"/"inferno" ->
      # viridis options; the rest ("RdYlBu"/"Spectral"/"RdBu"/"BrBG") ->
      # RColorBrewer palettes via scale_fill_distiller. `limits` is passed
      # through so the scale matches the marker domain; `guide` defaults to
      # "none" since the focal-site plots don't show a fill legend.
      marker_fill_scale <- function(limits, guide = "none") {
            p <- input$rast_palette
            if (p == "default" || !marker_follows_raster_palette()) {
                  ggplot2::scale_fill_gradientn(colors = c("gray90", green),
                                                limits = limits, guide = guide)
            } else if (p %in% c("viridis", "magma", "plasma", "inferno")) {
                  ggplot2::scale_fill_viridis_c(option = p,
                                                limits = limits, guide = guide)
            } else {
                  # RColorBrewer diverging/sequential palettes
                  ggplot2::scale_fill_distiller(palette = p,
                                                limits = limits, guide = guide)
            }
      }


      # Relative weight each analog plot carries in a prediction.
      #
      # analog_search() returns distances but not weights (pairs mode gives
      # env_dist and geog_dist only), and the search in this panel is called with
      # cutoff-only kernels -- so nothing upstream knows about weighting. It has
      # to be recomputed here.
      #
      # *** THIS FORMULA MUST MATCH THE KERNELS IN aim_kernels():
      #     kernel(max = theta * 3, weight = "gaussian", theta = theta)
      # Assumed Gaussian with theta as the standard deviation, and the joint
      # weight the product across the two dimensions. That reading comes from the
      # app's own help text ("bandwidths specify the standard deviations") and
      # from the cutoff sitting at three of them. If analogs uses a different
      # constant, the point sizes here are subtly wrong while everything still
      # looks reasonable -- so this is the line to check against the package. ***
      analog_weight <- function(geog_dist, env_dist) {
            w <- exp(-0.5 * (geog_dist / input$theta_geog)^2) *
                  exp(-0.5 * (env_dist / input$theta_clim)^2)
            # The cutoff is hard, not a tail: kernel(max = theta * 3) drops a
            # plot entirely rather than down-weighting it. Without this a plot
            # just past the cutoff gets exp(-4.5), which is small but not zero,
            # and is drawn as though it still counts for something.
            w[geog_dist > input$theta_geog * 3 | env_dist > input$theta_clim * 3] <- 0
            w
      }

      # Value column and fill scale for the analog plots, which color each plot
      # by the response at that plot. For a categorical layer that response is an
      # entity code, and handing codes to a continuous ramp paints arbitrary
      # integers on a gradient -- the points come out gray-to-green and mean
      # nothing. Convert to a labeled factor and use the entity colors instead,
      # the same ones the map and the composition panel use.
      analog_fill <- function(d) {
            if (categorical_layer()) {
                  cc <- comm_categories()
                  d$value <- factor(cc$labels[match(d$value, cc$codes)],
                                    levels = cc$labels)
                  # guide "none": the map legend beside these plots already keys
                  # the same colors, and a second copy would crowd two small panels.
                  list(d = d,
                       scale = ggplot2::scale_fill_manual(
                             values = stats::setNames(cc$colors, cc$labels),
                             drop = FALSE, guide = "none", na.value = "gray85"))
            } else {
                  dom <- if (!response_is_binary()) seq_domain(d$value) else 0:1
                  list(d = d, scale = marker_fill_scale(limits = dom))
            }
      }


      ## DATA SETUP -----------------------

      # Padded bbox: includes plots just outside the model bounding box so they
      # can still serve as analogs for cells near the bbox edge.
      # Pad tracks theta_geog * 3 (the hard cutoff used by analog_impact),
      # converted from km to degrees with a conservative ratio. 1° lat ≈ 111 km;
      # 1° lon at 50° N ≈ 71 km. Using /100 overestimates the pad by 10–30% at
      # CONUS latitudes, which is harmless.
      pool_bbox <- reactive({
            pad_deg <- input$theta_geog * 3 / 100
            list(x = c(input$bbox_x[1] - pad_deg, input$bbox_x[2] + pad_deg),
                 y = c(input$bbox_y[1] - pad_deg, input$bbox_y[2] + pad_deg))
      })

      # all plots within padded bbox (the analog pool)
      fia_plots <- reactive({
            pb <- pool_bbox()
            fia_all %>%
                  select(x, y) %>%
                  distinct() %>%
                  filter(between(x, pb$x[1], pb$x[2]),
                         between(y, pb$y[1], pb$y[2])) %>%
                  mutate(index = 1:nrow(.))
      })

      fia_plots_species_all <- reactive({
            fia_all %>%
                  # Recode NA forest_type to the Non-forest sentinel (0) so
                  # downstream filters and joins treat non-forest plots
                  # uniformly with real forest types.
                  mutate(forest_type = ifelse(is.na(forest_type), 0L, forest_type)) %>%
                  left_join(select(fia_plots(), x, y, index),
                            by = join_by(x, y))
      })

      plot_occ <- reactive({
            pb <- pool_bbox()
            # Plot-level type->group lookup (one row per forest_type code)
            ft_to_fg <- distinct(forest_type_codes, forest_type, forest_group)

            fia_plots_species_all() %>%
                  filter(between(x, pb$x[1], pb$x[2]),
                         between(y, pb$y[1], pb$y[2])) %>%
                  mutate(focal = species == input$species & !is.na(species)) %>%
                  group_by(x, y) %>%
                  summarize(
                        samp_method_cd = mean(samp_method_cd),
                        pres = as.integer(any(focal)), # focal species presence/absence (any life stage)
                        pres_seed = as.integer(any(focal & n_seedlings > 0)), # focal species seedling presence
                        pres_sap = as.integer(any(focal & n_saplings > 0)), # focal species sapling presence
                        pres_tree = as.integer(any(focal & n_trees > 0)), # focal species adult (tree) presence
                        batot = sum(basal_area[focal]), # total focal species basal area
                        baprop = weighted.mean(focal, basal_area), #  proportion ba
                        baplot = sum(basal_area), # total plot ba, all species
                        nttot = sum(n_trees[focal]), # focal species tree count
                        ntplot = sum(n_trees), # total plot tree count, all species
                        satot = sum(n_saplings[focal]), # focal species sapling count
                        saplot = sum(n_saplings), # total plot sapling count, all species
                        nstot = sum(n_seedlings[focal]), # focal species seedling count
                        nsplot = sum(n_seedlings), # total plot seedling count, all species
                        # Plot species richness, by stage. A count of distinct
                        # species is an aggregation ACROSS species done at the
                        # plot, exactly like baplot and ntplot -- which is what
                        # makes it an All-trees variable rather than something
                        # only the multi-entity path can produce. Because AIM is
                        # a weighted mean (linear) with entity-independent
                        # weights, interpolating this equals summing the
                        # interpolated per-entity occurrence surfaces. The two
                        # routes should agree to floating point, which makes each
                        # a check on the other.
                        nspp      = dplyr::n_distinct(species[!is.na(species)]),
                        nspp_seed = dplyr::n_distinct(species[!is.na(species) & n_seedlings > 0]),
                        nspp_sap  = dplyr::n_distinct(species[!is.na(species) & n_saplings > 0]),
                        nspp_tree = dplyr::n_distinct(species[!is.na(species) & n_trees > 0]),
                        forest_type = dplyr::first(forest_type), # plot-level FORTYPCD
                        .groups = "drop") %>%
                  left_join(ft_to_fg, by = "forest_type") %>%
                  mutate(baprop = ifelse(baplot == 0, 0, baprop),
                         # Per-plot regeneration fractions (focal species), on an
                         # occurrence basis to match the interpolated surface (see
                         # ratio_pieces): 0 = adults only, 0.5 = stage and adults
                         # both present, 1 = stage only. NA where the species is
                         # absent as both that stage and adults.
                         seed_frac = ifelse((nstot + nttot) > 0,
                                            (nstot > 0) / ((nstot > 0) + (nttot > 0)), NA_real_),
                         sap_frac  = ifelse((satot + nttot) > 0,
                                            (satot > 0) / ((satot > 0) + (nttot > 0)), NA_real_))
      })

      # Single active response vector, aligned with plot_occ() rows.
      # Computed from the (track, fia_var, focal_id) tuple.
      response_values <- reactive({
            po <- plot_occ()
            req(isTruthy(input$track), isTruthy(input$fia_var))

            v <- if (input$track == "Species") {
                  switch(input$fia_var,
                         "presence probability"  = po$pres,
                         "seedling presence"      = po$pres_seed,
                         "sapling presence"       = po$pres_sap,
                         "tree presence"          = po$pres_tree,
                         "proportion basal area" = po$baprop,
                         "total basal area"      = po$batot,
                         "number of trees"       = po$nttot,
                         "number of saplings"    = po$satot,
                         "number of seedlings"   = po$nstot,
                         "seedling fraction"      = po$seed_frac,
                         "sapling fraction"       = po$sap_frac)
            } else if (input$track == "Forest type") {
                  req(isTruthy(input$forest_type))
                  as.integer(!is.na(po$forest_type) &
                                   po$forest_type == as.integer(input$forest_type))
            } else if (input$track == "Forest group") {
                  req(isTruthy(input$forest_group))
                  as.integer(!is.na(po$forest_group) &
                                   po$forest_group == as.integer(input$forest_group))
            } else if (input$track == "All trees") {
                  switch(input$fia_var,
                         "presence probability" = as.integer((po$ntplot + po$saplot) > 0 | po$nsplot > 0),
                         "seedling presence"    = as.integer(po$nsplot > 0),
                         "sapling presence"     = as.integer(po$saplot > 0),
                         "tree presence"        = as.integer(po$ntplot > 0),
                         "total basal area"     = po$baplot,
                         "number of trees"      = po$ntplot,
                         "number of saplings"   = po$saplot,
                         "number of seedlings"  = po$nsplot,
                         "species richness"          = po$nspp,
                         "seedling species richness" = po$nspp_seed,
                         "sapling species richness"  = po$nspp_sap,
                         "adult species richness"    = po$nspp_tree,
                         # Occurrence-based, matching ratio_pieces() and the
                         # Species-track fractions in plot_occ(): 0 = adults only,
                         # 0.5 = stage and adults both present, 1 = stage only.
                         "seedling fraction"    = ifelse((po$nsplot + po$ntplot) > 0,
                                                         (po$nsplot > 0) / ((po$nsplot > 0) + (po$ntplot > 0)), NA_real_),
                         "sapling fraction"     = ifelse((po$saplot + po$ntplot) > 0,
                                                         (po$saplot > 0) / ((po$saplot > 0) + (po$ntplot > 0)), NA_real_))
            }

            # Guard against transient NULL during track/fia_var transitions
            req(!is.null(v))
            v
      })

      # Whether the focal entity is present, per track -- the same quantity the
      # "presence probability" variable measures, and independent of which
      # variable is being displayed. Used by the "occupied FIA plots" filter:
      # keying that off `display > 0` means "has regeneration" for the fraction
      # variables, which hides plots that are occupied but simply aren't
      # regenerating. (It coincides with display > 0 for the other variables.)
      ## Per-plot value to paint the FIA markers with: the observed counterpart
      ## to whatever the raster is showing. In Single mode that is the focal
      ## entity's response; in All mode it is the same reduction applied across
      ## entities at each plot, so observed richness sits next to modelled
      ## richness on the same color scale. That comparison is the cheapest
      ## available check on whether Sum(p_s) is calibrated.
      display_values <- reactive({
            if (!multi_mode()) return(response_values())

            # Guard before touching multi_values(). The regeneration fractions
            # have no per-entity form, so multi_values() req()s out for them --
            # and because it would otherwise be passed to comm_reduce_rows() as
            # an unforced promise, the reactiveStop would surface from inside
            # as.matrix()'s S4 dispatch as an error rather than a silent cancel.
            # Forcing it here, in this reactive's own frame, keeps a cancellation
            # a cancellation.
            if (!isTruthy(input$fia_var) || !input$fia_var %in% multi_vars) return(NULL)
            y <- multi_values()
            req(!is.null(y), ncol(y) > 0)

            v <- comm_reduce_rows(y, input$comm_measure)
            # NULL means the measure has no per-plot counterpart (the paired
            # categoricals: a plot is surveyed at one point in time, so there is
            # no observed "biggest winner"). Return early -- subsetting NULL in
            # the lumping step below yields a ZERO-LENGTH vector rather than
            # NULL, which slips past a null check and only fails later at
            # `pd$display <- dv`.
            if (is.null(v)) return(NULL)

            # Lump the observed categories to exactly the set the raster drew, so
            # a plot and the cell under it get the same color for the same
            # entity, and anything outside the top N reads as "other" on both.
            if (input$comm_measure %in% categorical_measures) {
                  cc <- comm_categories()
                  v[!is.na(v) & !v %in% cc$codes] <- 0L
            }
            v
      })

      occupied_values <- reactive({
            po <- plot_occ()
            req(isTruthy(input$track))

            # In All mode "occupied" means any entity is present -- a forested
            # plot -- since there is no focal entity to be occupied by.
            #
            # Assigned before use for the same reason as display_values(): forcing
            # multi_values() inside an operator means a failed req() escapes from
            # dispatch as an error rather than cancelling.
            if (multi_mode()) {
                  if (!isTruthy(input$fia_var) || !input$fia_var %in% multi_vars)
                        return(rep(FALSE, nrow(po)))
                  y <- multi_values()
                  return(rowSums(y > 0, na.rm = TRUE) > 0)
            }

            v <- if (input$track == "Species") {
                  po$pres
            } else if (input$track == "Forest type") {
                  req(isTruthy(input$forest_type))
                  as.integer(!is.na(po$forest_type) &
                                   po$forest_type == as.integer(input$forest_type))
            } else if (input$track == "Forest group") {
                  req(isTruthy(input$forest_group))
                  as.integer(!is.na(po$forest_group) &
                                   po$forest_group == as.integer(input$forest_group))
            } else if (input$track == "All trees") {
                  as.integer((po$ntplot + po$saplot) > 0 | po$nsplot > 0)
            }

            req(!is.null(v))
            as.logical(v)
      })

      # Plot-level data behind the FIA marker layer, with the active response
      # attached as `display` and focal-entity occupancy as `occupied`. Shared by
      # the map observer and the CSV export so the downloaded table is exactly
      # the set of points drawn on the map.
      #
      # When the vector selector is on a non-point layer ("none", current
      # animation/field) there is nothing on screen to mirror, so this falls back
      # to the full analog pool -- the export stays available and meaningful
      # regardless of what the map happens to be showing.
      vector_points <- reactive({
            req(plot_occ())
            # display_values() is NULL if a community measure ever lacks a
            # row-wise counterpart. Catching it here degrades to an empty layer;
            # letting it through assigns NULL to pd$display, which silently drops
            # the column and crashes inside addCircles() instead.
            dv <- display_values()
            # Length, not just non-NULL: a zero-length vector would fail later
            # inside the data-frame assignment with a recycling error rather
            # than cancelling here.
            req(!is.null(dv), length(dv) == nrow(plot_occ()))

            if(input$vector == "focal site analogs"){
                  pd <- switch(input$analog_direction,
                               "Reverse analogs" = site_analogs()$rev_analogs,
                               "Forward analogs" = site_analogs()$fwd_analogs,
                               "Contemporary analogs" = site_analogs()$bsl_analogs) %>%
                        mutate(x = analog_x, y = analog_y)
                  # Index the per-plot value by analog_index rather than
                  # reading a `response` y-column: pool()$values carries
                  # occ/num/denom (not `response`) for the pooled-ratio
                  # variables, so display_values() is the uniform source.
                  pd$display <- dv[pd$analog_index]
                  pd$occupied <- occupied_values()[pd$analog_index]
            } else {
                  # Show the full analog pool, including plots in the padded
                  # ring just outside the modeling bbox — these are still
                  # contributing as analogs to edge cells in the prediction.
                  pd <- plot_occ()
                  pd$display <- dv
                  pd$occupied <- occupied_values()
            }

            # Occupancy is the focal entity's presence, not `display > 0`:
            # for the regeneration fractions the latter would mean "has
            # regeneration" and hide adults-only plots, which are occupied.
            if(input$vector == "occupied FIA plots") pd <- filter(pd, occupied)

            pd
      })

      # Header-only version of the export schema, used when export_points()
      # can't be built (see the download handler). Keeps the file a valid CSV
      # with the documented columns rather than a Shiny error page.
      export_points_empty <- function(){
            tibble::tibble(x = numeric(), y = numeric(),
                           entity = character(), variable = character(),
                           value = numeric(), occupied = logical(),
                           sampled = logical())
      }

      # Curated export schema for the FIA point download.
      #
      # Deliberately narrow, and deliberately FIXED: the header is identical no
      # matter which variable is displayed, because the variable's name lives in
      # a `variable` column rather than in the column name of its values. That
      # means several exports can be bind_rows()'d, and a script can read `value`
      # without knowing what was on screen when the file was made.
      #
      # Everything else in plot_occ() (batot, nsplot, samp_method_cd, ...) is
      # internal shorthand. Emitting it would turn those names into a public
      # interface that has to be documented here and kept stable across refactors,
      # which is a bigger commitment than the columns are worth.
      export_points <- reactive({
            pd <- vector_points()

            out <- tibble::tibble(
                  x = pd$x,
                  y = pd$y,
                  entity = focal_label(),
                  # The forest-type/group tracks ignore fia_var: their response is
                  # membership in the type named by `entity`.
                  variable = if(input$track %in% c("Forest type", "Forest group")) "presence"
                  else input$fia_var,
                  value = pd$display,
                  # Focal entity presence, independent of `variable`. Not the same
                  # question as value > 0: for the regeneration fractions that
                  # would mean "has regeneration", which is false at occupied
                  # adults-only plots.
                  occupied = pd$occupied,
                  # FALSE at plots that exist in the frame but were never surveyed
                  # (samp_method_cd 2). Without this their zeros read downstream as
                  # true absences.
                  sampled = pd$samp_method_cd == 1
            )

            # The analog layer additionally carries the kernel geometry that
            # defines the analog set. Weight matches the plot-level decomposition
            # code: the product of the two gaussian kernel densities, unnormalized.
            if(all(c("env_dist", "geog_dist") %in% names(pd))){
                  out <- out %>%
                        mutate(climate_distance = pd$env_dist,
                               geographic_distance = pd$geog_dist,
                               weight = dnorm(pd$env_dist, sd = input$theta_clim) *
                                     dnorm(pd$geog_dist, sd = input$theta_geog))
            }

            out
      })

      # Pooled-ratio variables (proportion basal area, seedling/sapling fraction)
      # are interpolated as a ratio-of-means rather than a mean-of-ratios: we
      # interpolate the numerator and denominator separately and divide the
      # surfaces (then gate by occurrence in aim()). This returns the per-plot
      # numerator/denominator for those variables, or NULL for everything else
      # (which uses the single-response path). Occurrence is denom > 0.
      ## MULTI-ENTITY RESPONSE MATRIX -------------------
      ##
      ## plot_occ() collapses to one row per plot around a single focal entity
      ## (`focal = species == input$species`), so there is nowhere for a
      ## multi-column `y` to come from. This builds the wide counterpart: one
      ## row per plot, one column per entity, holding the same quantity
      ## response_values() would return if that entity were the focal one.
      ##
      ## Row alignment with plot_occ() is the whole ballgame. pool()$index_hst is
      ## built from plot_occ()[, c("x","y")], so a mismatch here silently attaches
      ## every entity's response to the wrong plot and produces maps that look
      ## entirely plausible. Alignment is guaranteed by construction (the join
      ## runs FROM plot_occ()'s coordinates) and then asserted on the coordinates
      ## themselves rather than on row counts.
      multi_values <- reactive({
            req(isTruthy(input$track))
            po <- plot_occ()
            xy <- po[, c("x", "y")]

            if (input$track %in% c("Forest type", "Forest group")) {
                  # Membership is plot-level -- one type per plot -- so the wide
                  # form is an indicator matrix and fia_var does not apply. Same
                  # reason response_values() ignores fia_var for these tracks.
                  key <- if (input$track == "Forest type") po$forest_type else po$forest_group
                  lvls <- sort(unique(key[!is.na(key)]))
                  req(length(lvls) > 0)
                  w <- as.data.frame(
                        vapply(lvls, function(k) as.integer(!is.na(key) & key == k),
                               integer(nrow(po))),
                        check.names = FALSE)

                  # Name the columns for the classes rather than their FORTYPCD
                  # integers. The column names ARE the entity labels everywhere
                  # downstream -- legend, category table, and any future ranking
                  # -- so resolving them here fixes all of those at once rather
                  # than translating codes at each display site.
                  nm <- if (input$track == "Forest type") ft_name(as.integer(lvls))
                  else fg_name(as.integer(lvls))
                  nm[is.na(nm)] <- paste(input$track, lvls[is.na(nm)])
                  names(w) <- make.unique(as.character(nm))

            } else {
                  req(isTruthy(input$fia_var), input$fia_var %in% multi_vars)
                  pb <- pool_bbox()

                  # Same summaries plot_occ() computes, but grouped by entity as
                  # well as plot. All of them at once: it is a single pass over
                  # the same rows either way, and computing them together keeps
                  # the definitions visibly identical to the single-entity path.
                  long <- fia_plots_species_all() %>%
                        filter(between(x, pb$x[1], pb$x[2]),
                               between(y, pb$y[1], pb$y[2]),
                               !is.na(species)) %>%
                        group_by(x, y, species) %>%
                        summarize(
                              # Presence is "the entity was recorded at this
                              # plot", matching plot_occ()'s any(focal) -- it is
                              # not conditioned on any stage count.
                              pres      = 1L,
                              pres_seed = as.integer(any(n_seedlings > 0)),
                              pres_sap  = as.integer(any(n_saplings > 0)),
                              pres_tree = as.integer(any(n_trees > 0)),
                              batot     = sum(basal_area),
                              nttot     = sum(n_trees),
                              satot     = sum(n_saplings),
                              nstot     = sum(n_seedlings),
                              .groups = "drop")

                  val_col <- switch(input$fia_var,
                                    "presence probability"  = "pres",
                                    "seedling presence"     = "pres_seed",
                                    "sapling presence"      = "pres_sap",
                                    "tree presence"         = "pres_tree",
                                    "total basal area"      = "batot",
                                    "proportion basal area" = "batot",
                                    "number of trees"       = "nttot",
                                    "number of saplings"    = "satot",
                                    "number of seedlings"   = "nstot")

                  wide <- long %>%
                        select(x, y, species, value = dplyr::all_of(val_col)) %>%
                        tidyr::pivot_wider(names_from = species, values_from = value,
                                           values_fill = 0)

                  # Join FROM plot_occ()'s coordinates, so row order is inherited
                  # rather than assumed. Plots where no entity was recorded are
                  # dropped by the !is.na(species) filter above and must come
                  # back as all-zero rows: absence is data, and dropping those
                  # plots would bias every entity's surface upward.
                  j <- xy %>% left_join(wide, by = dplyr::join_by(x, y))
                  stopifnot(nrow(j) == nrow(xy),
                            identical(j$x, xy$x), identical(j$y, xy$y))

                  w <- as.data.frame(j[, setdiff(names(j), c("x", "y")), drop = FALSE],
                                     check.names = FALSE)
                  w[is.na(w)] <- 0

                  if (identical(input$fia_var, "proportion basal area")) {
                        # Denominator is the plot total across ALL entities, so
                        # it is shared. Plots with no basal area get 0 rather
                        # than NaN, matching plot_occ()'s baprop guard.
                        denom <- po$baplot
                        w[] <- lapply(w, function(v) ifelse(denom == 0, 0, v / denom))
                  }
            }

            # Occupied-plot count per entity: the pool-level support signal, and
            # the basis for deciding a rare-entity floor. Attached rather than
            # filtered on, so that decision can be made against real numbers
            # instead of a guessed threshold.
            attr(w, "n_plots") <- vapply(w, function(v) sum(v > 0, na.rm = TRUE), integer(1))
            w
      })

      ratio_pieces <- reactive({
            req(isTruthy(input$track), isTruthy(input$fia_var))
            po <- plot_occ()
            # Regeneration fractions are a fraction of OCCURRENCE weight, not of
            # stem density. Seedlings and saplings are tallied on the microplot
            # (each stem = 74.965/acre) while adults come off the subplot (each =
            # 6.018/acre), so a density ratio applies a ~12.5x amplification to
            # every regen stem before dividing: one seedling + one adult reads
            # 0.93 rather than 0.5, squashing the balanced core up against the
            # pure-colonization end. Using presence indicators puts both stages on
            # the same footing and gives a readable axis: 0 = adults only
            # (relict/trailing), ~0.5 = both present (core), 1 = regen only
            # (colonizing/leading).
            #
            # Proportion basal area is deliberately left on a density basis: it is
            # focal BA over stand BA, with each stem's own TPA applied to both
            # sides, so it never crosses the microplot/subplot boundary and isn't
            # distorted the way the regen fractions are.
            if (input$track == "Species") {
                  switch(input$fia_var,
                         "proportion basal area" = list(num = po$batot, denom = po$baplot),
                         "seedling fraction"     = list(num = as.integer(po$nstot > 0),
                                                        denom = as.integer(po$nstot > 0) +
                                                              as.integer(po$nttot > 0)),
                         "sapling fraction"      = list(num = as.integer(po$satot > 0),
                                                        denom = as.integer(po$satot > 0) +
                                                              as.integer(po$nttot > 0)),
                         NULL)
            } else if (input$track == "All trees") {
                  switch(input$fia_var,
                         "seedling fraction" = list(num = as.integer(po$nsplot > 0),
                                                    denom = as.integer(po$nsplot > 0) +
                                                          as.integer(po$ntplot > 0)),
                         "sapling fraction"  = list(num = as.integer(po$saplot > 0),
                                                    denom = as.integer(po$saplot > 0) +
                                                          as.integer(po$ntplot > 0)),
                         NULL)
            } else NULL
      })

      clim_scales <- readRDS("data/clim_scales.rds")

      # Tracks whether the user has clicked Clear after the most recent upload.
      # The Clear button sets this to TRUE; a fresh upload resets it to FALSE.
      # polygon_upload_active() returns input$polygon_upload only when not
      # cleared, so downstream polygon reactives see "no polygon" after Clear
      # without needing to actually reset the fileInput widget.
      polygon_cleared <- reactiveVal(FALSE)

      observeEvent(input$polygon_upload, {
            polygon_cleared(FALSE)
      })

      observeEvent(input$polygon_clear, {
            polygon_cleared(TRUE)
      })

      polygon_upload_active <- reactive({
            if (isTRUE(polygon_cleared())) NULL else input$polygon_upload
      })

      user_polygon <- reactive({
            req(polygon_upload_active())
            files <- polygon_upload_active()
            exts <- tolower(tools::file_ext(files$name))

            tryCatch({
                  if("shp" %in% exts){
                        tmpdir <- tempfile()
                        dir.create(tmpdir)
                        file.copy(files$datapath, file.path(tmpdir, files$name))
                        shp_path <- file.path(tmpdir, files$name[exts == "shp"])
                        st_read(shp_path, quiet = TRUE)
                  } else if(any(exts %in% c("kml", "kmz", "gpkg", "geojson"))){
                        st_read(files$datapath[1], quiet = TRUE)
                  } else {
                        showNotification("Unsupported file format.", type = "error")
                        NULL
                  }
            }, error = function(e){
                  showNotification(paste("Error reading polygon:", conditionMessage(e)), type = "error")
                  NULL
            }) %>%
                  st_transform(4326)
      })


      ## focal site
      site <- reactiveValues(x = NULL, y = NULL)



      # REACTIVE UI COMPONENTS ---------------------------

      ## Variable selection is a (track x metric x stage) cascade that resolves to
      ## a single canonical `fia_var` string, which the rest of the app switches
      ## on. These helpers define the valid choices and the resolution.
      ## Repopulate variable choices when the track or life stage changes.
      ## Selection is preserved when it survives the new choice set, so moving
      ## between two stages that both offer "occurrence" does not bounce the
      ## user back to the first option.
      observeEvent(list(input$track, input$stage), {
            mc <- metric_choices(input$track, input$stage)
            cur <- if (isTruthy(input$metric)) input$metric else ""
            updateSelectInput(session, "metric", choices = mc,
                              selected = if (cur %in% mc) cur else mc[1])
      }, ignoreInit = TRUE)

      ## Entity-scope option text follows the track ("single species" vs "single
      ## forest type"); only the labels change, the values stay single/multi so
      ## nothing downstream is affected.
      ##
      ## observeEvent rather than observe: the handler body is isolated, so
      ## reading input$entity_scope here does not create a dependency on the
      ## input this writes to. A bare observe() would re-fire on its own write.
      observeEvent(input$track, {
            req(input$track != "All trees")
            one <- switch(input$track,
                          "Species"      = "single species",
                          "Forest type"  = "single forest type",
                          "Forest group" = "single forest group")
            ch <- setNames(c("single", "multi"), c(one, "community composition"))
            updateSelectInput(session, "entity_scope", choices = ch,
                              selected = if (isTruthy(input$entity_scope)) input$entity_scope else "single")
      }, ignoreInit = TRUE)

      ## Resolve (track, metric, stage) -> hidden canonical fia_var.
      observeEvent(list(input$track, input$metric, input$stage), {
            req(input$metric)
            v <- resolve_var(input$track, input$metric, input$stage)
            req(!is.null(v))
            updateSelectInput(session, "fia_var", selected = v)
      })

      ## Seed the climate-model chart's x/y axis variables from the user's
      ## left-sidebar climate-variable selection. Done once, on the first
      ## populated value of input$clim_vars, so it provides a sensible default
      ## (the first two analog variables) without clobbering a deliberate axis
      ## choice the user makes later. The dropdowns still offer all six
      ## variables -- only the initial default tracks the sidebar.
      .model_delta_axes_seeded <- reactiveVal(FALSE)
      observeEvent(input$clim_vars, {
            if (isTRUE(.model_delta_axes_seeded())) return()
            cv <- input$clim_vars
            req(length(cv) >= 1)
            vx <- cv[1]
            vy <- if (length(cv) >= 2) cv[2] else cv[1]
            updateSelectInput(session, "model_delta_x", selected = vx)
            updateSelectInput(session, "model_delta_y", selected = vy)
            .model_delta_axes_seeded(TRUE)
      })

      ## Keep the climate-relationships module's variable menus in sync with the
      ## climate variables selected in the left sidebar. Its plot indexes into
      ## clim()$hst, which only holds the selected variables, so the options must
      ## be exactly that set. Preserve the current pick when it is still
      ## available, otherwise fall back to the first selected variable.
      ## cr_v1 offers every display variable, not just the modeled ones, so its
      ## choices no longer depend on input$clim_vars and the observer that used
      ## to narrow them is gone. The UI declares the full set directly.

      ## Spatial-scope status lines. These exist so the left panel states its
      ## own preconditions: a right-sidebar panel that renders empty because no
      ## site or polygon is set now has a visible cause rather than looking
      ## broken. Both are uiOutput rather than textOutput so the "not set" case
      ## can carry its own instruction.
      output$site_status <- renderUI({
            if (is.null(site$x) || is.null(site$y)) {
                  helpText("None selected \u2014 click the map to set one.")
            } else {
                  helpText(sprintf("lat %.3f, lon %.3f", site$y, site$x))
            }
      })

      output$polygon_status <- renderUI({
            f <- polygon_upload_active()
            if (!isTruthy(f)) {
                  helpText("None uploaded.")
            } else {
                  # A shapefile arrives as several sibling files; name the layer
                  # rather than listing .dbf / .shx / .prj back at the user.
                  nm <- f$name[tolower(tools::file_ext(f$name)) %in%
                                     c("shp", "kml", "kmz", "gpkg", "geojson")]
                  helpText(if (length(nm)) nm[1] else paste(length(f$name), "files"))
            }
      })

      ## Dynamic snap-to-range button (label reflects focal entity)
      output$snap_focal_btn <- renderUI({
            switch(input$track,
                   "Species"      = actionButton("focal_domain", "Snap to species range"),
                   "Forest type"  = actionButton("focal_domain", "Snap to forest type range"),
                   "Forest group" = actionButton("focal_domain", "Snap to forest group range"),
                   "All trees"    = tags$button(class = "btn btn-default",
                                                disabled = NA,
                                                style = "opacity: 0.5; cursor: not-allowed;",
                                                title = "No focal entity in this mode",
                                                "Snap to focal range"))
      })

      ## Snap domain to focal species/forest-type/forest-group range
      observeEvent(input$focal_domain, {
            pad <- 1
            bb <- switch(input$track,
                         "Species" = fia_plots_species_all() %>%
                               filter(species == input$species),
                         "Forest type" = fia_plots_species_all() %>%
                               filter(forest_type == as.integer(input$forest_type)),
                         "Forest group" = {
                               ft_in_group <- forest_type_codes %>%
                                     filter(forest_group == as.integer(input$forest_group)) %>%
                                     pull(forest_type)
                               fia_plots_species_all() %>%
                                     filter(forest_type %in% ft_in_group)
                         },
                         NULL)
            req(!is.null(bb), nrow(bb) > 0)
            bb <- bb %>%
                  summarize(xmin = floor(min(x)) - pad,
                            xmax = ceiling(max(x)) + pad,
                            ymin = floor(min(y)) - pad,
                            ymax = ceiling(max(y)) + pad)
            # bbox_x / bbox_y are sliderInputs, so updateSliderInput(value=) is
            # the matching call. These previously used updateSelectInput(selected=)
            # and worked, but that pairing is not the documented one and is not
            # something to rely on.
            updateSliderInput(session, "bbox_x", value = c(bb$xmin, bb$xmax))
            updateSliderInput(session, "bbox_y", value = c(bb$ymin, bb$ymax))
      })

      ## Snap domain to current map view
      observeEvent(input$window_domain, {
            bb <- input$map_bounds
            updateSliderInput(session, "bbox_x", value = c(bb$west, bb$east))
            updateSliderInput(session, "bbox_y", value = c(bb$south, bb$north))
      })

      # for circuitscape rasters, time period must be "delta"
      observe({
            if(input$stat == "connectivity"){
                  updateSelectInput(session, "time", selected = "delta")
            }
      })

      ## Show or hide the connectivity entries in the left-sidebar "Map layers"
      ## menus, driven by the "Enable connectivity" switch in the right-sidebar
      ## Connectivity panel. Connectivity is a niche analysis, so its options are
      ## kept out of the core menus until the user asks for them.
      ##
      ## updateSelectInput() resets the selection to the first choice whenever
      ## `choices` changes, so `selected` is always passed explicitly. On enable,
      ## the current pick is preserved; on disable, a connectivity pick is no
      ## longer offered and falls back to the menu default.
      ##
      ## ignoreInit: the UI is built with the connectivity-off choices and the
      ## switch defaults to FALSE, so the two already agree at startup.
      ## Explains the withdrawal rather than leaving a switch that has silently
      ## stopped doing anything.
      output$cs_note <- renderUI({
            req(multi_mode())
            div(class = "ts-notice",
                "Connectivity models the movement of one species or forest type between its
                 current and future range, so it does not apply to community composition. For
                 connectivity of the forest as a whole, set the entity to 'all trees'.")
      })

      ## Also fires on multi_mode(): connectivity is a single-entity analysis.
      ## Its electrodes are built from one entity's historical and future
      ## surfaces and then solved, so under community composition it would
      ## silently describe whichever species the hidden picker still holds --
      ## producing something that looks like community connectivity and is not.
      ##
      ## There is no community version worth building here. The one that would
      ## be well defined is connectivity of an aggregate surface, and the All
      ## trees track already provides exactly that.
      observeEvent(list(input$cs_enable, multi_mode()), {
            on <- isTRUE(input$cs_enable) && !multi_mode()

            stat_choices <- if(on) stat_choices_cs else stat_choices_base
            vector_choices <- if(on) vector_choices_cs else vector_choices_base

            keep <- function(cur, choices, default){
                  if(isTruthy(cur) && cur %in% choices) cur else default
            }

            updateSelectInput(session, "stat", choices = stat_choices,
                              selected = keep(input$stat, stat_choices, "FIA variable"))
            updateSelectInput(session, "vector", choices = vector_choices,
                              selected = keep(input$vector, vector_choices, "occupied FIA plots"))
      }, ignoreInit = TRUE)

      # Swap eval_stat choices to the metric set cv_performance() produces for
      # the current response type. Binary (presence-type) responses get
      # ranking/threshold metrics; continuous responses get error metrics.
      #
      # Keyed on response_is_binary() rather than on fia_var directly: in
      # community mode fia_var is still "presence probability" -- that is the
      # per-entity variable being summed -- but the response being scored is
      # richness, a count. Testing fia_var offered AUC and TSS for it.
      #
      # observeEvent so the handler is isolated; ignoreInit because the UI is
      # declared with the binary set, which is right for the default variable.
      observeEvent(response_is_binary(), {
            if (response_is_binary()) {
                  updateSelectInput(session, "eval_stat",
                                    choices = c("AUC", "TSS", "Brier"),
                                    selected = "AUC")
            } else {
                  updateSelectInput(session, "eval_stat",
                                    choices = c("R-squared", "RMSE"),
                                    selected = "R-squared")
            }
      }, ignoreInit = TRUE)

      # highlight focal site analogs: activate correct vector layer; zoom to site
      observeEvent(input$highlight_site, {
            updateSelectInput(session, "vector", selected = "focal site analogs")

            if(!is.null(site$y)){
                  circle <- geo_circle(cbind(site$x, site$y), input$theta_geog * 3) %>%
                        as.data.frame()

                  leafletProxy("map") %>%
                        fitBounds(lng1 = min(circle$lon),
                                  lat1 = min(circle$lat),
                                  lng2 = max(circle$lon),
                                  lat2 = max(circle$lat))
            }
      })

      # add user-provided polygon (and its buffered version)
      observe({
            proxy <- leafletProxy("map") %>%
                  clearGroup("user_polygon")

            # Re-draw after any raster change by taking dependency on display_layer
            display_layer()

            # Visibility toggle (Cartography panel). When off, the group is
            # already cleared above; req() stops before redrawing. Toggling
            # back on re-fires this observer and redraws. Mirrors show_site.
            req(input$show_polygon)

            req(user_polygon())
            poly <- user_polygon()

            proxy %>%
                  addPolygons(
                        data = user_polygon_buffered(),
                        color = lighten_color(polygon_color()),
                        weight = 2,
                        dashArray = "5, 5",
                        fill = FALSE,
                        group = "user_polygon"
                  ) %>%
                  addPolygons(
                        data = poly,
                        color = polygon_color(),
                        weight = 2,
                        fillOpacity = 0.1,
                        fillColor = polygon_color(),
                        group = "user_polygon",
                        label = if(!is.null(poly$Name)) ~Name else NULL
                  )
      })


      ## CLIMATE DATA -----------------------------------------

      # Raw (uncropped, unaggregated) historical and future climate rasters,
      # restricted to the selected climate variables. The future raster depends
      # on the scenario selection; the historical raster is shared across
      # scenarios. Used as the source for both the modeling grid (via clim())
      # and for plot-level climate extraction (in pool(), site_analogs(), and
      # eval_var_set()).

      # Active scenario filename(s), reconstructed from the three scenario
      # inputs. scenario_model may hold multiple GCMs; this returns one path
      # per selected model (the future raster is later averaged across them
      # to form a model ensemble). Strict: if ANY selected combination is
      # missing a file on disk, surfaces a user-visible message via validate()
      # naming the missing models, rather than silently averaging a subset.
      # Downstream outputs (the map, the focal-site plots, etc.) display the
      # message instead of erroring.
      scenario_file <- reactive({
            req(isTruthy(input$scenario_ssp),
                isTruthy(input$scenario_period),
                isTruthy(input$scenario_model))
            f <- paste0(input$scenario_ssp, "_",
                        input$scenario_period, "_",
                        input$scenario_model, ".tif")
            paths <- file.path("data/clim_future", f)
            exists <- file.exists(paths)

            shiny::validate(need(
                  all(exists),
                  if (length(paths) == 1) {
                        "Selected scenario combination is not available. Try a different emissions pathway, era, or climate model."
                  } else {
                        paste0("No data for these model(s) at the selected pathway/era: ",
                               paste(toupper(input$scenario_model[!exists]), collapse = ", "),
                               ". Deselect them or choose a different pathway/era.")
                  }))
            paths
      })

      clim_pool <- reactive({
            # Resolve scenario_file() at the top level so a req()/validate()
            # short-circuit (e.g., no selected (ssp, period, model) combination
            # exists on disk) propagates as a clean reactive cancellation
            # rather than an rast() argument-evaluation error.
            sf <- scenario_file()
            r_hst <- rast("data/clim_historical.tif")[[input$clim_vars]]

            # scenario_file() may return multiple paths (one per selected GCM).
            # Read each, subset to the active climate variables, and average
            # layer-wise to form a model ensemble. Subsetting *before* the sum
            # guarantees layers align by position == name across models, since
            # raster addition is positional. Reduce() (vs. do.call("+", ...))
            # is used so the single-model case is a clean no-op rather than an
            # error from calling `+` with one argument.
            fut_list <- lapply(sf, function(p) rast(p)[[input$clim_vars]])
            r_fut <- Reduce(`+`, fut_list) / length(fut_list)
            names(r_fut) <- input$clim_vars   # Reduce()/`+` can drop layer names

            list(hst = r_hst, fut = r_fut)
      })

      clim <- reactive({
            cp <- clim_pool()
            bb_ext <- ext(c(input$bbox_x, input$bbox_y))
            clim_hst <- cp$hst %>%
                  crop(bb_ext) %>%
                  aggregate(input$clim_agg, na.rm = TRUE)
            clim_fut <- cp$fut %>%
                  crop(bb_ext) %>%
                  aggregate(input$clim_agg, na.rm = TRUE)
            list(hst = clim_hst, fut = clim_fut)
      })


      ## DISPLAY-ONLY COVARIATES ------------------------------------
      ##
      ## The x-axis variable for the climate-relationships plot, resolved
      ## independently of clim_pool(). Nothing here feeds the analog model: the
      ## point is that a user can plot against a variable the model was not
      ## fitted on, and the isolation is what makes that safe.
      ##
      ## The cost is one extra crop/aggregate per variable change. Cheap relative
      ## to a prediction, and it only recomputes when the axis variable, the
      ## bbox, the grain, or the scenario changes.
      cr_grid <- reactive({
            v <- input$cr_v1
            req(isTruthy(v), v %in% display_vars())
            bb_ext <- ext(c(input$bbox_x, input$bbox_y))

            grid_of <- function(r) r %>% crop(bb_ext) %>%
                  aggregate(input$clim_agg, na.rm = TRUE)

            if (!is.null(display_extra_vars[[v]])) {
                  # A static covariate: the same raster for both timeframes,
                  # since elevation does not change between them.
                  g <- grid_of(rast(display_extra_vars[[v]]$file))
                  names(g) <- v
                  return(list(hst = g, fut = g))
            }

            # Averaged across the selected GCMs, matching clim_pool()'s ensemble.
            sf <- scenario_file()
            fut_list <- lapply(sf, function(pth) rast(pth)[[v]])
            r_fut <- Reduce(`+`, fut_list) / length(fut_list)
            names(r_fut) <- v

            list(hst = grid_of(rast("data/clim_historical.tif")[[v]]),
                 fut = grid_of(r_fut))
      })

      ## Same variable at the FIA plot coordinates, row-aligned with plot_occ().
      cr_plot_values <- reactive({
            v <- input$cr_v1
            req(isTruthy(v), v %in% display_vars())
            xy <- as.data.frame(plot_occ()[, c("x", "y")])

            src <- if (!is.null(display_extra_vars[[v]])) {
                  rast(display_extra_vars[[v]]$file)
            } else if (isTRUE(input$cr_clim == "future")) {
                  sf <- scenario_file()
                  fl <- lapply(sf, function(pth) rast(pth)[[v]])
                  Reduce(`+`, fl) / length(fl)
            } else {
                  rast("data/clim_historical.tif")[[v]]
            }
            as.numeric(terra::extract(src, xy, ID = FALSE, method = "bilinear")[, 1])
      })

      ## Outline color for the focal polygon. selectizeInput with a swatch
      ## renderer rather than a text box: the swatches make the choice visual
      ## without adding a color-picker dependency, and `create = TRUE` still
      ## accepts any R color name or hex typed straight in.
      output$polygon_color_ui <- renderUI({
            req(isTruthy(polygon_upload_active()))
            # isolate: this renderUI CREATES polygon_color, so reading it
            # reactively would rebuild the widget on every change -- resetting
            # the dropdown and stealing focus mid-selection. The read is only to
            # carry the selection across a polygon swap.
            overlay_color_input("polygon_color", "Focal polygon color",
                                isolate(if (isTruthy(input$polygon_color)) input$polygon_color
                                        else "darkred"))
      })

      ## Matching control for the focal site, shown only once one is set.
      output$site_color_ui <- renderUI({
            req(!is.null(site$x), !is.null(site$y))
            overlay_color_input("site_color", "Focal site color",
                                isolate(if (isTruthy(input$site_color)) input$site_color
                                        else "#4c32a8"))
      })

      ## Resolved once so both the outline and its buffer ring agree, and so a
      ## typed value that R cannot parse falls back instead of erroring inside
      ## leaflet.
      valid_color <- function(v, fallback) {
            if (!isTruthy(v)) return(fallback)
            ok <- tryCatch({ grDevices::col2rgb(v); TRUE }, error = function(e) FALSE)
            if (ok) v else fallback
      }
      polygon_color <- reactive(valid_color(input$polygon_color, "darkred"))
      site_color    <- reactive(valid_color(input$site_color, "#4c32a8"))

      ## Restrict a raster to the selected area. Applied to the COMBINED raster
      ## rather than to each layer, because cropping them separately would leave
      ## the response and the covariate on different extents.
      cr_scope_mask <- function(r) {
            sc <- if (isTruthy(input$cr_scope)) input$cr_scope else "domain"
            if (identical(sc, "domain")) return(r)

            if (sc %in% c("polygon", "buffered")) {
                  req(isTruthy(polygon_upload_active()))
                  pv <- terra::vect(if (identical(sc, "polygon")) user_polygon()
                                    else user_polygon_buffered())
                  return(terra::mask(terra::crop(r, pv), pv))
            }

            # Focal site: the neighborhood the model draws analogs from, not the
            # single cell -- one cell would give one point and no relationship.
            req(!is.null(site$x), !is.null(site$y))
            b <- terra::buffer(terra::vect(cbind(site$x, site$y), crs = "EPSG:4326"),
                               width = input$theta_geog * 3 * 1000)
            terra::mask(terra::crop(r, b), b)
      }

      ## Is the x-axis variable one the model was actually fitted on? A
      ## relationship against anything else is a property of the landscape rather
      ## than of the model, and the plot says so rather than leaving it to be
      ## inferred.
      cr_var_in_model <- reactive({
            isTruthy(input$cr_v1) && input$cr_v1 %in% input$clim_vars
      })

      ## CLIMATE MODEL SELECTOR -------------------------------------

      # Per-(scenario x era x model) climate deltas (future - historical) for
      # the model-chooser module. Independent of input$clim_vars, the scenario
      # selection, and the AIM pipeline: reads ALL climate variables straight
      # from the raw rasters for *every* future file on disk, so the chart
      # shows the full ensemble as time trajectories rather than a single
      # selection snapshot. Extracting every variable (rather than just the
      # current axis pair) lets the user re-pick chart axes with no recompute.
      #
      # Deltas are back-converted from standardized (z-score) units -- how the
      # rasters are stored -- to native climate units by multiplying by each
      # variable's SD; see clim_units / clim_scales.
      #
      # The region over which each delta is summarized is user-controlled
      # (input$model_delta_region): the modeling-domain bbox, the focal site,
      # or the uploaded study-area polygon.
      #
      # Cost: one raw future-raster read per future file -- the full ensemble,
      # so this is heavy. It is therefore gated behind the "Update chart"
      # button via bindEvent(): nothing recomputes until the user asks, even
      # if the summary region changes in the meantime. Axis-variable changes
      # never recompute -- all variables are already extracted, so the plot
      # just reselects columns. bindCache() underneath still earns its place:
      # re-clicking Update after returning to a previously-seen region is a
      # cache hit rather than a re-read. (reactive + bindCache + bindEvent is
      # the composable equivalent of an eventReactive that also caches.)

      # Label used for the synthetic zero-delta "historical" anchor that roots
      # every model trajectory at the origin.
      .hist_era_label <- "1991-2010 (historical)"

      model_deltas <- reactive({
            region <- input$model_delta_region
            req(isTruthy(region))

            shiny::validate(need(nrow(scenario_meta) > 0,
                                 "No future climate files found in data/clim_future/."))

            if (region %in% c("cell", "site")) {
                  shiny::validate(need(!is.null(site$x) && !is.null(site$y),
                                       "Click the map to select a focal site, or choose a different summary region."))
            } else if (region %in% c("polygon", "buffered")) {
                  shiny::validate(need(isTruthy(polygon_upload_active()),
                                       "Upload a focal polygon, or choose a different summary region."))
            }

            # Extract ALL climate variables (not just the current axis pair),
            # so the user can re-pick chart axes without triggering a recompute
            # -- the plot just reselects columns from the cached result. One
            # extra cost is reading 6 raster layers per file instead of 2,
            # which is minor next to the per-file I/O.
            all_vars <- clim_vars()
            hst <- rast("data/clim_historical.tif")[[all_vars]]

            # Region geometry resolved once, outside the per-file loop. The
            # buffered polygon and the site neighborhood differ from their
            # unbuffered counterparts only slightly for a climate delta, but they
            # are offered anyway: a wide buffer around a small polygon can span a
            # real elevation gradient, and a user should not have to learn which
            # module offers which areas.
            bb_ext <- ext(c(input$bbox_x, input$bbox_y))
            region_vect <- switch(
                  region,
                  "polygon"  = terra::vect(st_transform(user_polygon(), crs(hst))),
                  "buffered" = terra::vect(st_transform(user_polygon_buffered(), crs(hst))),
                  "site"     = terra::buffer(
                        terra::vect(cbind(site$x, site$y), crs = "EPSG:4326"),
                        width = input$theta_geog * 3 * 1000),
                  NULL)

            # Summarize a (future - historical) delta SpatRaster to a named
            # numeric vector, one value per climate variable.
            summarize_delta <- function(delta) {
                  if (region == "cell") {
                        v <- terra::extract(delta, cbind(site$x, site$y),
                                            method = "bilinear")
                        stats::setNames(as.numeric(v[1, ]), names(delta))
                  } else if (!is.null(region_vect)) {
                        v <- terra::extract(delta, region_vect, fun = "mean",
                                            na.rm = TRUE, ID = FALSE)
                        stats::setNames(as.numeric(v[1, ]), names(delta))
                  } else {
                        terra::global(terra::crop(delta, bb_ext), "mean",
                                      na.rm = TRUE)[, 1] %>%
                              stats::setNames(names(delta))
                  }
            }

            # Per-variable SD for back-converting standardized deltas to native
            # units. mean cancels in a difference, so only SD is needed. A
            # missing scale entry falls back to 1 (delta stays standardized).
            var_sd <- vapply(all_vars, function(v) {
                  s <- clim_scales[[v]]$sd
                  if (is.null(s) || !is.finite(s)) 1 else s
            }, numeric(1))

            rows <- lapply(seq_len(nrow(scenario_meta)), function(i) {
                  m <- scenario_meta[i, ]
                  fut <- rast(file.path("data/clim_future", m$file))[[all_vars]]
                  d <- summarize_delta(fut - hst)
                  # Back-convert each delta from z-score units to native units.
                  d <- d[all_vars] * var_sd
                  tibble::tibble(scenario = m$ssp, scenario_label = m$ssp_label,
                                 era = m$period, model = m$model,
                                 model_label = m$model_label) %>%
                        bind_cols(tibble::as_tibble(as.list(d)))
            })
            futures <- dplyr::bind_rows(rows)

            # Synthetic historical anchor: one zero-delta row per
            # (scenario x model) trajectory, so every path starts at the
            # origin. Era label is the constant defined above.
            zero_delta <- stats::setNames(as.list(rep(0, length(all_vars))),
                                          all_vars)
            anchors <- futures %>%
                  distinct(scenario, scenario_label, model, model_label) %>%
                  mutate(era = .hist_era_label) %>%
                  bind_cols(tibble::as_tibble(zero_delta)[rep(1, nrow(.)), ])

            bind_rows(anchors, futures) %>%
                  # Era ordering: synthetic historical first, then future
                  # periods chronologically (string sort works for "YYYY-YYYY").
                  mutate(era = factor(era,
                                      levels = c(.hist_era_label,
                                                 sort(unique(futures$era)))))
      }) %>%
            # theta_geog is part of the key because the buffered polygon and the
            # site neighborhood are both sized by it -- without it, changing the
            # bandwidth would return a cached delta for the old geometry.
            bindCache(input$model_delta_region,
                      input$bbox_x, input$bbox_y,
                      site$x, site$y,
                      input$theta_geog,
                      polygon_upload_active()) %>%
            bindEvent(input$model_delta_update)

      output$model_delta_plot <- renderPlot({
            d <- model_deltas()

            # Before the first "Update chart" click, model_deltas() is NULL
            # (bindEvent gate). Show an instruction rather than a blank panel.
            if (is.null(d) || nrow(d) == 0) {
                  return(
                        ggplot() +
                              annotate("text", x = 0, y = 0,
                                       label = "Choose a summary region, then click \u201cUpdate chart\u201d.",
                                       color = "gray40", size = 5) +
                              theme_void()
                  )
            }

            # Axes: the per-variable delta columns selected by the dropdowns.
            # Changing axes re-renders from the cached data with no recompute.
            vx <- input$model_delta_x
            vy <- input$model_delta_y
            req(vx %in% names(d), vy %in% names(d))
            d$dx <- d[[vx]]
            d$dy <- d[[vy]]

            # Axis labels use the shared unit_lab() helper; deltas are already
            # back-converted from z-scores in model_deltas.

            # Selected models (current GCM ensemble) are drawn at full opacity;
            # non-selected trajectories are faded to 0.5 so the active ensemble
            # stands out without needing in-plot text labels.
            d$selected <- d$model %in% input$scenario_model

            # Shape palette: one glyph per era, generated to fit however many
            # eras exist (filled square / diamond / circle / triangle, then
            # cycling) rather than a hardcoded 4-value scale.
            eras <- levels(d$era)
            shape_vals <- stats::setNames(
                  rep(c(15, 18, 16, 17), length.out = length(eras)), eras)

            ggplot(d, aes(dx, dy,
                          color = model_label,
                          group = interaction(scenario, model))) +
                  geom_hline(yintercept = 0, color = "gray70", linewidth = .3) +
                  geom_vline(xintercept = 0, color = "gray70", linewidth = .3) +
                  geom_path(aes(alpha = selected)) +
                  geom_point(aes(shape = era, size = era, alpha = selected)) +
                  facet_wrap(~ scenario_label) +
                  scale_shape_manual(values = shape_vals) +
                  scale_size_manual(values = stats::setNames(
                        seq(1.5, 3.5, length.out = length(eras)), eras)) +
                  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5),
                                     guide = "none") +
                  coord_cartesian(clip = "off") +
                  theme(legend.position = "top",
                        legend.direction = "vertical") +
                  strip_style +
                  labs(x = paste0("\u0394 ", vx, " (", unit_lab(vx), ")"),
                       y = paste0("\u0394 ", vy, " (", unit_lab(vy), ")"),
                       color = "Climate model", shape = "Era", size = "Era")
      })


      pool <- reactive({
            req(plot_occ())
            xy <- plot_occ()[, c("x", "y")]
            # Extract climate from the raw (uncropped, unaggregated) rasters at
            # each plot's actual location. This gives the climate the plot is
            # actually experiencing rather than the spatial average over a
            # coarsened cell, and means pool() doesn't depend on the modeling
            # grid (input$clim_agg) or the bbox crop.
            cp <- clim_pool()
            clim_hst <- terra::extract(cp$hst, xy, ID = FALSE, method = "bilinear")
            clim_fut <- terra::extract(cp$fut, xy, ID = FALSE, method = "bilinear")

            # # code used interactively during dev to determine optimal index res
            # tune_index_res(x = cbind(xy, clim_fut), pool = cbind(xy, clim_hst),
            #                env = kernel("gaussian", input$theta_clim, input$theta_clim*3),
            #                geog = kernel("gaussian", input$theta_geog, input$theta_geog*3),
            #                stat = "sum_weights")

            # Response(s) for AIM smoothing. Pooled-ratio variables pack the
            # numerator/denominator/occurrence pieces so aim() can form the
            # surface as occurrence-gated ratio-of-means; everything else uses a
            # single response column.
            rp <- ratio_pieces()
            vals <- if (is.null(rp)) {
                  tibble::tibble(response = response_values(),
                                 samp_method_cd = plot_occ()$samp_method_cd)
            } else {
                  tibble::tibble(occ = as.integer(rp$denom > 0),
                                 num = rp$num,
                                 denom = rp$denom,
                                 samp_method_cd = plot_occ()$samp_method_cd)
            }

            # The same plots, in the same rows, in AIM coordinates (see
            # aim_surfaces()), as a table for the FFT method and an index for
            # the exact one. index_hst/index_fut stay lon/lat for the focal-site
            # analog search.
            aim_pts <- data.frame(aim_xy(xy), clim_hst)

            list(index_hst = build_analog_index(cbind(xy, clim_hst), coord_type = "lonlat",
                                                geog_res_adj = 4),
                 index_fut = build_analog_index(cbind(xy, clim_fut), coord_type = "lonlat",
                                                geog_res_adj = 4),
                 aim_pts = aim_pts,
                 aim_index = build_analog_index(aim_pts, coord_type = "projected",
                                                geog_res_adj = 4),
                 clim_hst = clim_hst,
                 clim_fut = clim_fut,
                 # Active response(s) for AIM smoothing (built as `vals` above),
                 # plus samp_method_cd carried through for plot styling.
                 values = vals)
      })


      ## ANALOG IMPACT MODEL -----------------------------------------

      ### full raster -----------------------
      aim <- reactive({

            is_ratio <- !is.null(ratio_pieces())

            # Stats that describe the WEIGHTS rather than an outcome column.
            # analog_impact() returns one layer for each of these regardless of
            # how many outcome columns `y` carries. "count" and "sum_weights"
            # were requested here but read nowhere in the app; they are dropped.
            # (sum_weights becomes useful again if the multi-species module
            # wants a per-species support signal: the weighted count of
            # contributing presence plots is weighted_mean_<sp> * sum_weights.)
            weight_stats <- "ess"

            # Outcome columns for AIM. pool()$values also carries
            # samp_method_cd, which exists for the site_analogs() / FIA-marker
            # path and is not a response -- interpolating it produced an extra
            # unread layer on every run.
            y_cols <- setdiff(names(pool()$values), "samp_method_cd")

            # analog_impact() names outcome layers `weighted_mean_<column>` when
            # `y` has two or more columns, but bare `weighted_mean` when it has
            # exactly one. Stripping the prefix with gsub() therefore yields
            # "weighted_mean" rather than "response" in the single-column case,
            # leaving aim()$bsl$response NULL. That never fired only because
            # samp_method_cd rode along as a spurious second column -- i.e. the
            # rename was load-bearing on an accident. Rename positionally
            # instead: outcome layers are whatever isn't a weight stat, returned
            # in y-column order.
            name_outcomes <- function(r){
                  is_outcome <- !names(r) %in% weight_stats
                  stopifnot(sum(is_outcome) == length(y_cols))
                  names(r)[is_outcome] <- y_cols
                  r
            }

            surf <- aim_surfaces(list(bsl = clim()$hst, fut = clim()$fut), pool(),
                                 y = pool()$values[, y_cols, drop = FALSE],
                                 theta_clim = input$theta_clim, theta_geog = input$theta_geog,
                                 stat = c("weighted_mean", weight_stats))
            aim_bsl <- name_outcomes(surf$bsl)
            aim_fut <- name_outcomes(surf$fut)

            # Pooled-ratio surface: ratio-of-means from the weighted_mean layers
            # occ/num/denom, combined into the standard `response` layer so all
            # downstream code is unchanged.
            #
            # num/denom is already conditional on occurrence -- absent plots
            # contribute 0 to both and drop out -- so it reads "the ratio among
            # plots where the focal entity occurs".
            #
            # The regeneration fractions stay conditional and are MASKED (NA)
            # where no occupied plot falls within kernel range. Multiplying by
            # occurrence instead would suppress unsuitable ground by driving the
            # value to 0 -- but 0 already means "occurs, adults only" (relict), so
            # the product would conflate two very different states. Masking keeps
            # them distinct: blank = doesn't occur, 0 = occurs but no regen. This
            # also matches the per-plot markers, which are already NA where the
            # species is absent and 0 where it's adults-only.
            #
            # Proportion basal area keeps the occurrence weighting: there 0 means
            # "the focal species is 0% of a real stand" vs "there is no stand" --
            # both read as "the species isn't here", so collapsing them costs
            # little, and the unconditional surface is better behaved.
            if (is_ratio) {
                  conditional <- grepl("fraction", input$fia_var)
                  combine_ratio <- function(r){
                        resp <- r[["num"]] / r[["denom"]]
                        if (!conditional) resp <- r[["occ"]] * resp
                        # denom == 0: no occupied plot in range. Undefined for the
                        # fractions; 0 for the occurrence-weighted proportion.
                        resp[!is.finite(resp[])] <- if (conditional) NA else 0
                        names(resp) <- "response"
                        c(r, resp)
                  }
                  aim_bsl <- combine_ratio(aim_bsl)
                  aim_fut <- combine_ratio(aim_fut)
            }

            list(fut = aim_fut,
                 bsl = aim_bsl)
      })

      ### focal site ------------------------
      site_analogs <- reactive({
            req(site$x, site$y)

            xy <- t(as.matrix(c(x = site$x, y = site$y), nrow = 1))
            # Extract from raw rasters to match pool() (consistent scaling).
            cp <- clim_pool()
            hst <- cbind(xy, terra::extract(cp$hst, xy, method = "bilinear"))
            fut <- cbind(xy, terra::extract(cp$fut, xy, method = "bilinear"))

            max_clim <- if(input$filter_clim) input$theta_clim * 3 else NULL

            # Hoist the pool reactive's result so the inner function doesn't
            # shadow it with its `index` argument.
            pool_data <- pool()
            # Climate columns in clim_hst and clim_fut share the same names
            # (e.g., "aet", "cwd") since the post-refactor rasters have no
            # 1/2 suffix. Tag them explicitly so downstream consumers
            # (clim_plot) can pull historical vs. future without ambiguity.
            clim_hst_tagged <- setNames(pool_data$clim_hst,
                                        paste0(names(pool_data$clim_hst), "_hst"))
            clim_fut_tagged <- setNames(pool_data$clim_fut,
                                        paste0(names(pool_data$clim_fut), "_fut"))
            pool_y <- cbind(pool_data$values, clim_hst_tagged, clim_fut_tagged)

            analog_query <- function(x, index){
                  suppressWarnings(analogs::analog_search(
                        x = x,
                        pool = index,

                        geog = analogs::kernel(max = input$theta_geog * 3),
                        env  = analogs::kernel(max = max_clim),
                        select = "all",
                        stat = NULL,
                        y = pool_y
                  ))
            }

            rev <- analog_query(fut, pool_data$index_hst)
            fwd <- analog_query(hst, pool_data$index_fut)
            bsl <- analog_query(hst, pool_data$index_hst)

            list(hst_site = hst,
                 fut_site = fut,
                 rev_analogs = rev,
                 fwd_analogs = fwd,
                 bsl_analogs = bsl
            )
      })




      ## CONNECTIVITY -------------------------------------------------

      wind <- reactive({
            rast("data/wind_quantized.tif") %>%
                  resample(aim()$bsl[[1]], method = "bilinear", threads = TRUE)
      })

      ### Sources and targets ----------
      electrodes <- reactive({
            if(input$cs_sources == "grid cells"){
                  sources <- aim()$bsl$response
            }else if(input$cs_sources == "selected site" && !is.null(site$y)){
                  xy <- matrix(c(site$x, site$y), nrow = 1)
                  sources <- aim()$bsl$response
                  sources <- rasterize(xy, sources) # * extract(sources, xy)[1,1]
            }else{
                  # fia plots
                  v <- plot_occ()[, c("x", "y")]
                  v$z <- response_values()
                  v <- v[v$z > 0,]
                  v <- vect(v, geom = c("x", "y"))
                  v <- rasterize(v, aim()$bsl$response, field = "z")
                  sources <- v
            }
            sources[is.na(sources[])] <- 0
            sources <- sources / global(sources, "sum")$sum

            if(input$cs_targets == "grid cells"){
                  destinations <- aim()$fut$response
            }else if(input$cs_targets == "selected site" && !is.null(site$y)){
                  xy <- matrix(c(site$x, site$y), nrow = 1)
                  destinations <- aim()$fut$response
                  destinations <- rasterize(xy, destinations)
            }else{
                  # fia plots
                  v <- plot_occ()[, c("x", "y")]
                  v$z <- response_values()
                  v <- v[v$z > 0,]
                  v <- vect(v, geom = c("x", "y"))
                  v <- rasterize(v, aim()$fut$response, field = "z")
                  destinations <- v
            }
            destinations[is.na(destinations[])] <- 0
            destinations <- destinations / global(destinations, "sum")$sum

            c(sources %>% setNames("sources"),
              destinations %>% setNames("destinations"))
      })

      suitability <- reactive({
            suit <- sqrt(aim()$bsl$response * aim()$fut$response)
            suit[is.na(suit)] <- 0
            suit
      })

      ### Resistance --------------
      resistance <- reactive({

            if(grepl("wind", input$cs_conductance_layer)){
                  if(input$cs_sources == "grid cells"){
                        bearing <- gradient_bearing(electrodes()$destinations)
                  }else{
                        # convert points to smooth surface for bearing computations
                        bearing <- electrodes()$destinations
                        bearing[bearing == 0] <- NA
                        bearing <- distance(bearing)
                        bearing <- gradient_bearing(-bearing)
                  }
            }

            if(input$cs_conductance_layer == "wind"){
                  conductance <- aligned_flux(bearing, wind())
                  # note: consider adding a minimum conductance to account for local diffusion
            } else if(input$cs_conductance_layer == "suitability") {
                  conductance <- suitability()
            } else if(input$cs_conductance_layer == "suitability * wind") {
                  conductance <- aligned_flux(bearing, wind()) * suitability()
            } else { # uniform conductance
                  # v <- 10^input$cs_conductance
                  conductance <- setValues(electrodes()$sources, 1)
            }

            mean_cond <- 10 ^ input$cs_conductance_scale
            conductance <- conductance / global(conductance, "mean", na.rm = T)$mean * mean_cond
            conductance <- setNames(conductance, "conductance")
            resistance <- (1 / (conductance + 1e-3)) %>% setNames("resistance")
            c(conductance, resistance)
      })

      ### Loss ---------------------
      leakage <- reactive({
            loss <- switch(input$cs_loss_layer,
                           "none" = setValues(electrodes()$sources, 1e-6),
                           "uniform" = setValues(electrodes()$sources, 1),
                           "mean suitability" = (1 - sqrt(electrodes()$sources * electrodes()$destinations)),
                           "max suitability" = (1 - max(electrodes())),
                           "min suitability" = (1 - min(electrodes()))) %>%
                  setNames("loss")
            mean_loss <- 10 ^ input$cs_loss_scale
            if(input$cs_loss_layer != "none") loss <- loss / global(loss, "mean", na.rm = T)$mean * mean_loss
            loss
      })

      ### Circuitscape --------------
      connectivity <- reactive({
            # Belt and braces: the stat option is withdrawn in community mode,
            # but a stale selection during a flush must not reach a solve.
            req(!multi_mode())

            # First connectivity computation in this R process pays the Julia /
            # Circuitscape warmup; subsequent calls are a no-op. See the
            # ensure_cs_setup() definition in global setup.
            ensure_cs_setup()

            # Ground nodes
            grounds <- (electrodes()$destinations + leakage()) %>% setNames("grounds")

            # Run circuitscape
            fit_cs <- function(correction = 1e-12){
                  cs_advanced(resistance()$resistance, electrodes()$sources,
                              grounds + correction, # add numerical stability for cholmod solver
                              resistance_is = "resistances", ground_is = "conductances",
                              solver = "cholmod") # cholmod is way faster than default
            }
            for(corr in seq(-12, -6)){ # larger factors are safer but introduce more spurious grounds, so start optimistic
                  cs <- try(fit_cs(10 ^ corr))
                  if(class(cs) != "try-error") break
            }

            # Derived variables
            # cs_advanced() returns a SpatRaster with layers 'current' and
            # 'voltage'. (Advanced mode solves a single circuit, so there is
            # no 'cumulative_current' layer -- that name applies only to the
            # multi-focal-node modes.)
            current <- cs$current %>% setNames("current_flow")
            voltage <- cs$voltage
            curr_ground <- grounds %>% "*"(voltage) %>% setNames("current_ground")
            curr_dest <- electrodes()$destinations %>% "*"(voltage) %>% setNames("current_destinations")
            curr_loss <- leakage() %>% "*"(voltage) %>% setNames("current_loss")

            # Power dissipation (pinch points)
            curr_power <- (current^2 * resistance()$resistance) %>% setNames("current_power")
            curr_dir <- gradient_bearing(-voltage) %>% setNames("current_direction")

            # Return results
            c(current, voltage, curr_dest, curr_loss, curr_power,
              resistance()$conductance, leakage(), curr_dir) %>%
                  mask(aim()$bsl$response)
      })



      # MAP LAYERS ---------------------------------------

      ### Basemap -----------------------------------------------
      # No reactive dependencies — renders once at startup. All subsequent
      # map updates go through leafletProxy in dedicated observers.
      output$map <- renderLeaflet({
            leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
                  addProviderTiles(
                        providers$OpenTopoMap,
                        options = providerTileOptions(className = "basemap-gray")
                  ) %>%

                  # zoom to FIA extent
                  fitBounds(lng1 = fia_xy_bounds$xmin,
                            lat1 = fia_xy_bounds$ymin,
                            lng2 = fia_xy_bounds$xmax,
                            lat2 = fia_xy_bounds$ymax) %>%

                  addScaleBar(
                        position = "bottomleft",
                        options = scaleBarOptions(
                              metric = TRUE, # Show kilometers/meters
                              imperial = FALSE, # Hide miles/feet
                              updateWhenIdle = TRUE # Update only after the user stops zooming/panning
                        )
                  )
      })


      ### Focal site -----------------------------
      # Update site coordinates on click only
      observeEvent(input$map_click, {
            site$x <- input$map_click$lng
            site$y <- input$map_click$lat
      })

      # Clear the focal site. Sets coords back to NULL; downstream observers
      # and reactives short-circuit on this, so the map marker, focal-site
      # rings, density-plot vertical line, and analog plots all clear.
      observeEvent(input$site_clear, {
            site$x <- NULL
            site$y <- NULL
      })

      # Draw focal site whenever relevant inputs change
      observe({
            # Dependencies that should trigger redraw
            input$bbox_x
            input$bbox_y
            input$theta_geog
            input$show_site
            site_color()
            site$x
            site$y

            leafletProxy("map") %>%
                  clearGroup("focal_site")

            if(input$show_site && !is.null(site$x)){
                  r1 <- geo_circle(cbind(site$x, site$y), input$theta_geog)
                  r3 <- geo_circle(cbind(site$x, site$y), input$theta_geog * 3)

                  leafletProxy("map") %>%
                        addCircles(lng = site$x, lat = site$y,
                                   radius = input$theta_geog / 10,
                                   opacity = 1, fillOpacity = .5, color = site_color(),
                                   group = "focal_site") %>%
                        addPolygons(lng = r1[,1], lat = r1[,2],
                                    color = site_color(), weight = 2, fill = FALSE,
                                    dashArray = "12, 5",
                                    group = "focal_site") %>%
                        addPolygons(lng = r3[,1], lat = r3[,2],
                                    color = site_color(), weight = 2, fill = FALSE,
                                    group = "focal_site")
            }
      })

      # raster layer of euclidean climatic distance to focal site
      site_diff <- reactive({
            req(site$x, site$y)

            s <- switch(input$analog_direction,
                        "Contemporary analogs" = site_analogs()$hst_site,
                        "Reverse analogs" = site_analogs()$fut_site,
                        "Forward analogs" = site_analogs()$hst_site)

            r <- switch(input$analog_direction,
                        "Contemporary analogs" = clim()$hst,
                        "Reverse analogs" = clim()$hst,
                        "Forward analogs" = clim()$fut)

            for(i in names(r)) r[[i]] <- abs(r[[i]] - s[[i]]) ^ 2
            r <- sqrt(app(r, sum))
            return(r)
      })


      ### Raster -------------------------------
      ## When the point layer is the focal site's analogs, which analogs it shows
      ## depends on analog_direction -- a control in a right-sidebar panel, with
      ## nothing on the map to say which setting is active. This is that line.
      ## NULL when the qualifier does not apply, so callers can drop it.
      analog_vector_note <- reactive({
            if (!isTruthy(input$vector) || input$vector != "focal site analogs") return(NULL)
            if (!isTruthy(input$analog_direction)) return(NULL)
            paste0("(", tolower(input$analog_direction), " of focal site)")
      })

      ## Plain-text form of the legend title, for ggplot axis labels.
      ##
      ## raster_title() emits HTML because the leaflet legend renders it. Callers
      ## used to strip the "<br/>" it once contained; when that was rewritten to
      ## <div> blocks the strippers silently stopped matching and started printing
      ## raw markup as an axis label. Stripping tags generically means the two
      ## cannot drift apart again.
      raster_title_text <- reactive({
            trimws(gsub("\\s+", " ", gsub("<[^>]*>", " ", raster_title())))
      })

      ## "delta" is the internal value; the control and the legend both say
      ## "change", which also covers the paired measures -- they are change
      ## quantities without being deltas of anything.
      time_label <- reactive({
            if (identical(input$time, "delta")) "change" else input$time
      })

      ## What the evaluation is scoring, for axis and legend labels: the
      ## community measure in community mode, the variable otherwise.
      eval_response_label <- reactive({
            if (multi_mode() && isTruthy(input$comm_measure)) comm_label()
            else tolower(input$fia_var)
      })

      ## Display name for the community measure, with the life stage folded in
      ## rather than repeated as a separate variable line. "sapling richness"
      ## says everything "richness (sapling presence)" said, on one line.
      comm_label <- reactive({
            req(isTruthy(input$comm_measure))
            stage <- switch(input$fia_var,
                            "seedling presence" = "seedling ",
                            "sapling presence"  = "sapling ",
                            "tree presence"     = "adult ",
                            "")
            paste0(stage, input$comm_measure)
      })

      ## MULTI-ENTITY RUN AND REDUCTION -----------------
      ##
      ## One analog_impact() call per timeframe with every entity as an outcome
      ## column. Benchmarks: ~0.19 s for one entity and ~2 s for all 143 at the
      ## default domain and 10 km grain, because kernel weights are computed once
      ## and reused across columns. Cost is linear in cells, so this is cheap at
      ## coarse grain and minutes-scale at 1 km -- no button yet, deliberately,
      ## until a cost estimate says one is needed.
      ## Size and rough duration of a community run at the current settings.
      ## Both are read before the run rather than after, which is the point:
      ## the expensive configuration has to be refusable, not merely slow.
      ##
      ## Timing coefficients are measured, not guessed -- about 0.10 s fixed plus
      ## 0.012 s per entity at 31,624 cells, per timeframe, scaling linearly in
      ## cells. Marginal cost per entity climbs somewhat above ~25 entities, so
      ## this runs a little optimistic at large entity counts.
      comm_run_cost <- reactive({
            req(multi_mode())
            n_ent <- ncol(multi_values())
            r <- clim()$hst
            n_cells <- terra::ncell(r)
            gb <- n_cells * n_ent * 8 * 2 / 1024^3
            secs <- (n_cells / 31624) * (0.10 + 0.012 * n_ent) * 2

            # Smallest grain that would fit, so a refusal can say what to do
            # rather than leaving the user to bisect the slider. Cells scale as
            # the inverse square of grain, so gb(g) = gb_now * (agg_now/g)^2.
            agg <- input$clim_agg
            fits <- ceiling(agg / sqrt(comm_max_gb / gb))
            list(entities = n_ent, cells = n_cells, gb = gb, secs = secs,
                 # NA when even the coarsest grain the slider offers is not
                 # enough -- then the bounding box is the thing to change.
                 grain_that_fits = if (is.finite(fits) && fits <= 10) max(fits, 1) else NA_integer_)
      })

      multi_aim <- reactive({
            req(multi_mode())
            y <- multi_values()
            req(ncol(y) > 0)

            # req(), not validate(): this reactive is reached from the map's
            # leaflet observer, and a validation condition raised outside a
            # render context propagates as an app error rather than a message.
            # The refusal is explained by comm_cost_note in the sidebar, which
            # is a render context and can say it properly.
            cost <- comm_run_cost()
            req(cost$gb <= comm_max_gb)

            # sum_weights is retained even though nothing reads it yet: it is one
            # extra layer for the whole run (a weight property, not per-entity),
            # and it is the second half of the per-cell support signal. The
            # weighted count of contributing presence plots for entity s is
            # weighted_mean_s * sum_weights, which is what any reliability
            # filter or rare-entity floor will be built on.
            weight_stats <- "sum_weights"

            # Same positional rename as aim(): outcome layers are whatever isn't
            # a weight stat, returned in y-column order.
            name_outcomes <- function(r){
                  is_outcome <- !names(r) %in% weight_stats
                  stopifnot(sum(is_outcome) == ncol(y))
                  names(r)[is_outcome] <- names(y)
                  r
            }

            surf <- aim_surfaces(list(bsl = clim()$hst, fut = clim()$fut), pool(), y = y,
                                 theta_clim = input$theta_clim, theta_geog = input$theta_geog,
                                 stat = c("weighted_mean", weight_stats))
            lapply(surf, name_outcomes)
      })

      ## Reduce the per-entity stack to one community surface. Reduction reads
      ## from the stack rather than accumulating during the run, which is what
      ## keeps every additional measure cheap -- re-reducing is a re-read, not a
      ## re-run.
      ## One entity -> color lookup for the whole app.
      ##
      ## Colors have to be keyed to entity IDENTITY, not to rank within a plot.
      ## The map ranks entities by cell count and the composition panel ranks
      ## them by value over the selected area, so assigning palette color 1 to
      ## each view's own first entity gives two entirely different mappings for
      ## the same species.
      ##
      ## The canonical order is domain-wide historical abundance, so the entities
      ## most likely to appear in any view get the first and most distinct colors.
      ## Past the palette length colors recycle, which can put two same-colored
      ## entities in one plot -- but only for entities rare enough to be far down
      ## this ordering, and a visible repeat is better than a mapping that
      ## silently disagrees between views.
      comm_entity_colors <- reactive({
            m <- multi_aim()
            ent <- setdiff(names(m$bsl), "sum_weights")
            rank_val <- as.numeric(terra::global(m$bsl[[ent]], "mean", na.rm = TRUE)[, 1])
            ord <- ent[order(rank_val, decreasing = TRUE)]
            stats::setNames(cat_palette_colors(length(ord), last_cat_palette()), ord)
      })

      ## Categorical measures: the argmax raster AND the code/label/color
      ## mapping, produced together.
      ##
      ## The mapping is carried here rather than as terra factor levels because
      ## assigning levels<- renumbers the categories to a sequential index --
      ## terra::cats() then returns 0..n-1 while the per-plot values are still
      ## column indices into the entity set, so almost nothing matches and the
      ## markers all collapse to "other".
      comm_category_result <- reactive({
            req(categorical_layer())
            m <- multi_aim()
            ent <- setdiff(names(m$bsl), "sum_weights")

            if (input$comm_measure %in% paired_measures) {
                  # Biggest winner / loser: argmax and argmin of the per-entity
                  # DIFFERENCE. These cannot be recovered from the reduced
                  # surfaces -- which entity gained most is not a function of how
                  # much the community gained -- so they are the categorical case
                  # that genuinely needs the stack.
                  va <- m$bsl[[ent]]
                  vb <- m$fut[[ent]]
                  dv <- vb - va
                  idx <- if (identical(input$comm_measure, "biggest winner"))
                        terra::which.max(dv) else terra::which.min(dv)

                  # Support over both eras together, so a cell that is empty now
                  # but gains entities still qualifies. Requiring historical
                  # support would mask exactly the arrivals this layer is for.
                  support <- (terra::app(va, fun = "sum", na.rm = TRUE) +
                                    terra::app(vb, fun = "sum", na.rm = TRUE)) / 2
            } else {
                  # Dominant entity: per timeframe.
                  r <- if (identical(input$time, "future")) m$fut else m$bsl
                  v <- r[[ent]]
                  idx <- terra::which.max(v)
                  support <- terra::app(v, fun = "sum", na.rm = TRUE)
            }

            # Support floor FIRST. which.max always names a winner, including in
            # cells where essentially nothing is expected to be present -- so
            # without this the map shows a dominant tree species across
            # unforested country. Masking before the frequency count also keeps
            # those cells out of the top-N ranking, which otherwise promotes
            # whichever entity happens to win the empty ground.
            idx <- apply_support_floor(idx, support, input$support_floor)

            # Lump all but the most frequent N. Ranking by cell count rather
            # than mean value answers "what covers this map", which is the
            # question a discrete palette can address.
            n_top <- if (isTruthy(input$comm_top_n)) input$comm_top_n else 10
            f <- terra::freq(idx)
            f <- f[order(-f$count), , drop = FALSE]
            top <- utils::head(f$value, max(1, n_top))
            other <- setdiff(f$value, top)
            if (length(other)) idx <- terra::subst(idx, other, 0)
            names(idx) <- input$comm_measure

            codes  <- c(if (length(other)) 0, top)
            labels <- c(if (length(other)) "other", ent[top])
            # Looked up by entity name, so these match the composition panel
            # regardless of how the two views happen to rank entities.
            cols <- unname(comm_entity_colors()[labels[labels != "other"]])
            cols <- c(if (length(other)) "gray80", cols)

            list(raster = idx, codes = codes, labels = labels, colors = cols)
      })

      ## COMMUNITY COMPOSITION MODULE -------------------
      ##
      ## The map-level measures reduce across entities; this reduces across
      ## SPACE instead, so entity identity survives. Same stack, other axis --
      ## which is why it costs nothing beyond a global() over layers already
      ## computed.
      comm_profile <- reactive({
            req(multi_mode())
            m <- multi_aim()
            ent <- setdiff(names(m$bsl), "sum_weights")
            scope <- if (isTruthy(input$comm_scope)) input$comm_scope else "domain"

            reduce_space <- function(r) {
                  v <- r[[ent]]

                  # One cell, not an average: the composition AT a place.
                  if (scope == "cell") {
                        req(site$x, site$y)
                        return(as.numeric(terra::extract(v, cbind(site$x, site$y))))
                  }

                  mean_over <- function(g)
                        as.numeric(terra::global(terra::mask(terra::crop(v, g), g),
                                                 "mean", na.rm = TRUE)[, 1])

                  if (scope %in% c("polygon", "buffered")) {
                        req(isTruthy(polygon_upload_active()))
                        return(mean_over(terra::vect(
                              if (scope == "polygon") user_polygon() else user_polygon_buffered())))
                  }

                  if (scope == "site") {
                        req(site$x, site$y, isTruthy(input$theta_geog))
                        return(mean_over(terra::buffer(
                              terra::vect(cbind(site$x, site$y), crs = "EPSG:4326"),
                              width = input$theta_geog * 3 * 1000)))
                  }

                  as.numeric(terra::global(v, "mean", na.rm = TRUE)[, 1])
            }

            hst <- reduce_space(m$bsl)
            fut <- reduce_space(m$fut)
            req(length(hst) == length(ent), length(fut) == length(ent))

            tibble::tibble(entity = ent, hst = hst, fut = fut) %>%
                  # Entities absent in both timeframes are not part of this
                  # community: they contribute zero to every reduction and only
                  # crowd the label spread and the ranking.
                  filter(hst > 0 | fut > 0) %>%
                  mutate(change = fut - hst,
                         # Proportional change is undefined where nothing is
                         # there now, and wildly unstable just above that -- the
                         # same denominator problem the support floor addresses
                         # on the maps. Entities below the floor are dropped from
                         # the proportional ranking rather than topping it.
                         prop = ifelse(hst > 0, change / hst, NA_real_)) %>%
                  filter(!is.na(hst), !is.na(fut))
      })

      ## The entity set the panel names and colors, shared by the alluvial and
      ## the scatterplot so an entity has one color across both and the two show
      ## the same named set. Everything outside it is the "other" pool -- which
      ## is also exactly what the scatterplot's gray points are.
      comm_bands <- reactive({
            d <- comm_profile()
            req(nrow(d) > 0)
            n <- if (isTruthy(input$comm_rank_n)) input$comm_rank_n else 15

            # Rank by the larger of the two values, so an entity arriving from
            # near zero is not excluded for having a small historical one.
            #
            # with_ties = FALSE is essential, not tidying: slice_max() otherwise
            # returns EVERY row tied at the cutoff value, and at a focal site
            # most entities sit at exactly 0. The n-th value is then 0, the tie
            # set is the whole pool, `rest` comes back empty, and nothing is ever
            # pooled into "other".
            d <- mutate(d, top = pmax(hst, fut, na.rm = TRUE))
            keep <- dplyr::slice_max(d, top, n = n, with_ties = FALSE)
            rest <- filter(d, !entity %in% keep$entity)

            # Largest first, so it sits at the bottom of the alluvial column and
            # the eye reads dominant entities against a fixed baseline.
            bands <- keep %>% arrange(desc(top)) %>% select(entity, hst, fut)

            # The tail is pooled rather than dropped: column height in the
            # alluvial is the community sum, and discarding ~128 entities would
            # make that aggregate silently wrong.
            if (nrow(rest))
                  bands <- bind_rows(bands,
                                     tibble::tibble(entity = "other",
                                                    hst = sum(rest$hst, na.rm = TRUE),
                                                    fut = sum(rest$fut, na.rm = TRUE)))

            # Same palette family as the dominant-entity map, so an entity keeps
            # a consistent identity across the map and both panel plots.
            cols <- comm_entity_colors()[bands$entity]
            cols[bands$entity == "other"] <- "gray80"
            names(cols) <- bands$entity

            list(bands = bands, cols = cols,
                 named = bands$entity[bands$entity != "other"])
      })

      ## Scope preconditions, stated rather than left as an empty panel.
      comm_profile_ready <- reactive({
            # Most fundamental precondition first: someone in single-entity mode
            # with the scope on "focal polygon" should be told about the entity
            # setting, not sent off to upload a polygon that would not help.
            if (!multi_mode())
                  return("This panel doesn't work in single-entity mode. Under 'Entity' in the left sidebar, switch from a single species or type to 'community composition'.")
            scope <- if (isTruthy(input$comm_scope)) input$comm_scope else "domain"
            if (scope %in% c("polygon", "buffered") && !isTruthy(polygon_upload_active()))
                  return("No focal polygon yet. Upload one under 'Study area' in the left sidebar, or summarize over the whole bounding box instead.")
            if (scope %in% c("site", "cell") && (is.null(site$x) || is.null(site$y)))
                  return("No focal site yet. Click anywhere on the map to set one.")
            NULL
      })

      output$comm_scatter <- renderPlot({
            shiny::validate(need(is.null(comm_profile_ready()), comm_profile_ready()))
            d <- comm_profile()
            req(nrow(d) > 0)
            cb <- comm_bands()

            # Named entities are the same set the alluvial draws, in the same
            # colors, so the two plots can be read against each other. The gray
            # points are literally the alluvial's "other" band, unpooled.
            d <- mutate(d, band = ifelse(entity %in% cb$named, entity, "other"))
            lab <- filter(d, band != "other") %>% arrange(desc(abs(change)))

            ggplot(d, aes(hst, change)) +
                  geom_hline(yintercept = 0, color = "gray70", linewidth = .3) +
                  geom_point(data = filter(d, band == "other"),
                             color = "gray80", size = 1.5) +
                  geom_point(data = lab, aes(color = entity), size = 2.4) +
                  # check_overlap rather than a repel dependency: rows are
                  # pre-sorted by |change| so the most extreme entities win the
                  # collisions.
                  geom_text(data = lab, aes(label = entity, color = entity),
                            size = 3.1, hjust = -0.12, vjust = 0.5,
                            check_overlap = TRUE) +
                  scale_color_manual(values = cb$cols, guide = "none") +
                  scale_x_sqrt() +
                  ggplot2::expand_limits(x = max(d$hst, na.rm = TRUE) * 1.35) +
                  theme_bw() +
                  labs(x = paste("historical", tolower(input$fia_var)),
                       y = "change",
                       subtitle = names(area_names)[match(
                             if (isTruthy(input$comm_scope)) input$comm_scope else "domain",
                             area_names)])
      })

      output$comm_ranking <- renderPlot({
            shiny::validate(need(is.null(comm_profile_ready()), comm_profile_ready()))
            d <- comm_profile()
            req(nrow(d) > 0)

            prop <- identical(input$comm_rank_by, "proportional")
            d$value <- if (prop) d$prop else d$change
            d <- filter(d, !is.na(value), value != 0)
            shiny::validate(need(nrow(d) > 0, "Nothing changes measurably over this area."))

            n <- if (isTruthy(input$comm_rank_n)) input$comm_rank_n else 15
            # Both tails, not the top of one: the losers are as much the answer
            # as the winners, and showing only gains hides half the turnover.
            # with_ties = FALSE for the same reason as comm_bands(): a tie at
            # the cutoff would return the whole pool.
            d <- bind_rows(dplyr::slice_max(d, value, n = ceiling(n / 2), with_ties = FALSE),
                           dplyr::slice_min(d, value, n = floor(n / 2), with_ties = FALSE)) %>%
                  distinct(entity, .keep_all = TRUE) %>%
                  mutate(entity = factor(entity, levels = entity[order(value)]))

            ggplot(d, aes(value, entity, fill = value > 0)) +
                  geom_col() +
                  geom_vline(xintercept = 0, color = "gray30", linewidth = .3) +
                  scale_fill_manual(values = c(`TRUE` = "dodgerblue4", `FALSE` = "darkred"),
                                    guide = "none") +
                  theme_bw() +
                  labs(y = NULL,
                       x = if (prop) "proportional change" else "change")
      })

      composition_layer <- reactive({
            req(multi_mode(), isTruthy(input$comm_measure))
            if (input$comm_measure %in% categorical_measures)
                  return(comm_category_result()$raster)
            req(input$comm_measure %in% comm_measure_choices(input$track, input$fia_var, input$time))
            m <- multi_aim()

            reduce <- function(r){
                  ent <- setdiff(names(r), "sum_weights")
                  v <- r[[ent]]

                  # Shared pieces. Everything below is layer-wise terra
                  # arithmetic rather than an app() callback over cells, which
                  # matters at 143 layers.
                  s <- terra::app(v, fun = "sum", na.rm = TRUE)

                  # Shannon entropy of the composition, normalized to
                  # proportions q = v / sum(v):
                  #   H = -sum(q log q) = log(s) - sum(v log v) / s
                  # ifel() guards log(0); a zero entity contributes nothing to
                  # the entropy, and 0 * -Inf would otherwise be NaN.
                  shannon <- function(){
                        vlv <- terra::app(terra::ifel(v > 0, v * log(v), 0),
                                          fun = "sum", na.rm = TRUE)
                        log(s) - vlv / s
                  }

                  out <- switch(input$comm_measure,
                                # Expected richness. Summing occurrence
                                # probabilities is deliberate: thresholding at
                                # p > 0.5 would discard the difference between
                                # 0.51 and 0.99, and Sum(p) is the unbiased
                                # estimator if p is calibrated. Whether it IS
                                # calibrated across ~143 mostly-rare entities is
                                # the open question this layer exposes.
                                "richness" = s,

                                # Hill N1 -- the effective number of entities, in
                                # the same units as richness, so the two are
                                # directly comparable. A cell where one species
                                # dominates has low diversity however many rare
                                # species carry small probabilities, which is
                                # exactly the behavior richness lacks.
                                "diversity" = exp(shannon()),

                                # Pielou's J, against the FULL entity set rather
                                # than a count of entities present. With
                                # continuous probabilities almost every entity is
                                # nonzero somewhere, so a "number present" divisor
                                # would be degenerate. Note this makes evenness
                                # depend on the size of the entity set, which
                                # moves with theta_geog and the bbox.
                                "evenness" = shannon() / log(terra::nlyr(v)),

                                # Share of the community held by its largest
                                # member. Bounded 0-1 and needs no normalization.
                                "dominance" = terra::app(v, fun = "max", na.rm = TRUE) / s)
                  req(!is.null(out))

                  # app(na.rm = TRUE) returns 0 where every layer is NA, i.e.
                  # outside analog range -- where a single-entity surface would
                  # correctly be blank. Restore nodata there, or the map grows a
                  # fringe of hard zeros that read as "no species" rather than
                  # "no prediction". Every entity layer shares one NA mask, since
                  # they come from one call with one set of weights and
                  # multi_values() contains no NAs, so masking by any single
                  # entity layer is enough.
                  out <- terra::mask(out, v[[1]])

                  # Support floor, applied per timeframe so a delta is blank
                  # wherever EITHER timeframe is unsupported. `s` is the natural
                  # denominator here: expected richness for an occurrence
                  # variable, total basal area or stems for an abundance one.
                  # Richness is exempt because it IS the denominator.
                  if (!identical(input$comm_measure, "richness"))
                        out <- apply_support_floor(out, s, input$support_floor)

                  names(out) <- input$comm_measure
                  out
            }

            ## Paired reduction: both stacks at once, entity by entity. The
            ## per-entity differences cannot be recovered from the reduced
            ## surfaces -- Sum(max(f-h, 0)) is not a function of Sum(f) and
            ## Sum(h) -- which is what makes these genuinely require the stack
            ## rather than being obtainable from a cheap All-trees variable.
            ##
            ## `a` is the reference composition and `b` the comparison. For
            ## temporal pairing that is historical and future; the same formulas
            ## take a focal site's composition and every cell's for spatial
            ## pairing, which is why they are written in terms of a and b.
            reduce_paired <- function(a, b){
                  ent <- setdiff(names(a), "sum_weights")
                  va <- a[[ent]]
                  vb <- b[[ent]]

                  d      <- vb - va
                  gains  <- terra::app(terra::ifel(d > 0,  d, 0), fun = "sum", na.rm = TRUE)
                  losses <- terra::app(terra::ifel(d < 0, -d, 0), fun = "sum", na.rm = TRUE)
                  sa     <- terra::app(va, fun = "sum", na.rm = TRUE)
                  sb     <- terra::app(vb, fun = "sum", na.rm = TRUE)

                  out <- switch(input$comm_measure,
                                # Bray-Curtis dissimilarity: sum|a-b| over
                                # sum(a+b), which is (gains + losses) over the
                                # combined community. 0 is identical
                                # composition, 1 is no overlap. Threshold-free,
                                # unlike Jaccard, which needs presence/absence
                                # sets and so needs a cutoff.
                                "turnover" = (gains + losses) / (sa + sb),

                                # Absolute expected arrivals and departures, in
                                # the variable's own units. Asymmetric on
                                # purpose: what must arrive is a different
                                # question from what must leave, and for
                                # sourcing or connectivity only one of them is
                                # the relevant demand.
                                "expected gains"  = gains,
                                "expected losses" = losses,

                                # Gains as a share of the future community: the
                                # fraction of what will be there that is not
                                # there now. This is the colonization demand used
                                # to weight connectivity, and it deliberately
                                # ignores losses.
                                "colonization demand" = gains / sb,

                                # Its mirror: losses as a share of the CURRENT
                                # community -- how much of what is here now is
                                # expected to depart. Normalizing by sa rather
                                # than sb is the whole difference, and it is what
                                # makes this the question a manager asks about
                                # the stand in front of them rather than about
                                # the one that replaces it.
                                #
                                # With the two fractions and turnover together,
                                # it is visible WHICH side drives turnover:
                                # Bray-Curtis is their average when the two
                                # communities are similar in size, so a high
                                # turnover with a low loss fraction is
                                # accumulation rather than replacement.
                                "loss fraction" = losses / sa)
                  req(!is.null(out))

                  out <- terra::mask(out, va[[1]])

                  # Ratios only. Turnover is normalized by the combined
                  # community, colonization demand by the future one, so each is
                  # floored against its own denominator. The absolute sums stay
                  # meaningful at low support and are left alone.
                  if (identical(input$comm_measure, "turnover"))
                        out <- apply_support_floor(out, (sa + sb) / 2, input$support_floor)
                  if (identical(input$comm_measure, "colonization demand"))
                        out <- apply_support_floor(out, sb, input$support_floor)
                  if (identical(input$comm_measure, "loss fraction"))
                        out <- apply_support_floor(out, sa, input$support_floor)

                  names(out) <- input$comm_measure
                  out
            }

            if (input$comm_measure %in% paired_measures)
                  return(reduce_paired(m$bsl, m$fut))

            hst <- reduce(m$bsl)
            fut <- reduce(m$fut)
            switch(input$time,
                   "historical" = hst,
                   "future"     = fut,
                   "delta"      = fut - hst)
      })
      ## The support floor's units are the denominator's units, which change with
      ## the variable and with entity scope.
      ##
      ## The value is reset only when the UNITS change, not on every variable
      ## change: moving between two occurrence variables keeps whatever the user
      ## set, while moving from expected species to ft2/acre resets, because
      ## carrying a threshold of 1 across that boundary would mean something
      ## entirely different. Passing value = NULL leaves the field alone.
      support_floor_units <- reactiveVal(NULL)
      observeEvent(list(multi_mode(), input$fia_var, input$track), {
            req(isTruthy(input$fia_var))
            sp <- support_floor_spec(multi_mode(), input$fia_var, input$track)
            units_changed <- !identical(support_floor_units(), sp$label)
            support_floor_units(sp$label)
            updateNumericInput(session, "support_floor",
                               label = sp$label, step = sp$step,
                               value = if (units_changed) sp$default else NULL)
      })

      ## Repopulate community measures for the current track and variable.
      ## Invalid combinations are removed from the menu rather than flagged
      ## after selection. observeEvent so the handler's read of
      ## input$comm_measure is isolated and cannot re-trigger its own write.
      observeEvent(list(input$track, input$fia_var, input$time), {
            req(multi_mode())
            ch <- comm_measure_choices(input$track, input$fia_var, input$time)
            req(length(ch) > 0)
            cur <- if (isTruthy(input$comm_measure)) input$comm_measure else ""
            updateSelectInput(session, "comm_measure", choices = ch,
                              selected = if (cur %in% ch) cur else ch[1])
      })

      ## A modal when the run first goes over the limit.
      ##
      ## The inline notice is not enough on its own: it renders in the FIA
      ## variable panel, while the grain slider that pushes a run over the limit
      ## is in the Study area panel -- and with one accordion panel open at a
      ## time, the notice is suspended at exactly the moment it matters. So the
      ## user coarsens nothing, sees the map quietly go blank, and has no idea
      ## why.
      ##
      ## Fires once per crossing rather than on every slider tick: the flag is
      ## read with isolate() so the observer does not depend on what it writes,
      ## and clears when the settings come back under the limit.
      comm_limit_warned <- reactiveVal(FALSE)

      observe({
            req(multi_mode())
            cost <- comm_run_cost()

            if (cost$gb <= comm_max_gb) {
                  if (isolate(comm_limit_warned())) comm_limit_warned(FALSE)
                  return()
            }
            if (isolate(comm_limit_warned())) return()
            comm_limit_warned(TRUE)

            fix <- cost$grain_that_fits
            showModal(modalDialog(
                  title = "Too large to map",
                  tags$p(sprintf(
                        "This query would model %s species/forest types over %s cells, which would need about %.1f GB of memory, exceeding this app's %.1f GB limit. The operation will not be run at these settings.",
                        format(cost$entities, big.mark = ","),
                        format(cost$cells, big.mark = ","),
                        cost$gb, comm_max_gb)),
                  tags$p(if (is.na(fix))
                        "Even the coarsest grain size will not fit this bounding box. Shrink the box under 'Study area' in the left sidebar."
                        else sprintf("You could increase grain size to %d km or coarser, shrink the bounding box, or choose a single-species variable instead of `community composition`.", fix)),
                  easyClose = TRUE,
                  footer = tagList(
                        # Acting on the advice directly, rather than sending the
                        # user to find a slider in a panel they are not in.
                        if (!is.na(fix))
                              actionButton("comm_apply_grain",
                                           sprintf("Set grain to %d km", fix)),
                        modalButton("Close")
                  )
            ))
      })

      observeEvent(input$comm_apply_grain, {
            cost <- comm_run_cost()
            req(!is.na(cost$grain_that_fits))
            updateSliderInput(session, "clim_agg", value = cost$grain_that_fits)
            removeModal()
      })

      ## What the current settings will cost, shown before it is incurred rather
      ## than only when it is refused. Silent for runs that are quick enough not
      ## to warrant a warning: a notice on every community run would be noise,
      ## and noise is what gets ignored when it matters.
      output$comm_cost_note <- renderUI({
            req(multi_mode())
            cost <- comm_run_cost()
            if (cost$gb <= comm_max_gb && cost$secs < 10) return(NULL)

            if (cost$gb > comm_max_gb) {
                  div(class = "ts-notice",
                      sprintf("This query would model %s species/forest types over %s cells, which would need about %.1f GB of memory, exceeding this app's %.1f GB limit. The operation will not be run at these settings.",
                              format(cost$entities, big.mark = ","),
                              format(cost$cells, big.mark = ","),
                              cost$gb, comm_max_gb,
                              if (is.na(cost$grain_that_fits))
                                    "Even the coarsest grain size will not fit this bounding box; shrink the box under 'Study area' in the left sidebar."
                              else sprintf("Set the grain size to %d km or coarser, or shrink the bounding box, under 'Study area' in the left sidebar.",
                                           cost$grain_that_fits)))
            } else {
                  div(class = "ts-notice",
                      sprintf("Mapping %s entities over %s cells takes roughly %s and about %.1f GB. Coarsening the grain size is the quickest way to speed this up.",
                              format(cost$entities, big.mark = ","),
                              format(cost$cells, big.mark = ","),
                              if (cost$secs < 90) sprintf("%.0f seconds", cost$secs)
                              else sprintf("%.0f minutes", cost$secs / 60),
                              cost$gb))
            }
      })

      ## Only one case remains that the menu cannot express: no measure at all
      ## applies, so there is nothing to select and the reason needs saying.
      output$comm_note <- renderUI({
            req(multi_mode(), isTruthy(input$fia_var))
            req(length(comm_measure_choices(input$track, input$fia_var, input$time)) == 0)

            msg <- if (!input$fia_var %in% multi_vars) {
                  sprintf("'%s' has no per-entity form yet \u2014 the regeneration fractions are pooled ratios, which need three coupled columns per entity. Choose another variable, or switch back to a single entity.",
                          input$fia_var)
            } else {
                  sprintf("Every FIA plot is assigned a single %s, so counting how many occur together at one plot always gives one. Set Timeframe to 'change' for the turnover measures, which compare how the mix of classes shifts over time and do work here.",
                          tolower(input$track))
            }
            div(class = "ts-notice", msg)
      })

      ## Paired measures are their own quantity, not a delta of a per-timeframe
      ## one, so they must not pick up the diverging zero-centered delta palette
      ## even when input$time happens to be "delta". Expected gains and losses
      ## are strictly non-negative; turnover runs 0 (identical) to 1 (no
      ## overlap). None of them center on zero.
      paired_layer <- reactive({
            isTruthy(input$stat) && input$stat == "FIA variable" &&
                  multi_mode() && isTruthy(input$comm_measure) &&
                  input$comm_measure %in% paired_measures
      })

      ## Categorical layers: the value is an entity index, not a quantity.
      ## Everything that assumes continuity has to branch on this -- palette,
      ## reprojection method, scale transform, and the panels that plot value
      ## distributions.
      categorical_layer <- reactive({
            isTruthy(input$stat) && input$stat == "FIA variable" &&
                  multi_mode() && isTruthy(input$comm_measure) &&
                  input$comm_measure %in% categorical_measures
      })

      ## Category codes, labels and colors for the categorical layers. Derived
      ## from the raster's own factor levels so the markers cannot drift from the
      ## surface they sit on: one mapping, two consumers.
      comm_categories <- reactive({
            req(categorical_layer())
            comm_category_result()[c("codes", "labels", "colors")]
      })

      ## Layers bounded on [0, 1] that read better on a diverging palette than a
      ## sequential one, with the domain pinned to the full range so the
      ## midpoint stays put rather than floating with the observed max.
      ##
      ## Diversity is deliberately NOT here: it is exp(H) in effective-entity
      ## units, unbounded, and shares units with richness -- keeping the two on
      ## the same sequential palette is what makes them comparable.
      bounded_ratio_layer <- reactive({
            if (!isTruthy(input$stat) || input$stat != "FIA variable") return(FALSE)
            if (multi_mode()) return(isTruthy(input$comm_measure) &&
                                           input$comm_measure %in% c("evenness", "dominance"))
            isTruthy(input$fia_var) && grepl("fraction", input$fia_var)
      })

      ## Layers whose NAs carry meaning -- a support mask or an entity not found
      ## within range -- rather than merely marking the edge of the analog
      ## domain. Filling these with 0 on a delta map would read as "no change"
      ## instead of "no information", and for a ratio 0 is a real value at one
      ## end of the scale.
      masked_layer <- reactive({
            if (!isTruthy(input$stat) || input$stat != "FIA variable") return(FALSE)
            if (multi_mode()) return(isTruthy(input$comm_measure) &&
                                           !identical(input$comm_measure, "richness"))
            isTruthy(input$fia_var) && grepl("fraction", input$fia_var)
      })


      display_layer <- reactive({
            req(focal_id())

            # In All mode the community reduction takes the place of the
            # per-entity response. The model diagnostics (ESS, focal site
            # similarity, connectivity) are entity-independent and keep working
            # unchanged, so only the "FIA variable" branch is replaced -- and
            # aim() is not required, so All mode does not also pay for a
            # single-entity run it will not display.
            if (multi_mode() && input$stat == "FIA variable") return(composition_layer())

            req(aim())  # only run once we have data

            if(input$stat == "ESS"){
                  r <- switch(input$time,
                              "delta" = 2 / (1/aim()$fut$ess + 1/aim()$bsl$ess), # harmonic mean
                              "historical" = aim()$bsl$ess,
                              "future" = aim()$fut$ess)
            }else if(input$stat == "focal site similarity"){
                  r <- site_diff()
            }else if(input$stat == "connectivity"){
                  st <- cs_stats[cs_stat_names == input$cs_stat]
                  r <- connectivity()[[st]]
            }else{ # input$stat == "FIA variable"
                  # For the pooled-ratio variables, aim() carries the occurrence
                  # layer alongside the response; it is the denominator these
                  # fractions are conditional on. Mask per timeframe, so a delta
                  # is blank wherever either end is unsupported.
                  is_frac <- isTruthy(input$fia_var) && grepl("fraction", input$fia_var)
                  floor_at <- function(tf) {
                        x <- tf$response
                        if (is_frac && !is.null(tf$occ))
                              x <- apply_support_floor(x, tf$occ, input$support_floor)
                        x
                  }
                  r <- switch(input$time,
                              "delta" = floor_at(aim()$fut) - floor_at(aim()$bsl),
                              "historical" = floor_at(aim()$bsl),
                              "future" = floor_at(aim()$fut))
            }

            r
      })

      raster_title <- reactive({
            ## Legend title. Each semantic unit is its own block, so it starts on
            ## a new line but WRAPS internally rather than being split at a
            ## hard-coded point. The old version broke "presence probability"
            ## after the first word regardless of available width; now that the
            ## title sits in a flexible column beside the color bar, forced
            ## breaks fight the layout instead of helping it.
            lines <- character()
            muted <- logical()
            add <- function(txt, dim = FALSE){
                  if (!length(txt) || is.na(txt) || !nzchar(txt)) return(invisible(NULL))
                  lines <<- c(lines, txt)
                  muted <<- c(muted, dim)
            }

            if(input$stat == "ESS"){
                  add(input$stat)

            }else if(input$stat == "focal site similarity"){
                  add(gsub(" analogs", "", input$analog_direction))
                  add("climate distance from focal site")

            }else if(multi_mode() && input$stat == "FIA variable"){
                  add(focal_label())
                  add(time_label())
                  add(comm_label())

                  # Which per-entity variable the measure was computed FROM.
                  # Only richness implies it -- richness is defined over
                  # occurrence and nothing else, and comm_label() already folds
                  # the life stage into its name. Every other measure can be
                  # built from any variable, so "diversity" alone names two
                  # different maps: the diversity of occurrence probabilities and
                  # the diversity of basal area.
                  if (!identical(input$comm_measure, "richness"))
                        add(paste0("of ", tolower(input$fia_var)), dim = TRUE)

                  # The measures that keep the variable's units. Most reduce to a
                  # dimensionless quantity -- diversity and richness are effective
                  # entity counts, and evenness, dominance, turnover and the two
                  # change fractions are ratios in which the units cancel. The
                  # expected gains and losses are sums of per-entity change, so
                  # they are still in whatever the variable measures.
                  if (input$comm_measure %in% c("expected gains", "expected losses"))
                        add(switch(input$fia_var,
                                   "total basal area"    = "(ft2 BA / acre)",
                                   "number of trees"     = "(trees / acre)",
                                   "number of saplings"  = "(saplings / acre)",
                                   "number of seedlings" = "(seedlings / acre)",
                                   # Summed change in occurrence probability is a
                                   # count of entities, like richness.
                                   if (identical(input$track, "Species")) "(species)"
                                   else "(classes)"),
                            dim = TRUE)

            }else{
                  add(focal_label())
                  add(time_label())
                  add(tolower(input$fia_var))
                  add(switch(input$fia_var,
                             "total basal area"    = "(ft2 BA / acre)",
                             "number of trees"     = "(trees / acre)",
                             "number of saplings"  = "(saplings / acre)",
                             "number of seedlings" = "(seedlings / acre)",
                             NULL),
                      dim = TRUE)
                  if(input$stat == "connectivity") add(input$cs_stat)
            }

            # When the layer is categorical the marker legend is suppressed --
            # the points reuse these same entity colors, so a second copy would
            # repeat it verbatim. That leaves this legend serving both layers,
            # which means it also has to say what the points are.
            if (categorical_layer() && !is.null(analog_vector_note())) {
                  add(paste("points:", analog_vector_note()), dim = TRUE)
            }

            # Support qualifier. Stated in the legend rather than only in a
            # tooltip because the legend is what ends up in a screenshot in
            # someone's planning document -- and "expected richness 8" over a
            # 10 km cell will otherwise be read as eight species in 100 km2,
            # which is wrong by a large factor rather than a subtlety.
            #
            # Not applied to the model diagnostics: ESS is a count of plots,
            # focal site similarity is a climate distance, connectivity is a
            # circuit quantity. None are estimated at plot support.
            # In All mode this applies whatever the underlying variable: every
            # community reduction describes the composition of a plot-sized
            # footprint, and diversity in particular rises with area regardless
            # of whether it was built from occurrence or basal area.
            if (input$stat == "FIA variable" &&
                (multi_mode() ||
                 (isTruthy(input$fia_var) && input$fia_var %in% support_dependent_vars))) {
                  add("(per ~0.07 ha plot)", dim = TRUE)
            }

            paste0("<div", ifelse(muted, " style='font-weight:400; opacity:.75;'", ""),
                   ">", lines, "</div>", collapse = "")
      })

      observe({
            r <- display_layer()
            req(r)

            # Scale transforms are meaningless on category codes, and
            # setValues() would drop the factor levels the legend reads from.
            if (!categorical_layer()) {
                  trans <- switch(input$rast_transform,
                                  "identity" = identity,
                                  "log10" = log10,
                                  "sqrt" = sqrt,
                                  "uniform" = function(x) scales::rescale(rank(x, na.last = "keep")))
                  r <- setValues(r, trans(r[]))
                  r[!is.finite(r[])] <- NA
            }

            proxy <- leafletProxy("map") %>%
                  clearImages() %>%
                  removeControl("raster_legend")

            ## palettes ##
            proj_method <- "bilinear"
            mm <- minmax(r)
            dom <- if(input$rast_transform == "identity") c(0, mm[2]) else mm

            choose_palette <- function(x) if(input$rast_palette == "default") x else input$rast_palette

            if(categorical_layer()){
                  # Nearest neighbour: bilinear interpolation of category codes
                  # invents indices that mean nothing. Same reason current
                  # direction uses it.
                  proj_method <- "ngb"
                  cc <- comm_categories()
                  # `levels`, not `domain`: colorFactor() builds its levels with
                  # sort(unique(domain)), so passing the codes as a domain pairs
                  # color k with the k-th SMALLEST entity index while the legend
                  # labels it by frequency rank. Every name shifts. Supplying
                  # levels explicitly keeps the given order.
                  rast_pal <- colorFactor(palette = cc$colors, levels = cc$codes,
                                          na.color = "transparent")
                  dom <- cc$codes
            } else if(input$time == "delta" && input$stat == "FIA variable" && !paired_layer()){
                  # units are deltas
                  r <- display_layer()
                  # Fill non-analog NAs -- except where the NAs are a deliberate
                  # mask (species not found within range, or below the support
                  # floor). Filling those would read as "no change" on the
                  # diverging scale rather than "no information".
                  if(!masked_layer()) r[is.na(r[])] <- 0
                  mm <- minmax(r)
                  dom <- max(abs(mm)) * c(-1, 1)
                  rast_pal <- colorNumeric(
                        palette = c("darkred", "orange", "gray90", "dodgerblue", "darkblue") %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "ESS"){
                  # units are sample sizes
                  rast_pal <- colorNumeric(
                        palette = c("gray90", "darkblue") %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "FIA variable") {
                  # historical or future estimate: domain is unit except for total BA
                  if(bounded_ratio_layer()){
                        # Bounded [0, 1] ratios -- the regeneration fractions
                        # (relict <-> colonizing), and community evenness and
                        # dominance -- default to RdBu. Fixing the domain to the
                        # full [0, 1] range keeps the midpoint put rather than
                        # floating with the observed max, so the same color means
                        # the same value across species and scenarios.
                        if(input$rast_transform == "identity") dom <- c(0, 1)
                        rast_pal <- colorNumeric(
                              palette = "RdBu" %>% choose_palette(),
                              domain = dom,
                              na.color = "transparent"
                        )
                  } else {
                        rast_pal <- colorNumeric(
                              palette = c("gray90", green) %>% choose_palette(),
                              domain = dom,
                              na.color = "transparent"
                        )
                  }
            } else if(input$stat == "connectivity" && input$cs_stat == "current direction") {
                  proj_method <- "ngb" # because bilinear interpolation breaks for a circular variable
                  dom <- if(input$rast_transform == "identity") c(0, 360) else mm
                  rast_pal <- colorNumeric(
                        palette = rainbow(20)[c(1:20, 1)],
                        domain = dom,
                        na.color = "transparent"
                  )
            } else if(input$stat == "focal site similarity") {
                  rast_pal <- colorNumeric(
                        palette = "magma" %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent",
                        reverse = TRUE
                  )
            } else {
                  rast_pal <- colorNumeric(
                        palette = c("gray90", "darkorchid4") %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent"
                  )
            }

            proxy <- proxy %>%
                  addRasterImage(
                        r,
                        colors  = rast_pal,
                        opacity = input$rast_opacity,
                        method = proj_method,
                        group = "AIM results",
                        maxBytes = 10 * 1024 * 1024 # increase limit to 10 MB
                  )

            if(input$rast_opacity > 0){
                  # addLegend() derives its labels from `values` whenever `pal`
                  # is supplied, and silently ignores a `labels` argument -- which
                  # is why the categorical legend has to use the manual
                  # colors/labels form instead. Entity names, not code numbers.
                  proxy <- if (categorical_layer()) {
                        cc <- comm_categories()
                        proxy %>%
                              addLegend(
                                    position = "topright",
                                    colors   = cc$colors,
                                    labels   = cc$labels,
                                    title    = HTML(raster_title()),
                                    layerId  = "raster_legend",
                                    opacity  = input$rast_opacity,
                                    className = "info legend legend-discrete"
                              )
                  } else {
                        proxy %>%
                              addLegend(
                                    position = "topright",
                                    pal      = rast_pal,
                                    values   = dom,
                                    title    = HTML(raster_title()),
                                    layerId  = "raster_legend",
                                    opacity  = input$rast_opacity
                              )
                  }
            }
      })


      ### Vectors ----------------------------
      observe({
            req(plot_occ())

            proxy <- leafletProxy("map")

            ## clear old layers ##
            proxy %>%
                  # clearMarkers() %>%
                  clearGroup("FIA plots") %>%
                  clearGroup("current_flow") %>%
                  clearGroup("current_field") %>%
                  removeControl("point_legend")

            if(input$vector %in% c("all FIA plots", "occupied FIA plots", "focal site analogs")){

                  pd <- vector_points()

                  # Regeneration fractions get the variable's full [0, 1] range so
                  # 0.5 stays on the diverging palette's midpoint. Not in All
                  # mode, where the fractions have no per-entity form and the
                  # displayed value is a count.
                  cat_layer <- categorical_layer()
                  cc <- if (cat_layer) comm_categories() else NULL

                  dom <- if(cat_layer) cc$codes
                  else if(response_is_binary()) 0:1
                  else if(bounded_ratio_layer()) c(0, 1)
                  else seq_domain(pd$display)

                  # Markers follow the raster palette only when the raster is
                  # itself showing FIA response on a sequential scale (see
                  # marker_follows_raster_palette); otherwise the canonical
                  # gray90 -> green ramp, except for the regeneration fractions,
                  # which default to the same diverging RdBu as their raster (on
                  # the sequential ramp their 0 = relict end renders near-white).
                  # colorNumeric accepts either a color vector or a palette name.
                  marker_palette <- if(input$rast_palette != "default" &&
                                       marker_follows_raster_palette()){
                        input$rast_palette
                  } else if(bounded_ratio_layer()){
                        "RdBu"
                  } else {
                        c("gray90", green)
                  }
                  # Categorical layers key on entity codes, so the markers need
                  # the raster's own colorFactor rather than a numeric ramp --
                  # and its breaks are the categories themselves.
                  plot_pal <- if(cat_layer){
                        # levels, not domain -- see the raster palette above.
                        colorFactor(palette = cc$colors, levels = cc$codes)
                  } else {
                        colorNumeric(palette = marker_palette, domain = dom)
                  }
                  # Only feeds the marker legend, which categorical layers do
                  # not draw -- guarded so pretty() is not handed a code vector.
                  breaks <- if(cat_layer) numeric(0)
                  else pretty(seq(dom[1], dom[2], length.out = 100), n = 5)

                  proxy <- proxy %>%
                        addCircles(
                              data = filter(pd, samp_method_cd == 2),
                              lng = ~x,
                              lat = ~y,
                              radius = point_radius()/2, # meters
                              stroke = FALSE, fillColor = "black",
                              group = "FIA plots"
                        ) %>%
                        addCircles(
                              # Drop plots with no defined value: the regeneration
                              # fractions are NA wherever the species is absent as
                              # both the stage and adults, and colorNumeric would
                              # otherwise paint them in its na.color (a mid gray
                              # that reads darker than a true zero). No-op for the
                              # other variables, which are 0 rather than NA when
                              # the focal entity is absent.
                              data = filter(pd, samp_method_cd == 1, !is.na(display)),
                              lng = ~x,
                              lat = ~y,
                              radius = point_radius(), # meters
                              fillColor = ~plot_pal(display),
                              fillOpacity = 1,
                              stroke = TRUE, color = "black", weight = 1,
                              group = "FIA plots"
                        )

                  # legend
                  breaks <- pretty(seq(dom[1], dom[2], length.out = 100), n = 5)
                  legend_labels <- c(breaks, "Non-sampled plot")
                  legend_colors <- c(plot_pal(breaks), "gray") # gray, to match the semi-transparent black markers
                  legend_radii <- c(rep(6, length(breaks)), 3)
                  # Categorical layers: markers use the raster's own colors and
                  # categories, so a second legend would repeat it verbatim. The
                  # raster legend serves both.
                  if (!cat_layer) proxy <- proxy %>%
                        addControl(
                              circleLegend(
                                    # In All mode the markers show observed
                                    # richness, not the per-entity variable.
                                    title  = HTML(paste(c("FIA plots: historical",
                                                          if (multi_mode()) paste("observed", comm_label())
                                                          else input$fia_var,
                                                          analog_vector_note()),
                                                        collapse = "<br/>")),
                                    labels = legend_labels,
                                    colors = legend_colors,
                                    radii = legend_radii
                              ),
                              position = "topright",
                              layerId  = "point_legend",
                              className = "circle-legend leaflet-control"
                        )
            }

            if(input$vector == "current animation"){
                  v_data <- format_velocity_data(connectivity()[["current_flow"]],
                                                 connectivity()[["current_direction"]])
                  proxy <- proxy %>%
                        addVelocity(
                              content = v_data,
                              options = velocityOptions(
                                    colorScale = c("#000000", "#000000"),
                                    lineWidth = 2,
                                    velocityScale = 10^input$vector_scale * 25,
                                    particleMultiplier = 0.005,
                                    displayValues = FALSE,
                                    velocityType = "Current",
                                    emptyString = "No current data"
                              ),
                              group = "current_flow"
                        )
            }

            if(input$vector == "current field"){
                  mag <- connectivity()[["current_flow"]]
                  mag <- mag / minmax(mag)[2]
                  spokes <- make_spokes(mag, connectivity()[["current_direction"]],
                                        subsample = 11 - input$clim_agg,
                                        arrow_scale = 10^input$vector_scale * .1)

                  line_pal <- colorNumeric(palette = c("gray90", "black"), domain = sqrt(spokes$lines$mag))
                  point_pal <- colorNumeric(palette = c("gray90", "black"), domain = sqrt(spokes$points$mag))

                  proxy <- proxy %>%
                        addGlPolylines(data = spokes$lines, color = line_pal(sqrt(spokes$lines$mag)), weight = .25,
                                       group = "current_field") %>%
                        addCircles(
                              data = spokes$points %>% bind_cols(st_coordinates(spokes$points)),
                              lng = ~X, lat = ~Y,
                              radius = point_radius()/2 / 10^input$vector_scale,
                              fillOpacity = 1, stroke = FALSE, fillColor = point_pal(sqrt(spokes$points$mag)),
                              group = "current_field"
                        )
            }

      })


      ## FOCAL SITE PLOTS ----------------------------

      eras <- reactive({
            switch(input$analog_direction,
                   "Reverse analogs" = c("historical", "future"),
                   "Forward analogs" = c("future", "historical"),
                   "Contemporary analogs" = c("historical", "historical"))
      })

      ### climate space --------------------------------
      output$clim_plot <- renderPlot({
            req(site_analogs())
            a <- site_analogs()

            v1 <- input$clim_vars[1]
            v2 <- tail(input$clim_vars, 1) # if only one var selected, v2 == v1
            scl1 <- clim_scales[[v1]]
            scl2 <- clim_scales[[v2]]

            if(input$analog_direction == "Reverse analogs"){
                  d <- a$rev_analogs
                  d$v1 <- d[[paste0(v1, "_hst")]]
                  d$v2 <- d[[paste0(v2, "_hst")]]
                  s <- a$fut_site
            }else if(input$analog_direction == "Forward analogs"){
                  d <- a$fwd_analogs
                  d$v1 <- d[[paste0(v1, "_fut")]]
                  d$v2 <- d[[paste0(v2, "_fut")]]
                  s <- a$hst_site
            }else if(input$analog_direction == "Contemporary analogs"){
                  d <- a$bsl_analogs
                  d$v1 <- d[[paste0(v1, "_hst")]]
                  d$v2 <- d[[paste0(v2, "_hst")]]
                  s <- a$hst_site
            }else{
                  stop("bug")
            }

            d$v1 <- d$v1 * scl1$sd + scl1$mean
            d$v2 <- d$v2 * scl2$sd + scl2$mean

            # display_values(), not response_values(): the analog SET does not
            # depend on the entity -- the same plots are near and climatically
            # similar whichever species is selected -- but the point coloring
            # does, and in community mode it should be the community measure, as
            # on the map markers. Indexed via analog_index because pool()$values
            # has no `response` column for the pooled-ratio variables.
            dv <- display_values()
            shiny::validate(need(
                  !is.null(dv),
                  "This community measure has no per-plot value, so the analog plots cannot be colored by it."))
            req(length(dv) >= max(d$analog_index, na.rm = TRUE))
            d$value <- dv[d$analog_index]
            d <- arrange(d, value) %>% mutate(panel = "Climate space")


            # search perimeters
            perims <- rbind(clim_circle(c(s[[v1]], s[[v2]]), input$theta_clim,
                                        c(scl1$mean, scl2$mean), c(scl1$sd, scl2$sd)) %>%
                                  mutate(r = "b"),
                            clim_circle(c(s[[v1]], s[[v2]]), input$theta_clim * 3,
                                        c(scl1$mean, scl2$mean), c(scl1$sd, scl2$sd)) %>%
                                  mutate(r = "a"))

            fs <- analog_fill(d)
            d <- fs$d

            p <- ggplot() +
                  geom_path(data = perims,
                            aes(v1, v2, group = r, linetype = r),
                            color = site_color(), linewidth = .5, alpha = .5) +
                  geom_point(data = filter(d, samp_method_cd == 2),
                             aes(v1, v2, fill = value),
                             color = "gray70", size = 1.5) +
                  geom_point(data = filter(d, samp_method_cd == 1),
                             aes(v1, v2, fill = value,
                                 shape = paste0("analog FIA plot,\n", eras()[1], " climate")),
                             color = "black", size = 3) +
                  geom_point(aes(x = s[[v1]] * scl1$sd + scl1$mean,
                                 y = s[[v2]] * scl2$sd + scl2$mean),
                             shape = 22,
                             color = site_color(), fill = site_color(), size = 4, alpha = .5) +
                  fs$scale +
                  scale_shape_manual(values = 21) +
                  scale_linetype(guide = "none") +
                  facet_wrap(~ panel) +
                  theme(legend.position = "top") +
                  strip_style +
                  labs(y = v1,
                       x = v2,
                       shape = NULL)

            if(input$analog_direction != "Contemporary analogs"){
                  h <- a$hst_site
                  f <- a$fut_site
                  p <- p +
                        annotate(geom = "segment",
                                 x = h[[v1]] * scl1$sd + scl1$mean,
                                 y = h[[v2]] * scl2$sd + scl2$mean,
                                 xend = f[[v1]] * scl1$sd + scl1$mean,
                                 yend = f[[v2]] * scl2$sd + scl2$mean,
                                 arrow = arrow(angle = 15, type = "closed"),
                                 color = site_color(), alpha = .5)
            }

            if(nrow(filter(d, env_dist <= input$theta_clim * 3)) == 0){
                  p <- p +
                        annotate(geom = "text", color = "red",
                                 label = "No FIA plots found\nwithin geographic & climatic range\nof focal site",
                                 x = mean(perims$v1), y = mean(perims$v2))
            }

            p
      })


      ### distance space --------------------------------


      output$dist_plot <- renderPlot({
            req(site_analogs())

            d <- switch(input$analog_direction,
                        "Reverse analogs" = site_analogs()$rev_analogs,
                        "Forward analogs" = site_analogs()$fwd_analogs,
                        "Contemporary analogs" = site_analogs()$bsl_analogs)
            # display_values(), not response_values(): the analog SET does not
            # depend on the entity -- the same plots are near and climatically
            # similar whichever species is selected -- but the point coloring
            # does, and in community mode it should be the community measure, as
            # on the map markers. Indexed via analog_index because pool()$values
            # has no `response` column for the pooled-ratio variables.
            dv <- display_values()
            shiny::validate(need(
                  !is.null(dv),
                  "This community measure has no per-plot value, so the analog plots cannot be colored by it."))
            req(length(dv) >= max(d$analog_index, na.rm = TRUE))
            d$value <- dv[d$analog_index]
            d <- arrange(d, value) %>% mutate(panel = "Distance space")

            fs <- analog_fill(d)
            d <- fs$d

            ymax <- if(input$filter_clim) input$theta_clim * 3.3 else NA

            # Both axes here ARE the two distances, so the weight is a function
            # of position and can be drawn as a field. That says what the scheme
            # is regardless of which points are shown -- point size alone goes
            # silent when the displayed plots all fall outside the cutoff, which
            # is exactly when someone has turned the filter off to look at them.
            #
            xmax_g <- input$theta_geog * 3.3
            ymax_g <- if (is.na(ymax)) max(d$env_dist, na.rm = TRUE) else ymax
            wfield <- expand.grid(
                  geog_dist = seq(0, xmax_g, length.out = 140),
                  env_dist  = seq(0, ymax_g, length.out = 140))
            wfield$w <- analog_weight(wfield$geog_dist, wfield$env_dist)

            ggplot(d, aes(geog_dist, env_dist, fill = value)) +
                  # Now the only representation of the weighting, so it spans a
                  # wide value range rather than a light tint. Stopping short of
                  # pure black and pure white at the ends keeps the black-outlined
                  # points legible over the darkest corner and leaves the panel
                  # background distinguishable from a zero-weight region.
                  #
                  # alpha rather than fill, because fill already carries the
                  # response value and a second fill scale would need another
                  # package.
                  geom_raster(data = wfield, aes(geog_dist, env_dist, alpha = w),
                              fill = "black", inherit.aes = FALSE) +
                  scale_alpha_continuous(range = c(0.05, 0.95), limits = c(0, 1),
                                         guide = "none") +
                  geom_vline(xintercept = input$theta_geog, linetype = "longdash", color = site_color(), alpha = .5) +
                  geom_vline(xintercept = input$theta_geog * 3, color = site_color(), alpha = .5) +
                  geom_hline(yintercept = input$theta_clim, linetype = "longdash", color = site_color(), alpha = .5) +
                  geom_hline(yintercept = input$theta_clim * 3, color = site_color(), alpha = .5) +
                  scale_x_continuous(expand = c(0, 0), limits = c(0, input$theta_geog * 3.3)) +
                  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +

                  geom_point(data = filter(d, samp_method_cd == 2), color = "gray70", size = 1.5) +
                  geom_point(data = filter(d, samp_method_cd == 1), shape = 21, color = "black", size = 3) +

                  geom_point(data = tibble::tibble(geog_dist = 0, env_dist = 0,
                                                   shape_label = paste0("selected site,\n",
                                                                        eras()[2], " climate")),
                             aes(x = geog_dist, y = env_dist, shape = shape_label),
                             color = site_color(), fill = site_color(),
                             size = 10, alpha = .5, inherit.aes = FALSE) +

                  fs$scale +
                  scale_shape_manual(values = 22) +
                  facet_wrap(~ panel) +
                  theme(legend.position = "top") +
                  strip_style +
                  labs(y = "Climate difference from focal site (z-score)",
                       x = "Geographic distance to focal site (km)",
                       shape = NULL)
      })

      ### forest composition --------------------------------
      ## Composition profile: what the community IS, historical and future.
      ##
      ## Informationally this is the same two numbers per entity as the
      ## scatterplot -- (hst, fut) versus (hst, change) -- but the encoding
      ## carries something neither of the plots above does: because the bands are
      ## stacked, the total height of each column is the COMMUNITY SUM. So
      ## per-entity change and aggregate change are legible in one figure, and a
      ## community that reshuffles without changing size looks different from one
      ## that shrinks.
      ##
      ## Built from the model stack rather than from the analog plots, so it
      ## follows the "Summarize over" control like the rest of the panel. The
      ## previous version was site-only and silently ignored it.
      ##
      ## Forward analogs are gone: "where does this site's climate go" is a
      ## question about the analog machinery, not about composition, and it has
      ## no counterpart in the AIM surfaces.
      output$comm_composition <- renderPlot({
            shiny::validate(need(is.null(comm_profile_ready()), comm_profile_ready()))
            cb <- comm_bands()
            cols <- cb$cols

            # Stack positions; "other" is last in the band order so it ends up
            # on top rather than displacing the dominant entities' baseline.
            bands <- cb$bands %>%
                  mutate(hst_top = cumsum(hst), hst_bot = hst_top - hst,
                         fut_top = cumsum(fut), fut_bot = fut_top - fut)

            # Smoothstep between the two columns, which is what gives the
            # ribbons their shape. Straight edges read as a chart of two bars;
            # the curve reads as flow.
            tt <- seq(0, 1, length.out = 60)
            sm <- tt * tt * (3 - 2 * tt)
            rib <- purrr::map_dfr(seq_len(nrow(bands)), function(i) {
                  b <- bands[i, ]
                  tibble::tibble(entity = b$entity, t = tt,
                                 ymin = b$hst_bot + (b$fut_bot - b$hst_bot) * sm,
                                 ymax = b$hst_top + (b$fut_top - b$hst_top) * sm)
            })
            rib$entity <- factor(rib$entity, levels = bands$entity)

            # Label side. Left is the default because the left column IS the
            # current-day composition, and reading it as a list is the
            # first-order thing this plot is for. Entities that grow move to the
            # right, where their band is thicker and the label has room -- an
            # arriving entity has almost no historical band to sit against.
            lab <- bands %>%
                  mutate(side = ifelse(fut > hst, "right", "left"),
                         y_anchor = ifelse(side == "right",
                                           (fut_bot + fut_top) / 2,
                                           (hst_bot + hst_top) / 2))

            # Every entity gets a label, so they have to be pushed apart rather
            # than dropped. One pass per side: sort by anchor, then walk upward
            # enforcing a minimum gap. The gap is sized so the labels always fit
            # in the panel, and a leader line ties each back to its band.
            y_max <- max(bands$hst_top, bands$fut_top, 0, na.rm = TRUE)
            spread <- function(y) {
                  if (!length(y)) return(y)
                  gap <- y_max / max(26, length(y) + 1)
                  o <- order(y)
                  ys <- y[o]
                  for (i in seq_along(ys)[-1]) ys[i] <- max(ys[i], ys[i - 1] + gap)
                  # A long run can push the top label past the panel; slide the
                  # whole stack down by the overshoot rather than let it clip.
                  over <- max(0, ys[length(ys)] - y_max)
                  ys <- ys - over
                  out <- numeric(length(y))
                  out[o] <- ys
                  out
            }
            lab <- lab %>%
                  group_by(side) %>%
                  mutate(y_label = spread(y_anchor)) %>%
                  ungroup() %>%
                  mutate(x_anchor = ifelse(side == "right", 1, 0),
                         x_label  = ifelse(side == "right", 1.03, -0.03),
                         hjust    = ifelse(side == "right", 0, 1))

            ggplot(rib, aes(t, ymin = ymin, ymax = ymax, fill = entity, group = entity)) +
                  geom_ribbon(color = "white", linewidth = .2) +
                  # Leader lines, in the band's own color, so a displaced label
                  # is still unambiguously attached to its band.
                  geom_segment(data = lab, inherit.aes = FALSE,
                               aes(x = x_anchor, xend = x_label,
                                   y = y_anchor, yend = y_label, color = entity),
                               linewidth = .3, alpha = .6) +
                  geom_text(data = lab, inherit.aes = FALSE,
                            aes(x = x_label, y = y_label, label = entity,
                                hjust = hjust, color = entity),
                            size = 3.1) +
                  scale_fill_manual(values = cols, guide = "none") +
                  scale_color_manual(values = cols, guide = "none") +
                  scale_x_continuous(breaks = c(0, 1),
                                     labels = c("historical", "future"),
                                     expand = ggplot2::expansion(mult = c(0.42, 0.42))) +
                  theme_bw() +
                  theme(panel.grid.major.x = element_blank(),
                        panel.grid.minor = element_blank()) +
                  labs(x = NULL, y = tolower(input$fia_var))
      })


      ## POLYGON RASTER VALUES ------------------------

      # Equal-area CRS for buffering in meters (NAD83 / Conus Albers).
      # Adequate for any CONUS-scale buffer; would need to change for
      # other regions.
      .buffer_crs <- "EPSG:5070"

      # Polygon buffered by theta_geog * 3 (km), the analog search's geographic
      # cutoff. Features are unioned before buffering so multi-feature polygons
      # produce a single dissolved buffer rather than per-feature buffers.
      # Returned in WGS84 to match the user polygon's CRS so map and raster
      # operations are consistent.
      user_polygon_buffered <- reactive({
            req(user_polygon(), isTruthy(input$theta_geog))
            poly <- user_polygon()
            st_transform(poly, .buffer_crs) %>%
                  st_union() %>%
                  st_buffer(dist = input$theta_geog * 3 * 1000) %>%   # km -> m
                  st_transform(4326)
      })

      ## Raster values within the focal site's analog neighborhood. The site was
      ## the only area this module treated as a single cell rather than a
      ## distribution, which made it the odd one out against every other module
      ## and gave it a vertical line where the others get a curve. The single
      ## cell is still marked, so "where my site sits within its neighborhood"
      ## stays readable.
      site_raster_values <- reactive({
            req(!is.null(site$x), !is.null(site$y), display_layer(),
                isTruthy(input$theta_geog))
            b <- terra::buffer(terra::vect(cbind(site$x, site$y), crs = "EPSG:4326"),
                               width = input$theta_geog * 3 * 1000)
            terra::extract(display_layer(), b, na.rm = TRUE) %>%
                  pull(2) %>%
                  na.omit()
      })

      polygon_raster_values <- reactive({
            req(user_polygon(), display_layer())

            poly <- user_polygon()
            r <- display_layer()

            # reproject polygon to match raster CRS
            poly_proj <- st_transform(poly, crs(r))

            terra::extract(r, vect(poly_proj), na.rm = TRUE) %>%
                  pull(2) %>%  # column 1 is ID, column 2 is values
                  na.omit()
      })

      polygon_buffered_raster_values <- reactive({
            req(user_polygon_buffered(), display_layer())

            poly <- user_polygon_buffered()
            r <- display_layer()

            poly_proj <- st_transform(poly, crs(r))

            terra::extract(r, vect(poly_proj), na.rm = TRUE) %>%
                  pull(2) %>%
                  na.omit()
      })

      output$polygon_density <- renderPlot({
            req(display_layer())

            # A categorical layer has no distribution to draw a density of --
            # entity codes are labels, not magnitudes. The question transposes
            # instead: what share of each area does each entity account for.
            # Same areas, same colors, a composition bar apiece.
            if (categorical_layer()) {
                  cc <- comm_categories()
                  share <- function(v) {
                        v <- v[!is.na(v)]
                        if (!length(v)) return(NULL)
                        tibble::tibble(label = cc$labels[match(v, cc$codes)]) %>%
                              filter(!is.na(label)) %>%
                              count(label) %>%
                              mutate(share = n / sum(n))
                  }

                  areas <- list("Model bounding box" = na.omit(values(display_layer(), mat = FALSE)))
                  if (isTruthy(polygon_upload_active())) {
                        areas[["Buffered polygon"]] <- polygon_buffered_raster_values()
                        areas[["Focal polygon"]]    <- polygon_raster_values()
                  }
                  if (!is.null(site$x) && !is.null(site$y)) {
                        areas[["Focal site neighborhood"]] <- site_raster_values()
                        # The single clicked cell, as its own bar. It is one
                        # entity so the bar is solid, which is the point: it says
                        # which entity holds the site, in the same colors as
                        # every other row.
                        areas[["Focal site"]] <- terra::extract(
                              display_layer(), cbind(site$x, site$y))[1, 1]
                  }

                  db <- bind_rows(lapply(names(areas), function(nm) {
                        sh <- share(areas[[nm]])
                        if (is.null(sh)) NULL else mutate(sh, source = nm)
                  }))
                  shiny::validate(need(nrow(db) > 0, "No values to summarize over these areas."))

                  db$source <- factor(db$source, levels = names(areas))
                  # Vertical bars, so the stacking order runs down the page in
                  # the same direction as the legend beside it -- the eye can
                  # match a band to its key entry by height alone. Reversed
                  # levels put the first legend entry at the top of the stack.
                  db$label <- factor(db$label, levels = rev(cc$labels))

                  return(
                        ggplot(db, aes(source, share, fill = label)) +
                              geom_col(width = 0.7) +
                              scale_fill_manual(values = stats::setNames(cc$colors, cc$labels),
                                                breaks = cc$labels, name = NULL) +
                              scale_y_continuous(labels = scales::percent, expand = expansion(0, 0)) +
                              theme_bw() +
                              theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
                              labs(x = NULL, y = "share of area",
                                   subtitle = raster_title_text())
                  )
            }

            domain_vals <- na.omit(values(display_layer(), mat = FALSE))

            # Polygon-dependent values: NULL when no polygon is uploaded. The
            # corresponding curves (focal polygon, buffer) are then simply absent
            # from the plot. Check the upload input directly to avoid
            # triggering user_polygon()'s internal req() short-circuit.
            has_polygon <- isTruthy(polygon_upload_active())
            poly_vals   <- if (has_polygon) polygon_raster_values() else NULL
            buffer_vals <- if (has_polygon) polygon_buffered_raster_values() else NULL

            has_site  <- !is.null(site$x) && !is.null(site$y)
            site_vals <- if (has_site) site_raster_values() else NULL

            # Follow the user's overlay colors, so a curve here and the boundary
            # it describes on the map are recognizably the same thing.
            study_color   <- polygon_color()
            buffer_color  <- lighten_color(polygon_color())
            domain_color  <- "gray60"
            # Follows the user's focal-site color, so the marker here and the
            # marker on the map are the same thing to look at.
            focal_color   <- site_color()

            # Focal site value: extract from the displayed raster at the click
            # location. NA if no site selected or if the point falls on a
            # non-modeled cell.
            focal_val <- if (!is.null(site$x) && !is.null(site$y)) {
                  v <- terra::extract(display_layer(),
                                      cbind(site$x, site$y),
                                      method = "bilinear")
                  v <- v[1, 1]
                  if (is.finite(v)) v else NA_real_
            } else {
                  NA_real_
            }

            # 'uniform' rank-transforms all three distributions against the
            # modeling-domain ECDF, so the three curves remain comparable on a
            # common [0, 1] axis. Done in data space (ggplot scale_x_*
            # transforms operate on the visual layer only, which would distort
            # density estimates by squishing them after they're computed).
            x_lab <- raster_title_text()
            if (input$density_transform == "uniform") {
                  ecdf_ref <- ecdf(domain_vals)
                  if (has_polygon) {
                        poly_vals   <- ecdf_ref(poly_vals)
                        buffer_vals <- ecdf_ref(buffer_vals)
                  }
                  if (has_site && length(site_vals)) site_vals <- ecdf_ref(site_vals)
                  domain_vals <- ecdf_ref(domain_vals)
                  if (!is.na(focal_val)) focal_val <- ecdf_ref(focal_val)
                  x_lab <- paste0(x_lab, " (uniform rank)")
                  d_bounds <- c(0, 1)
            } else {
                  d_bounds <- display_bounds()
            }

            # Build the long data frame, including only the curves we have data
            # for.
            parts <- list(
                  tibble::tibble(value = domain_vals, source = "Model bounding box")
            )
            if (has_polygon) {
                  parts <- c(parts, list(
                        tibble::tibble(value = buffer_vals, source = "Buffered polygon"),
                        tibble::tibble(value = poly_vals,   source = "Focal polygon")
                  ))
            }
            if (has_site && length(site_vals)) {
                  parts <- c(parts, list(
                        tibble::tibble(value = site_vals, source = "Focal site neighborhood")
                  ))
            }
            d <- bind_rows(parts) %>%
                  mutate(source = factor(source, levels = c("Model bounding box",
                                                            "Buffered polygon",
                                                            "Focal polygon",
                                                            "Focal site neighborhood",
                                                            "Focal site")))

            pal <- c("Model bounding box"      = domain_color,
                     "Buffered polygon"        = buffer_color,
                     "Focal polygon"           = study_color,
                     "Focal site neighborhood" = lighten_color(focal_color),
                     "Focal site"              = focal_color)
            ltypes <- c("Model bounding box"      = "solid",
                        "Buffered polygon"        = "dashed",
                        "Focal polygon"           = "solid",
                        "Focal site neighborhood" = "dashed",
                        "Focal site"              = "solid")

            # Optional focal-site reference line. Skipped (NULL → no-op when
            # added to ggplot) if no focal site is selected or the cell value
            # isn't finite. Built as a data-driven geom_vline so the color and
            # linetype scales pick it up for the legend.
            focal_layer <- if (!is.na(focal_val)) {
                  geom_vline(data = tibble::tibble(xint = focal_val,
                                                   source = factor("Focal site",
                                                                   levels = c("Model bounding box",
                                                                              "Buffered polygon",
                                                                              "Focal polygon",
                                                                              "Focal site neighborhood",
                                                                              "Focal site"))),
                             aes(xintercept = xint, color = source,
                                 linetype = source),
                             linewidth = 0.6,
                             alpha = 0.7,
                             inherit.aes = FALSE)
            } else {
                  NULL
            }

            # Legend break order: include Focal site only when one is active.
            legend_breaks <- c("Model bounding box", "Buffered polygon", "Focal polygon")
            if (has_site && length(site_vals))
                  legend_breaks <- c(legend_breaks, "Focal site neighborhood")
            if (!is.na(focal_val)) legend_breaks <- c(legend_breaks, "Focal site")

            # ECDF only. The paired density and ECDF panels showed the same
            # comparison twice, and the density was the weaker of the two here:
            # what the module is for is reading off how one area sits relative to
            # another at a given value, which is the cumulative curve. Two panels
            # also halved the width available for either.
            p <- ggplot(d, aes(x = value, color = source, linetype = source)) +
                  stat_ecdf(aes(group = source),
                            geom = "step", linewidth = 0.7, pad = FALSE) +
                  focal_layer +
                  scale_y_continuous(
                        breaks = c(0, 0.25, 0.5, 0.75, 1),
                        expand = expansion(0, 0)
                  ) +
                  scale_color_manual(values = pal,
                                     breaks = legend_breaks) +
                  scale_linetype_manual(values = ltypes,
                                        breaks = legend_breaks) +
                  theme(legend.position = "top") +
                  labs(x = x_lab,
                       y = "cumulative share of area",
                       color = NULL, linetype = NULL)

            # X-axis limits from the active variable's natural bounds. Infinite
            # bounds become NA so ggplot auto-scales that side from data.
            x_limits <- ifelse(is.finite(d_bounds), d_bounds, NA_real_)

            if (input$density_transform == "sqrt") {
                  p <- p + scale_x_sqrt(limits = x_limits,
                                        expand = expansion(0, 0))
            } else if (input$density_transform == "log10") {
                  # log10 lower limit can't be 0 (undefined); let data set it.
                  log_limits <- x_limits
                  if (!is.na(log_limits[1]) && log_limits[1] <= 0) log_limits[1] <- NA
                  p <- p + scale_x_log10(limits = log_limits,
                                         expand = expansion(0, 0))
            } else {
                  p <- p + scale_x_continuous(limits = x_limits,
                                              expand = expansion(0, 0))
            }

            p
      })


      ## MODEL EVALUATION ----------------------------

      ## What the evaluation cross-validates.
      ##
      ## In single-entity mode this is just the focal entity's response. In
      ## community mode it is the plot-level counterpart of the selected
      ## community measure -- observed richness at each plot, for instance --
      ## so the module scores the surface the map is actually showing rather
      ## than a species the user cannot even see.
      ##
      ## Only linear measures can be handled this way, and richness is the one
      ## that qualifies. Because AIM is a weighted mean with entity-independent
      ## weights, reducing then interpolating equals interpolating then
      ## reducing, so the held-out prediction of plot richness IS the held-out
      ## Sum(p_s). A non-linear measure would need each entity predicted
      ## separately at each held-out plot, and analog_cv() scores one response
      ## at a time -- 143 CV runs, not the cheap multi-outcome trick.
      eval_response <- reactive({
            if (!multi_mode()) return(response_values())
            req(identical(input$comm_measure, "richness"))
            req(isTruthy(input$fia_var), input$fia_var %in% multi_vars)
            comm_reduce_rows(multi_values(), "richness")
      })

      ## Measures the evaluation cannot score, and why. Returned as a message so
      ## the panel explains itself rather than rendering empty.
      eval_unsupported <- reactive({
            if (!multi_mode()) return(NULL)
            if (identical(input$comm_measure, "richness")) return(NULL)
            "Accuracy can only be measured for richness. Switch the community measure to richness, or switch back to a single species."
      })

      model_eval <- reactiveValues(comparison = NULL,
                                   scatter = NULL,
                                   stamp = NULL,
                                   swept = NULL)

      # Message shown when a scope is selected whose defining input isn't set.
      # The scope list is static rather than rebuilt as inputs appear: a
      # reactive choice list makes setting a focal site rebuild the selectize
      # control, whose round-trip re-emits input$eval_scope and re-renders this
      # section even under a scope that ignores the site. Stating the
      # precondition is also more useful than hiding the option.
      eval_scope_missing_msg <- function(scope){
            switch(scope,
                   "polygon"  = "Upload a focal polygon under 'Study area' in the left sidebar to evaluate over it.",
                   "buffered" = "Upload a focal polygon under 'Study area' in the left sidebar to evaluate over its buffered neighborhood.",
                   "site"     = "Click the map to set a focal site, then this scope evaluates the plots within its geographic cutoff.",
                   "This evaluation scope is not available.")
      }

      # Stamp components. Plain functions rather than reactives, so a caller
      # takes a dependency on the polygon or the site only when it actually
      # calls for that piece. Bundling them into one reactive made every map
      # click re-render the heatmap and its warning banner, even under a scope
      # that ignores the site entirely.
      eval_polygon_id <- function(){
            f <- polygon_upload_active()
            if (isTruthy(f)) paste(f$name, f$size, collapse = "|") else NA_character_
      }
      eval_site_id <- function(){
            if (is.null(site$x) || is.null(site$y)) NA_character_
            else sprintf("%.6f,%.6f", site$x, site$y)
      }

      # Settings that define the response and the analog pool, for every scope.
      eval_model_stamp <- function(){
            list(track      = input$track,
                 focal      = focal_id(),
                 # focal_id() already distinguishes community mode from a single
                 # entity, but not which measure is being scored.
                 measure    = if (multi_mode()) input$comm_measure else NA_character_,
                 fia_var    = input$fia_var,
                 bbox       = c(input$bbox_x, input$bbox_y),
                 clim_vars  = paste(sort(input$clim_vars), collapse = "_"),
                 theta_clim = input$theta_clim,
                 theta_geog = input$theta_geog)
      }

      # Full stamp, recorded when a run completes.
      eval_run_stamp <- function(){
            c(eval_model_stamp(),
              list(polygon = eval_polygon_id(), site = eval_site_id()))
      }

      # Stamped fields that have changed since the stored run, with two
      # exemptions applied.
      #
      # SWEPT DIMENSIONS. A bandwidth or variable-set sweep varies those inputs
      # by design, and "Apply optimal parameters" then moves them deliberately.
      # Counting them as staleness would make the heatmap mark itself stale the
      # instant the user acts on its own recommendation.
      #
      # SCOPE RELEVANCE. Handled by only assembling the fields the current scope
      # depends on: stamp_diff() intersects names, so a field left out of `now`
      # is simply not compared. theta_geog is the one field with two roles -- it
      # is a swept kernel parameter, but it also sets the buffer and neighborhood
      # radii, so under those scopes it is never exempt even when swept.
      eval_stale_fields <- reactive({
            req(model_eval$comparison)
            scope <- if (isTruthy(input$eval_scope)) input$eval_scope else "domain"

            exempt <- if (is.null(model_eval$swept)) character() else model_eval$swept
            if (scope %in% c("buffered", "site")) exempt <- setdiff(exempt, "theta_geog")

            now <- eval_model_stamp()
            if (scope %in% c("polygon", "buffered")) now$polygon <- eval_polygon_id()
            if (scope == "site") now$site <- eval_site_id()

            stamp_diff(now, model_eval$stamp, exempt)
      })

      # Logical membership vector per scope, for a set of pool plot coordinates.
      # `xy` must be the POST-filter coordinates so the vectors align with the
      # rows analog_cv() returns.
      #
      # `which` limits which scopes are computed, and that is load-bearing rather
      # than an optimization: reading site$x or the polygon here takes a reactive
      # dependency on it, so computing every scope unconditionally makes a map
      # click invalidate the scatterplot even when the domain scope is selected.
      # The `&&` short-circuits before touching those reactives.
      eval_scope_membership <- function(xy, which = unname(eval_scope_names)){
            keep <- list()
            if ("domain" %in% which) keep$domain <- rep(TRUE, nrow(xy))

            # Check the upload input directly rather than calling user_polygon(),
            # whose internal req() would silently cancel the whole evaluation
            # when no polygon is loaded. Same guard as polygon_density().
            if (any(c("polygon", "buffered") %in% which) &&
                isTruthy(polygon_upload_active())) {
                  pts <- st_as_sf(as.data.frame(xy), coords = c("x", "y"), crs = 4326)
                  if ("polygon" %in% which)
                        keep$polygon <- lengths(st_intersects(pts, st_union(user_polygon()))) > 0
                  if ("buffered" %in% which)
                        keep$buffered <- lengths(st_intersects(pts, user_polygon_buffered())) > 0
            }
            if ("site" %in% which && !is.null(site$x) && !is.null(site$y)) {
                  dkm <- geosphere::distHaversine(cbind(site$x, site$y),
                                                  as.matrix(xy[, c("x", "y")])) / 1000
                  keep$site <- dkm <= input$theta_geog * 3
            }
            keep
      }


      eval_var_set <- function(clim_vars, vary_bandwidths = TRUE){
            # Extract historical climate from the raw raster at each pool plot,
            # matching pool().
            clim_r <- rast("data/clim_historical.tif")[[clim_vars]]

            # Build pool (x/y + climate) and response from a SINGLE plot_occ()
            # snapshot so the two are row-aligned. Capture plot_occ() once and
            # derive the response from that same object rather than re-deriving
            # it, which could in principle see a different value.
            po <- plot_occ()
            xy <- as.data.frame(po[, c("x", "y")])
            clim_hst <- terra::extract(clim_r, xy, ID = FALSE, method = "bilinear")
            d <- cbind(xy, clim_hst)
            values <- eval_response()
            req(!is.null(values), length(values) == nrow(po))

            # analog_cv() coerces `pool` internally and requires `y` to have
            # exactly nrow(pool) rows. bilinear extraction yields NA climate
            # for plots off the raster's edge; drop those rows from BOTH pool
            # and response together so the row counts stay locked. (The old
            # two-fold path tolerated mismatches because build_analog_index()
            # never row-checked y; analog_cv() checks immediately.)
            complete <- stats::complete.cases(d) & !is.na(values)
            d <- d[complete, , drop = FALSE]
            values <- values[complete]
            stopifnot(nrow(d) == length(values))

            # The CV pool in AIM coordinates, so held-out predictions use the
            # same distances as the maps. `d` stays lon/lat for scopes and plots.
            d_aim <- data.frame(aim_xy(d), d[, clim_vars, drop = FALSE])

            # Scope membership, aligned with the filtered rows -- and therefore
            # with the rows analog_cv() returns. Computed once per variable set
            # rather than per bandwidth cell, since it depends on neither.
            #
            # Only on the summarizing path. eval() is called with summarize=TRUE
            # exactly when vary_bandwidths is TRUE, and the scatter path resolves
            # its own scope in the renderer. Computing this unconditionally reads
            # site$x and the polygon, which would make eval_scatter_data() depend
            # on them -- and a single map click would then re-run the whole LOO
            # cross-validation to produce an identical plot.
            scopes <- if (vary_bandwidths) eval_scope_membership(d[, c("x", "y")]) else NULL

            # Single-configuration evaluation: leave-one-out CV via the
            # analogs package's native analog_cv(), with performance metrics
            # from cv_performance(). analog_cv() predicts each plot from a
            # neighborhood excluding its own pool row (exclude_self, handled
            # internally); cv_performance() then computes held-out metrics
            # from the obs/residual columns.
            #
            # Note on cost: analog_cv() rebuilds the analog index internally
            # on every call (it rejects pre-built indices by design), so the
            # bandwidth grid below incurs one index build per cell. A direct
            # analog_impact() call with a hoisted index cannot avoid this --
            # exclude_self requires raw pool data, which forces an internal
            # rebuild regardless. analog_cv() is retained here for clarity and
            # correctness; the per-cell rebuild is an accepted cost.
            #
            # cv_performance() auto-classifies the response as binary
            # (obs all in {0,1}, both classes present) or continuous and
            # returns the corresponding metric set in long format
            # (columns: variable, type, metric, value):
            #   continuous -> n, rmse, mae, bias, r2
            #   binary     -> n, auc, tss, tss_threshold, brier
            # `r2` is out-of-sample R-squared (1 - SS_res/SS_tot), which can
            # be negative -- unlike the old cor(pred,obs)^2. There is no AIC
            # equivalent: AIC is an in-sample criterion and cv_performance()
            # is a held-out CV tool, so `delta AIC` was dropped from the UI.
            eval <- function(theta_clim, theta_geog, summarize = TRUE){
                  k <- aim_kernels(theta_clim, theta_geog)
                  cv <- suppressWarnings(analogs::analog_cv(
                        fun        = analogs::analog_impact,
                        pool       = d_aim,
                        y          = values,
                        cv_method  = "loo",
                        env        = k$env,
                        geog       = k$geog,
                        stat       = "weighted_mean",
                        coord_type = "projected"))

                  if(!summarize){
                        # Scatterplot path: return the cv object itself plus the
                        # plot coordinates. Scoping is deferred to the renderer
                        # for the same reason it is here -- membership depends on
                        # the polygon and the focal site, and resolving it inside
                        # this function would make every map click re-run a full
                        # LOO cross-validation. Returning cv (rather than just
                        # obs/pred) lets the renderer call cv_performance() on a
                        # row subset, exactly as the heatmap path does.
                        return(list(cv = cv, xy = d[, c("x", "y"), drop = FALSE]))
                  }

                  # Heatmap path: native CV performance metrics, computed once
                  # per scope from the SAME cv object. Scoping is a row subset,
                  # so every scope's metrics come out of one run for the cost of
                  # an extra aggregation -- which is why the scope selector can
                  # be a display control rather than something that recomputes.
                  #
                  # Rename `metric` -> `stat` so downstream consumers
                  # (eval_heatmap, optimal_params) keep their existing column
                  # vocabulary, and tag with the bandwidth coordinates.
                  purrr::imap(scopes, function(keep, scope_name){
                        if (sum(keep) < eval_min_n) return(NULL)
                        # A small subset can hold a single class, leaving AUC and
                        # TSS undefined. Drop that scope's cell rather than
                        # failing the whole run.
                        p <- tryCatch(analogs::cv_performance(cv[keep, , drop = FALSE]),
                                      error = function(e) NULL)
                        if (is.null(p)) return(NULL)
                        p %>%
                              dplyr::rename(stat = metric) %>%
                              dplyr::select(stat, value) %>%
                              dplyr::mutate(scope = scope_name,
                                            theta_clim = theta_clim,
                                            theta_geog = theta_geog)
                  }) %>%
                        purrr::list_rbind()
            }

            if(vary_bandwidths){
                  mult <- c(.25, .5, 1, 2, 4)
                  clim_mult <- if("clim" %in% input$eval_params) mult else 1
                  geog_mult <- if("geog" %in% input$eval_params) mult else 1
                  params <- tidyr::expand_grid(theta_clim = input$theta_clim * clim_mult,
                                               theta_geog = input$theta_geog * geog_mult)

                  e <- purrr::pmap(params, eval) %>%
                        purrr::list_rbind() %>%
                        mutate(clim_vars = paste(clim_vars, collapse = "_"))
            }else{
                  e <- eval(input$theta_clim, input$theta_geog,
                            summarize = FALSE)
            }

            return(e)
      }

      observeEvent(input$run_eval, {
            # Refuse rather than sweeping bandwidths against a response the
            # module cannot score.
            if (!is.null(eval_unsupported())) {
                  showNotification(eval_unsupported(), type = "warning", duration = 10)
                  return()
            }

            if("vars" %in% input$eval_params){
                  cv <- input$clim_vars
                  n <- min(input$eval_n_vars[1], length(cv)):min(input$eval_n_vars[2], length(cv))
                  v <- lapply(n, function(x) apply(combn(cv, x), 2, list))
                  v <- unlist(unlist(v, recursive = F), recursive = F)
            }else{
                  v <- list(input$clim_vars)
            }

            e <- v %>% lapply(eval_var_set) %>% bind_rows()

            model_eval$comparison <- e
            # Snapshot the settings these results describe, plus the dimensions
            # this run deliberately varied. Both are needed: the stamp detects
            # change, the swept list keeps the run from indicting itself for the
            # very dimensions it was exploring.
            model_eval$stamp <- eval_run_stamp()
            model_eval$swept <- c(
                  if ("clim" %in% input$eval_params) "theta_clim",
                  if ("geog" %in% input$eval_params) "theta_geog",
                  if ("vars" %in% input$eval_params) "clim_vars"
            )
      })

      ## Staleness marker. Per-scope metrics are aggregated at run time, so they
      ## cannot be recomputed for new settings without re-running the CV.
      ## Storing per-plot residuals instead would make rescoping free, but the
      ## eval grid can reach hundreds of configurations across ~88k pool plots --
      ## tens of millions of rows. Pre-aggregation is the right call; explicit
      ## invalidation is its cost.
      output$eval_stale <- renderUI({
            stale_banner(eval_stale_fields(), eval_stamp_labels,
                         action = "Click 'Run evaluations' to update them.")
      })

      observeEvent(input$optimal_params, {
            req(model_eval$comparison)

            # Refuse rather than silently applying bandwidths derived from
            # settings the user has since moved away from. The warning banner is
            # visible, but a button that ignores its own panel's warning is worse
            # than one that explains itself.
            stale <- eval_stale_fields()
            if (length(stale)) {
                  showNotification(
                        sprintf("These results predate a change to the %s. Re-run evaluations before applying optimal parameters.",
                                stale_phrase(stale, eval_stamp_labels)),
                        type = "warning", duration = 8)
                  return()
            }

            metric <- switch(input$eval_stat,
                             "RMSE" = "rmse",
                             "R-squared" = "r2",
                             "AUC" = "auc",
                             "TSS" = "tss",
                             "Brier" = "brier")

            # Metrics where larger is better; the rest are minimized.
            higher_better <- c("r2", "auc", "tss")

            # Optimize against the SELECTED scope: the best bandwidth for a
            # particular unit is not generally the best bandwidth domain-wide,
            # and the whole point of the scope selector is to let the user ask
            # for the first one.
            opt <- model_eval$comparison %>% ungroup() %>%
                  filter(scope == input$eval_scope, stat == metric) %>%
                  mutate(value = ifelse(stat %in% higher_better, value, -value)) %>%
                  filter(value == max(value))

            shiny::validate(need(nrow(opt) > 0,
                                 "No results for this area. Re-run the evaluations, or choose a larger area."))

            updateSliderInput(session, "theta_geog", value = opt$theta_geog)
            updateSliderInput(session, "theta_clim", value = opt$theta_clim)
            updateSelectizeInput(session, "clim_vars", selected = unlist(strsplit(opt$clim_vars, "_")))
      })

      output$eval_heatmap <- renderPlot({
            shiny::validate(need(is.null(eval_unsupported()), eval_unsupported()))
            req(model_eval$comparison)

            metric <- switch(input$eval_stat,
                             "RMSE" = "rmse",
                             "R-squared" = "r2",
                             "AUC" = "auc",
                             "TSS" = "tss",
                             "Brier" = "brier")

            # Metrics where larger is better; the rest are minimized.
            opt <- if(metric %in% c("r2", "auc", "tss")) 1 else -1

            scope_label <- names(eval_scope_names)[eval_scope_names == input$eval_scope]

            e <- model_eval$comparison %>% ungroup() %>%
                  filter(scope == input$eval_scope)

            # No rows means one of three things: the scope's input isn't set,
            # fewer than eval_min_n plots fall inside it, or the last run
            # predates the selection. Name the first case specifically, since it
            # is the only one the user can act on immediately.
            scope_set <- switch(input$eval_scope,
                                "polygon"  = ,
                                "buffered" = isTruthy(polygon_upload_active()),
                                "site"     = !is.null(site$x) && !is.null(site$y),
                                TRUE)
            shiny::validate(need(
                  nrow(e) > 0,
                  if (!scope_set) eval_scope_missing_msg(input$eval_scope)
                  else sprintf("No results for '%s'. Either the last run was done before you chose this area, or fewer than %d plots fall inside it.",
                               scope_label, eval_min_n)))

            # cv_performance() reports n alongside the requested metric. It can
            # vary across bandwidth cells, since a plot with no analog in range
            # yields no held-out prediction, so report the range rather than
            # implying a single number.
            n_vals <- e %>% filter(stat == "n") %>% pull(value)
            n_text <- if (!length(n_vals)) "n unknown"
            else if (min(n_vals) == max(n_vals)) sprintf("n = %d", as.integer(min(n_vals)))
            else sprintf("n = %d\u2013%d", as.integer(min(n_vals)), as.integer(max(n_vals)))

            e <- e %>%
                  filter(stat == metric) %>%
                  mutate(label = signif(value, 2),
                         optimal = value * opt == max(value * opt),
                         text_color = case_when(optimal ~ "red",
                                                scales::rescale(value) < .5 ~ "white",
                                                TRUE ~ "black"),
                         clim_vars = gsub("_", " + ", clim_vars),
                         n_vars = stringr::str_count(clim_vars, "\\+") + 1)

            shiny::validate(need(nrow(e) > 0,
                                 sprintf("'%s' cannot be calculated for this variable within '%s'.",
                                         input$eval_stat, scope_label)))

            subtitle <- paste0(scope_label, " \u2014 ", n_text,
                               if (length(n_vals) && min(n_vals) < 50)
                                     "  (small samples; treat these numbers as indicative)" else "")

            p <- ggplot(e, aes(theta_geog, theta_clim, fill = value, label = label)) +
                  facet_wrap(~clim_vars) +
                  geom_tile() +
                  geom_text(aes(color = I(text_color), fontface = "bold")) +
                  scale_x_log10(breaks = unique(e$theta_geog)) +
                  scale_y_log10(breaks = unique(e$theta_clim)) +
                  scale_fill_viridis_c() +
                  guides(fill = guide_colorbar(barwidth = 12)) +
                  theme(legend.position = "bottom") +
                  strip_style +
                  labs(y = "climate bandwidth (z)",
                       x = "geographic bandwidth (km)",
                       subtitle = subtitle,
                       fill = paste0(input$eval_stat, " (", eval_response_label(), ")"))

            # `%>%` binds tighter than `+`, so this cannot be chained onto the
            # ggplot expression above -- it would pipe only the labs() object.
            gg_stale(p, length(eval_stale_fields()) > 0)
      })


      eval_scatter_data <- reactive({
            eval_var_set(input$clim_vars, vary_bandwidths = FALSE)
      })

      # Scope resolution for the single-configuration section. Shared by the
      # scatter and its stat readout so the two can never disagree about which
      # plots they are describing.
      eval_scatter_scope <- reactive({
            e <- eval_scatter_data()
            req(e)
            scope <- if (isTruthy(input$eval_scope)) input$eval_scope else "domain"
            keep <- eval_scope_membership(e$xy, which = scope)[[scope]]
            # NULL means the scope's defining input isn't set. Falling back to
            # every plot would label domain-wide results with the selected
            # scope's name, which is worse than showing nothing.
            shiny::validate(need(!is.null(keep), eval_scope_missing_msg(scope)))
            list(scope = scope,
                 label = names(eval_scope_names)[eval_scope_names == scope],
                 keep  = keep)
      })

      output$eval_scatter <- renderPlot({
            shiny::validate(need(is.null(eval_unsupported()), eval_unsupported()))
            e <- eval_scatter_data()
            s <- eval_scatter_scope()
            req(e, s)

            d <- tibble::tibble(obs      = e$cv$obs,
                                pred     = e$cv$obs - e$cv$residual,
                                in_scope = s$keep) %>%
                  filter(!is.na(obs), !is.na(pred))

            shiny::validate(need(
                  sum(d$in_scope) >= eval_min_n,
                  sprintf("Only %d plot(s) fall within '%s' \u2014 too few to judge accuracy from.",
                          sum(d$in_scope), s$label)))

            # One calibration plot for every response type: bin the predictions,
            # and for each bin show what was actually observed.
            #
            # Predicted on the x axis, not observed, because that is the direction
            # the tool is used in. "Given a plot with 7 species, what does the
            # model say" is the inverse of anyone's question; "the map reads 3.2
            # here, what is really there" is the question, and this answers it
            # directly. It also reframes shrinkage as a usable correction rather
            # than an alarm: predictions that plateau at 3.2 against observations
            # running to 10 read as "when it says 3.2, expect about 5".
            #
            # Equal-COUNT bins, so every point carries the same standard error and
            # none is more trustworthy than another. That is also why the points
            # are not sized: sizing by bin count stretched a range of 1,142 to
            # 1,143 plots across the full radius scale, turning rounding noise
            # into apparently dramatic differences. The rug shows where the
            # predictions actually fall, which is the real density question.
            discrete_obs <- all(abs(d$obs - round(d$obs)) < 1e-8, na.rm = TRUE)
            binary_obs   <- all(d$obs %in% c(0, 1))

            # Equal-count bins via explicit quantile breaks, de-duplicated.
            # cut_number() fails outright when quantiles tie, and many responses
            # have a large spike of identical predictions -- proportion basal area
            # is zero wherever the species is absent. Falling back to equal-width
            # there produced one bin of 77,000 plots and others of 8, which is
            # what made the tail of the curve fly around.
            #
            # Merging tied breaks instead gives one wide bin covering the spike
            # and equal-count bins across the rest: counts still vary, but by a
            # factor rather than by four orders of magnitude.
            calib <- function(x) {
                  nb <- max(3, min(20, floor(nrow(x) / 200)))
                  qs <- unique(stats::quantile(x$pred, seq(0, 1, length.out = nb + 1),
                                               na.rm = TRUE))
                  x$bin <- if (length(qs) >= 4) cut(x$pred, qs, include.lowest = TRUE)
                  else cut(x$pred, min(nb, 5))
                  x %>%
                        group_by(bin) %>%
                        summarize(p  = mean(pred),
                                  o  = mean(obs),
                                  lo = stats::quantile(obs, 0.25),
                                  hi = stats::quantile(obs, 0.75),
                                  n  = dplyr::n(),
                                  .groups = "drop") %>%
                        # A mean of a handful of plots is noise dressed as a
                        # calibration point.
                        filter(n >= 30)
            }

            f_in  <- d[d$in_scope, , drop = FALSE]
            f_out <- d[!d$in_scope, , drop = FALSE]
            b <- calib(f_in)

            shiny::validate(need(nrow(b) >= 3,
                                 "Not enough spread in the predictions to assess calibration."))

            # Equal-count binning can fail when many predictions tie -- a
            # presence/absence response has thousands at effectively zero -- in
            # which case the fallback is equal-width and the counts vary widely.
            # Report which happened rather than claiming a per-point sample size
            # that is not true.
            n_per_bin <- if (diff(range(b$n)) <= 0.02 * stats::median(b$n))
                  sprintf("~%s plots", format(round(stats::median(b$n)), big.mark = ","))
            else sprintf("%s-%s plots",
                         format(min(b$n), big.mark = ","),
                         format(max(b$n), big.mark = ","))

            p <- ggplot(b, aes(p, o)) +
                  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                              color = "dodgerblue")

            # Raw points for continuous responses only. Their observed values are
            # spread out, so the cloud shows the joint distribution the bins
            # summarize; a discrete response would just stack them into stripes.
            if (!discrete_obs) p <- p +
                  geom_point(data = f_in, aes(pred, obs), inherit.aes = FALSE,
                             color = "gray55",
                             # Floor at 0.12, not 0.03: the low floor was carried
                             # over from a version where points stacked on top of
                             # each other and reinforced. Here they are spread
                             # across a cloud, so each one has to be visible on
                             # its own.
                             alpha = max(0.12, min(0.5, 8000 / max(nrow(f_in), 1))),
                             size = 0.8)

            p <- p +
                  geom_rug(data = f_in, aes(x = pred), inherit.aes = FALSE,
                           alpha = max(0.01, min(0.3, 2000 / nrow(f_in))),
                           sides = "b", color = "gray30")

            # When a scope narrows the evaluation, the same curve for everything
            # OUTSIDE it goes underneath in gray. That answers the question the
            # scope selector exists to raise -- is the model as well calibrated
            # here as it is elsewhere -- which a single curve cannot.
            if (nrow(f_out) >= 200) {
                  b_out <- calib(f_out)
                  if (nrow(b_out) >= 3) p <- p +
                              geom_line(data = b_out, color = "gray65", linewidth = 0.7)
            }

            # Interquartile spread of what was observed in each bin, as a band
            # rather than bars on the points. The point is the MEAN, which for a
            # skewed response legitimately falls outside the quartiles -- drawn as
            # bars attached to each point that reads as an error bar that has come
            # adrift. A band reads as what it is: the spread of outcomes behind
            # the curve, a separate quantity from the calibration estimate.
            #
            # Skipped for presence/absence, where the quartiles of a 0/1 variable
            # carry nothing the mean does not.
            if (!binary_obs) p <- p +
                  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "red", alpha = 0.12)

            p <- p +
                  geom_line(color = "red", linewidth = 0.9) +
                  geom_point(color = "red", size = 2.2)

            # sqrt axes tame the skew of abundance variables; on counts of 1-12 or
            # a probability they buy nothing.
            if (!discrete_obs) p <- p + scale_x_sqrt() + scale_y_sqrt()

            p +
                  coord_fixed() +
                  theme_bw() +
                  labs(x = paste0("predicted ", eval_response_label()),
                       y = if (binary_obs) "observed frequency"
                       else paste0("observed ", eval_response_label()),
                       title = paste0(
                             "Model calibration curve\n(", s$label, "; ",
                             n_per_bin, " per point)"))
      })

      # Held-out accuracy for the current settings, over the selected scope.
      # Reports the whole metric set cv_performance() returns rather than only
      # input$eval_stat: for a single configuration there is no cost to showing
      # all of them, and they answer different questions (bias and R-squared
      # disagree in useful ways). Deliberately NOT faded by gg_stale -- this
      # section is computed from current settings, so unlike the heatmap it is
      # never stale, and that contrast is the point.
      output$eval_scatter_stats <- renderUI({
            req(is.null(eval_unsupported()))
            e <- eval_scatter_data()
            s <- eval_scatter_scope()
            req(e, s)
            if (sum(s$keep) < eval_min_n) return(NULL)   # scatter says why

            p <- tryCatch(analogs::cv_performance(e$cv[s$keep, , drop = FALSE]),
                          error = function(err) NULL)
            if (is.null(p) || !nrow(p)) return(NULL)

            fmt <- function(metric, value) {
                  if (identical(metric, "n")) format(as.integer(value), big.mark = ",")
                  else signif(value, 3)
            }
            items <- p %>%
                  filter(metric %in% names(eval_metric_labels)) %>%
                  arrange(match(metric, names(eval_metric_labels)))

            div(class = "mb-2", style = "font-size: 0.9em;",
                tags$span(style = "color: #666;", paste0(s$label, ": ")),
                lapply(seq_len(nrow(items)), function(i) {
                      tags$span(style = "margin-right: 1.1em;",
                                tags$span(style = "color: #666;",
                                          eval_metric_labels[[items$metric[i]]]),
                                " ",
                                tags$strong(fmt(items$metric[i], items$value[i])))
                }))
      })


      # CLIMATE RESPONSE PLOT --------------------------------
      output$climate_response <- renderPlot({
            req(display_layer(), clim())

            req(isTruthy(input$cr_v1), input$cr_v1 %in% display_vars())
            v1 <- input$cr_v1

            cr_scope <- if (isTruthy(input$cr_scope)) input$cr_scope else "domain"
            scope_set <- switch(cr_scope,
                                "polygon"  = ,
                                "buffered" = isTruthy(polygon_upload_active()),
                                "site"     = !is.null(site$x) && !is.null(site$y),
                                TRUE)
            shiny::validate(need(scope_set, eval_scope_missing_msg(cr_scope)))
            scope_label <- names(eval_scope_names)[eval_scope_names == cr_scope]

            # Climate for the x-axis. Historical by default; 'future' lets you
            # pair a future prediction with the future climate that produced it,
            # rather than reading it against current climate.
            use_fut <- isTRUE(input$cr_clim == "future")

            if(input$cr_data == "modeled grid cells"){
                  cr <- if (use_fut) cr_grid()$fut else cr_grid()$hst
                  d <- rast(c(y = display_layer(), clim1 = cr)) %>%
                        cr_scope_mask() %>%
                        as.data.frame()
            }

            if(input$cr_data == "FIA plots"){
                  # display_values(), not response_values(): in community mode the
                  # per-plot counterpart of the community measure -- observed
                  # richness, say -- so this path describes the same quantity the
                  # grid-cell path does. response_values() would silently plot the
                  # focal species the entity picker still holds while hidden.
                  #
                  # It also covers the pooled-ratio variables, whose pool()$values
                  # carries occ/num/denom rather than a `response` column.
                  dv <- display_values()
                  shiny::validate(need(
                        !is.null(dv),
                        "This community measure has no per-plot counterpart. Switch the data to 'modeled grid cells'."))
                  cv <- cr_plot_values()
                  req(length(cv) == length(dv))
                  d <- data.frame(y = dv, clim1 = cv)

                  # Same membership test the evaluation module uses, so "focal
                  # polygon" selects the same plots in both.
                  keep <- eval_scope_membership(plot_occ()[, c("x", "y")],
                                                which = cr_scope)[[cr_scope]]
                  req(!is.null(keep), length(keep) == nrow(d))
                  d <- d[keep, , drop = FALSE]
            }

            # Climate rasters are stored standardized, so back-convert to native
            # units (z * sd + mean) before plotting -- this puts the x-axis in
            # mm / degrees C. Falls back to leaving a variable standardized if it
            # has no scale entry.
            back_transform <- function(x, v) {
                  s <- clim_scales[[v]]
                  if (is.null(s) || !is.finite(s$sd)) return(x)
                  x * s$sd + s$mean
            }
            d$clim1 <- back_transform(d$clim1, v1)

            # A categorical layer inverts the plot rather than adapting it. There
            # is no curve to fit through entity identities, but the question
            # transposes cleanly: instead of "how does the value change along this
            # gradient", it becomes "where along this gradient does each entity
            # win". Horizontal violins answer that -- one distribution of the
            # climate variable per entity, read down the page.
            if (categorical_layer()) {
                  cc <- comm_categories()
                  d$label <- cc$labels[match(d$y, cc$codes)]
                  d <- d[!is.na(d$label) & is.finite(d$clim1), , drop = FALSE]

                  shiny::validate(need(nrow(d) > 0,
                                       "No values to plot for this layer over the current area."))

                  # Ordered by median position along the climate gradient, so the
                  # violins step across the panel and the axis itself reads as the
                  # gradient. Frequency order left the entities scattered, which
                  # made the sequence -- the point of the plot -- something the
                  # reader had to reconstruct.
                  #
                  # Median rather than mean: these distributions are often skewed
                  # or bimodal for a species spanning two climate zones, and the
                  # mean can land between the modes, in a range the entity does
                  # not actually occupy.
                  #
                  # "other" stays pinned to the bottom regardless. It is a pooled
                  # remainder with no meaningful position on the gradient, and
                  # letting it sort among the entities would imply it had one.
                  med <- d %>%
                        filter(label != "other") %>%
                        group_by(label) %>%
                        summarize(m = stats::median(clim1, na.rm = TRUE), .groups = "drop") %>%
                        arrange(m)
                  lev <- c(if ("other" %in% d$label) "other", as.character(med$label))
                  d$label <- factor(d$label, levels = lev)

                  p_main <- ggplot(d, aes(clim1, label, fill = label)) +
                        geom_violin(scale = "width", color = NA) +
                        # Same entity colors as the map and the composition panel,
                        # so an entity is recognizable across all of them.
                        scale_fill_manual(values = stats::setNames(cc$colors, cc$labels),
                                          guide = "none") +
                        theme_bw() +
                        labs(x = paste0(input$cr_clim, " ", v1, " (", unit_lab(v1), ")"),
                             y = NULL,
                             subtitle = scope_label,
                             caption = if (!cr_var_in_model())
                                   paste0(v1, " is not one of the variables the model was fitted on; this relationship comes from its correlation with those that were."))

                  p_dens <- ggplot(d, aes(clim1)) +
                        geom_density(fill = "#2b6a4d", color = NA, alpha = .35) +
                        theme_void()

                  return(p_dens / p_main + plot_layout(heights = c(1, 8)))
            }

            # Smooth on a link scale so fitted curves respect each variable's
            # natural bounds: logit for [0, 1] variables, log for non-negative
            # variables, identity for differences. Uses a GLM with a natural-spline
            # basis; quasi-families accept the continuous grid-cell values without
            # integer/proportion warnings. The FIA-plots dataset is the raw
            # historical per-plot variable, so its family depends on the variable
            # alone; only the modeled-grid dataset honors the raster stat and time.
            # In community mode the plotted quantity is the community measure,
            # not the per-entity variable: richness and diversity are unbounded
            # counts however they were built, while evenness, dominance, turnover
            # and the two change fractions live on [0, 1].
            bounded01 <- if (multi_mode())
                  isTruthy(input$comm_measure) &&
                  input$comm_measure %in% c("evenness", "dominance", "turnover",
                                            "colonization demand", "loss fraction")
            else grepl("presence|fraction", input$fia_var) ||
                  input$fia_var == "proportion basal area"
            fam <- if (input$cr_data == "FIA plots") {
                  if (bounded01) quasibinomial() else quasipoisson()
            } else if (input$stat == "connectivity" &&
                       isTRUE(input$cs_stat == "current direction")) {
                  gaussian()                      # circular (0-360 deg): a log link is meaningless
            } else if (input$stat != "FIA variable") {
                  quasipoisson()                  # ESS / similarity / connectivity magnitudes: non-negative
            } else if (input$time == "delta") {
                  gaussian()                      # future - historical difference (can be negative)
            } else if (bounded01) {
                  quasibinomial()                 # bounded [0, 1]
            } else {
                  quasipoisson()                  # counts / basal area: non-negative
            }

            # Base the y-axis on the fitted lines, not the SE ribbons -- the
            # ribbons can blow up where a group's smooth extrapolates into a
            # sparse tail. Refit the same per-group smooth, take the range of the
            # fitted values, and clip the display there with coord_cartesian (the
            # geom_smooth fit itself still uses all the data; only the view is
            # bounded).
            dd <- d[is.finite(d$y) & is.finite(d$clim1), ]
            shiny::validate(need(
                  nrow(dd) > 0,
                  "No values to plot for this layer over the current area."))
            fits <- tryCatch(fitted(glm(y ~ splines::ns(clim1, 4), family = fam, data = dd)),
                             error = function(e) numeric(0))
            # range() of an all-NA vector returns c(Inf, -Inf), and handing that
            # to coord_cartesian() blanks the panel with no error -- the plot
            # renders, the axes render, and nothing is drawn.
            ylim <- if (length(fits) && any(is.finite(fits)))
                  range(fits[is.finite(fits)]) else NULL
            if (!is.null(ylim) && (!all(is.finite(ylim)) || diff(ylim) == 0)) ylim <- NULL

            p_main <- ggplot(d, aes(clim1, y)) +
                  geom_smooth(method = "glm",
                              formula = y ~ splines::ns(x, 4),
                              method.args = list(family = fam),
                              color = "#2b6a4d", fill = "#2b6a4d") +
                  coord_cartesian(ylim = ylim) +
                  # Label the climate axes with their timeframe (so historical vs
                  # future is never ambiguous) and their native units.
                  labs(x = paste0(input$cr_clim, " ", v1, " (", unit_lab(v1), ")"),
                       subtitle = scope_label,
                       # FIA-plots data is raw historical observation, so label it
                       # by the variable rather than the (stat/time-aware) raster title.
                       y = if (input$cr_data == "FIA plots")
                             paste("observed", if (multi_mode()) comm_label() else input$fia_var)
                       else raster_title_text(),
                       caption = if (!cr_var_in_model() && input$cr_data != "FIA plots")
                             paste0(v1, " is not one of the variables the model was fitted on; this relationship comes from its correlation with those that were."))

            # Marginal density along the x variable: shows where the data
            # actually is, which is what explains any wandering at the ends of the
            # smooth. It carries its own y-scale, so it does not interact with the
            # fitted-line ylim on the main panel.
            p_dens <- ggplot(d, aes(clim1)) +
                  geom_density(fill = "#2b6a4d", color = NA, alpha = .35) +
                  theme_void()

            # patchwork aligns the two panels on a shared x-axis; both are built
            # from the same clim1.
            p_dens / p_main + plot_layout(heights = c(1, 8))
      })


      # DOWNLOADS --------------------------------------
      output$download_raster <- downloadHandler(
            filename = function() {
                  paste(input$download_filename, ".tif", sep="")
            },

            content = function(file) {
                  terra::writeRaster(display_layer(), file)
            }
      )
      output$download_points <- downloadHandler(
            filename = function() {
                  paste(input$download_filename, "_points.csv", sep="")
            },

            content = function(file) {
                  # isolate() because this runs outside a render context, where a
                  # req() reaching up through the reactive chain surfaces as a 500
                  # (which the browser saves as download_points.html) rather than
                  # as a quiet no-op. The realistic path here is "focal site
                  # analogs" selected with no site clicked yet; a header-only file
                  # is the truthful answer to that, and always a valid CSV.
                  pd <- tryCatch(isolate(export_points()),
                                 error = function(e) export_points_empty())

                  write_csv(pd, file)
            }
      )
      output$download_metadata <- downloadHandler(
            filename = function() {
                  paste(input$download_filename, ".txt", sep="")
            },

            content = function(file) {
                  params <- names(input)
                  out <- c()
                  for(p in params){
                        out <- c(out, paste(p, input[[p]], sep = ": "))
                  }
                  if(input$stat != "connectivity"){
                        out <- out[!grepl("cs_", params)]
                  }
                  writeLines(out, file)
            }
      )


}

shinyApp(ui, server)
