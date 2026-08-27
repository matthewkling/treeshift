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
unit_lab <- function(v) {
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
eval_scope_names <- c("Modeling domain"         = "domain",
                      "Study area"              = "polygon",
                      "Buffered study area"     = "buffered",
                      "Focal site neighborhood" = "site")

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
                       fia_var    = "FIA variable",
                       bbox       = "modeling domain",
                       clim_vars  = "climate variables",
                       theta_clim = "climate bandwidth",
                       theta_geog = "geographic bandwidth",
                       polygon    = "study area polygon",
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
      div(class = "alert alert-warning py-2 px-3 mb-2",
          style = "font-size: 0.85em;",
          sprintf("Computed before the %s changed. %s",
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
            theme(plot.caption = element_text(color = "#b8860b", face = "bold",
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
                          and climatic similarity to FIA plots, optionally accounting for pedicted future climate change.
                          Use it to assess potential shifts in species distributions and other forest characteristics,
                          model connectivity between current and potential future distributions,
                          and explore properties of plots that are analogs to your locations of interest."),
                        br(),
                        helpText("Here on the lefthand sidebar, you'll find core controls for setting model parameters and
                                 choosing display variables. The righthand sidebar houses modules with their own
                                 dedicated controls and plots. Hover your cursor over controls for explanatory pop-ups."),
                        br(),
                        helpText("NOTE: This tool is a prototype and is still under development.
                          Please contact matthew.kling@umontana.edu with questions, suggestions, or bug reports.")
                  ),

                  accordion_panel(
                        "FIA variable",
                        selectInput("track", "Focal entity",
                                    choices = c("Species", "Forest type", "Forest group", "All trees"),
                                    selected = "Species") %>%
                              tooltip(HTML("Choose the focal entity for the model.<br><br>
                                          <dl style='margin:0;'>
                                              <dt>Species</dt> <dd>Model presence or abundance of a single focal tree species.</dd>
                                              <dt>Forest type</dt> <dd>Model probability that a site is classified as a focal FIA forest type (a fine-grained classification).</dd>
                                              <dt>Forest group</dt> <dd>Model probability that a site is classified within a focal FIA forest type group (a coarser classification).</dd>
                                              <dt>All trees</dt> <dd>Model community-wide forest properties (any-species presence or total basal area).</dd>
                                          </dl>")),
                        conditionalPanel(
                              condition = "input.track == 'Species'",
                              selectizeInput("species", "Species", choices = spp, selected = "Pinus edulis") %>%
                                    tooltip("Select a focal species by scrolling or typing")
                        ),
                        conditionalPanel(
                              condition = "input.track == 'Forest type'",
                              selectizeInput("forest_type", "Forest type",
                                             choices = forest_type_choices,
                                             selected = forest_type_choices[[1]]) %>%
                                    tooltip("Select a focal FIA forest type by scrolling or typing. Labels show the FORTYPCD followed by the type name.")
                        ),
                        conditionalPanel(
                              condition = "input.track == 'Forest group'",
                              selectizeInput("forest_group", "Forest group",
                                             choices = forest_group_choices,
                                             selected = forest_group_choices[[1]]) %>%
                                    tooltip("Select a focal FIA forest type group by scrolling or typing.")
                        ),
                        selectInput("metric", "Variable",
                                    c("occurrence", "abundance", "basal area",
                                      "proportion basal area", "regeneration fraction")) %>%
                              tooltip(HTML("Choose the quantity to model.<br><br>
                                          <dl style='margin:0;'>
                                              <dt>Occurrence</dt> <dd>Probability the focal entity is present within a single FIA plot footprint (~0.07 ha) randomly located within the grid cell.</dd>
                                              <dt>Abundance</dt> <dd>Population density of the selected life stage, in stems per acre.</dd>
                                              <dt>Basal area</dt> <dd>Live tree and sapling basal area, in ft&sup2; per acre.</dd>
                                              <dt>Proportion basal area</dt> <dd>The focal species' share of total stand basal area, multiplied by the probability that stand basal area is nonzero.</dd>
                                              <dt>Regeneration fraction</dt> <dd>Probability that the selected regeneration stage is present relative to adults, conditional on the species occuring on the plot. Calcualted as p(regen) / (p(regen) + p(adult)). Adults only = 0; regen only = 1; adults and regen together = intermediate value (not precisely 0.5 due to differing sample areas).</dd>
                                          </dl>")),
                        ## Life stage: shown only for metrics/tracks that have one
                        ## (occurrence / abundance / regeneration fraction on the
                        ## Species and All trees tracks). Basal-area metrics and the
                        ## forest type/group tracks have no stage.
                        conditionalPanel(
                              condition = "(input.track == 'Species' || input.track == 'All trees') && (input.metric == 'occurrence' || input.metric == 'abundance' || input.metric == 'regeneration fraction')",
                              selectInput("stage", "Life stage",
                                          c("seedling", "sapling", "adult", "any"),
                                          selected = "any") %>%
                                    tooltip(HTML("Choose the life stage (size class).<br><br>
                                                <dl style='margin:0;'>
                                                    <dt>Seedling</dt> <dd>&lt; 1 in. DBH/DRC.</dd>
                                                    <dt>Sapling</dt> <dd>1.0&ndash;4.9 in. DBH/DRC.</dd>
                                                    <dt>Adult</dt> <dd>&ge; 5.0 in. DBH/DRC.</dd>
                                                    <dt>Any</dt> <dd>All life stages combined.</dd>
                                                </dl>"))
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
                                               "seedling fraction", "sapling fraction")))
                  ),

                  ## Three nested spatial extents, one shared resolution.
                  ##
                  ## The modeling domain, the study area polygon, and the focal
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
                        helpText("Where the model runs, and the areas you compare against it:
                                 the modeling domain, an optional study area polygon, and an
                                 optional focal site. All three are evaluated at the grain set
                                 below."), hr(),

                        sliderInput("clim_agg", "Raster grain size", 1, 10, value = 10, step = 1, post = " km") %>%
                              tooltip("Resolution of the modeling grid, and of every summary drawn from it. Larger values are faster but less precise."),

                        hr(),

                        tags$strong("Modeling domain"),
                        helpText("The extent the model is computed over."),
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

                        tags$strong("Study area polygon"),
                        helpText("An area within the domain to summarize over. Optional."),
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
                        layout_column_wrap(
                              width = 1/2,
                              actionButton("highlight_site", "Highlight on map") %>%
                                    tooltip("Switch the vector layer to this site's analog plots and zoom the map to its geographic cutoff."),
                              actionButton("site_clear", "Clear") %>%
                                    tooltip("Remove the focal site from the map and plots.")
                        )
                  ),

                  accordion_panel(
                        "Climate data",
                        helpText("Select climate variables and future scenarios."), hr(),
                        selectInput("clim_vars", "Climate variables", clim_vars(), selected = c("cwd", "aet"), multiple = TRUE) %>%
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
                        helpText("Analog FIA plots are weighted using Gaussian kernels based on climatic and geographic distances from a site.
                                 Bandwidths specify the standard deviations of these kernels. Plots farther than 3 stdevs are considered non-analog and excluded."), hr(),
                        sliderInput("theta_geog", "Geographic bandwidth (km)", 5, 100, value = 30),
                        sliderInput("theta_clim", "Climate bandwidth (z)", .05, .5, value = .1) %>%
                              tooltip("Climate bandwidth units are in standard deviations relative to climate variation across all FIA plots in the contiguous US")
                  ),

                  accordion_panel(
                        "Map layers",
                        helpText("Choose which model results are drawn on the map."), hr(),
                        selectInput("stat", "Raster variable", stat_choices_base) %>%
                              tooltip(HTML("<dl style='margin:0;'>
                                              <dt>FIA variable</dt> <dd>Interpolated values for the plot variable selected above</dd>
                                              <dt>Effective sample size</dt> <dd>Effective number of FIA plots within climatic and geographic range, accounting for weights</dd>
                                              <dt>Focal site similarity</dt> <dd>Climatic distance to the site selected on map, based on the selected analog direction</dd>
                                              <dt>Connectivity</dt> <dd>Circuit-based connectivity outputs between historical and future surfaces (available when the connectivity module is enabled); choose the specific connectivity variable below</dd>
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
                        selectInput("time", "Raster variable timeframe", c("historical", "future", "delta")) %>%
                              tooltip("Choose whether to display historical (2000-2025) or future (scenario selected above) model predictions, or their difference (delta).
                                      Only relevant for `FIA variable` and `ESS` raster variables."),
                        selectInput("vector", "Vector variable", vector_choices_base) %>%
                              tooltip(HTML("Select which point or line layer is displayed:<br><br>
                                          <dl style='margin:0;'>
                                              <dt>Occupied FIA plots</dt> <dd>Plots where the selected species was recorded</dd>
                                              <dt>All FIA plots</dt> <dd>All plots including non-forested plots not surveyed</dd>
                                              <dt>Focal site analogs</dt> <dd>Plots within climatic and geographic range of the selected site</dd>
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
                        selectInput("rast_palette", "Color palette",
                                    choices = c("default", "viridis", "magma", "plasma", "inferno",
                                                "RdYlBu", "Spectral", "RdBu", "BrBG"),
                                    selected = "default") %>%
                              tooltip("Override the default palette for the current raster variable. When the raster shows the historical or future FIA variable, the FIA plot markers and focal-site scatterplots adopt the same palette."),
                        input_switch("show_site", "Show focal site", TRUE),
                        input_switch("show_polygon", "Show study area polygon", TRUE)
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

                              ## Named for its contents rather than its input:
                              ## all three plots are views of site_analogs(), the
                              ## analog SET of the focal site. The site itself is
                              ## now set under "Study area" in the left sidebar.
                              accordion_panel(
                                    "Focal site analogs",
                                    helpText(
                                          "This module lets you explore climate analogs for individual sites of interest.
                                          It does not influence the model results shown on the map except for the 'Focal site similarity' raster variable.
                                          Click the map to select a site; its coordinates and a clear button are under 'Study area' in the left sidebar.
                                          The plots below display properties of FIA plots that are both geographically close and climatically similar to your site.
                                          Purple reference lines illustrate bandwidths (dashed) and hard cutoffs (3x bandwidth) for geographic and climatic distances."), hr(),
                                    layout_column_wrap(
                                          width = 1/2,
                                          selectInput("analog_direction", NULL, c("Contemporary analogs", "Reverse analogs", "Forward analogs")) %>%
                                                tooltip(HTML("Choose which type of analog to visualize:
                                                <ul>
                                                    <li>Contemporary analogs have historical climates similar to the focal site's historical climate.</li>
                                                    <li>Reverse (inbound) analogs have historical climates similar to the focal site's future climate.</li>
                                                    <li>Forward (outbound) analogs have future climates similar to the focal site's historical climate.</li>
                                                </ul>")),
                                          input_switch("filter_clim", "Exclude non-analog plots", TRUE) %>%
                                                tooltip("Whether to include sites that are more than 3 climate bandwidths away from the focal site.
                                                        Sites more than 3 geographic bandwidths away are never shown.")
                                    ),

                                    layout_column_wrap(
                                          width = 1/2,
                                          plotOutput("clim_plot"),
                                          plotOutput("dist_plot")
                                    ),
                                    plotOutput("comp_plot") %>%
                                          tooltip("This plot shows the kernel-weighted composition of species or forest types across analog FIA plots.
                                                  Contemporary analogs give the current modeled composition of the focal site.
                                                  Reverse analogs give the potential future composition of the focal site.
                                                  Forward analogs show the current composition of locations that may be suitable for the focal sites's current inhabitants in the future."),
                              ),

                              ## Renamed from "Study area": that name now belongs
                              ## to the left panel that owns the polygon. This
                              ## panel is named for what it renders, matching the
                              ## rule that right-sidebar panels are local to one
                              ## output.
                              accordion_panel(
                                    "Value distribution",
                                    helpText("Compare the distribution of raster values inside your study area against the broader modeling context.
                                             Upload a polygon under 'Study area' in the left sidebar.
                                             The plot shows three curves: the modeling domain (all cells in the current bbox), your study area (cells inside the polygon), and the buffered area
                                             (cells within 3× the geographic bandwidth of the polygon — the predictive neighborhood used for analog matching)."), hr(),

                                    plotOutput("polygon_density"),
                                    selectInput("density_transform", "X-axis scale",
                                                choices = c("identity", "uniform", "sqrt", "log10"),
                                                selected = "identity") %>%
                                          tooltip("Optionally transform the x-axis for skewed variables. 'uniform' rank-transforms the values across the modeling domain, spreading skewed distributions evenly.")
                              ),

                              accordion_panel(
                                    "Climate relationships",
                                    helpText("Use this module to visualize how your selected FIA variable varies across climate gradients.
                                             The plot shows smoothed relationships, based on either FIA plot data or interpolated grid cell values."), hr(),

                                    layout_column_wrap(
                                          width = 1/3,
                                          # cr_v1 / cr_v2 choices are populated dynamically from the
                                          # selected climate variables (see observeEvent(input$clim_vars));
                                          # these initial values are just placeholders.
                                          selectInput("cr_v1", "Variable 1 (x-axis)", c("cwd", "aet")),
                                          selectInput("cr_clim", "Climate timeframe", c("historical", "future")) %>%
                                                tooltip("Which climate the x-axis and grouping variable are drawn from. Use 'future' to pair a future prediction with the future climate that produced it. FIA plot observations are historical, so pairing them with future climate relates an observed response to a climate the plot has not experienced."),
                                          selectInput("cr_data", "Dataset", c("modeled grid cells", "FIA plots")),
                                          selectInput("cr_v2", "Variable 2 (levels)", c("aet", "cwd")),
                                          sliderInput("cr_n", "Number of levels for variable 2", 1, 6, 3, step = 1)
                                    ),


                                    plotOutput("climate_response", height = "600px")
                              ),

                              accordion_panel(
                                    "Model evaluation",
                                    helpText("This module uses cross-validation to assess the accuracy of contemporary AIM predictions for the selected species and FIA variable.
                                             Choose the area the accuracy is measured over: predictions always draw on every available analog plot, but the metrics
                                             can be restricted to plots inside your study area or near your focal site."), hr(),

                                    selectInput("eval_scope", "Evaluation scope",
                                                choices = eval_scope_names,
                                                selected = "domain") %>%
                                          tooltip(HTML("Which held-out plots the accuracy metrics are computed from. The model itself is unchanged.<br><br>
                                                       <dl style='margin:0;'>
                                                           <dt>Modeling domain</dt> <dd>Every plot in the analog pool.</dd>
                                                           <dt>Study area</dt> <dd>Plots inside the uploaded polygon. Often few; check the reported n.</dd>
                                                           <dt>Buffered study area</dt> <dd>Plots within 3 geographic bandwidths of the polygon &mdash; the neighborhood predictions there are actually drawn from. Usually the better choice for a small unit.</dd>
                                                           <dt>Focal site neighborhood</dt> <dd>Plots within the geographic cutoff of the clicked site.</dd>
                                                       </dl>
                                                       Scopes whose input is not yet set under 'Study area' in the left sidebar will say so when selected.")),
                                    uiOutput("eval_stale"),

                                    helpText("The first plot shows predicted versus observed values for the currently-selected model parameters.
                                             Out-of-scope plots are drawn in grey for context."), br(),
                                    plotOutput("eval_scatter"),
                                    uiOutput("eval_scatter_stats"), hr(),

                                    helpText("The tool below helps you explore how different parameter values affect predictive accuracy.
                                             For bandwidth parameters, it tests five values bracketing the value selected under 'Analog bandwidths'.
                                             For climate variables, it tests all combinations of the variables selected under 'Climate data'.
                                             The scatterplot above always reflects your current settings; the heatmap below shows the last completed
                                             run, and is faded and flagged if those settings have since changed."), br(),

                                    layout_column_wrap(
                                          width = 1/2,
                                          checkboxGroupInput("eval_params", NULL,
                                                             choiceNames = c("Test geographic bandwidths", "Test climate bandwidths", "Test climate variables"),
                                                             choiceValues = c("geog", "clim", "vars")),
                                          selectInput("eval_stat", "Evaluation statistic",  c("AUC", "TSS", "Brier", "R-squared", "RMSE")),
                                          sliderInput("eval_n_vars", "Max n climate vars", 1, 6, value = c(1, 2), step = 1),
                                          actionButton("optimal_params", "Apply optimal parameters"),
                                          actionButton("run_eval", "Run evaluations")
                                    ),

                                    plotOutput("eval_heatmap", height = "600px")
                              ),

                              accordion_panel(
                                    "Climate models",
                                    helpText("This module helps you choose among climate models. For every
                                             available scenario and GCM, it plots the projected climate change
                                             (future minus historical) as a trajectory across time periods, so
                                             you can see which models are hot/dry, warm/wet, and so on, and how
                                             the spread fans out over eras and emissions scenarios. Models in your
                                             current ensemble selection (set in the left sidebar) are drawn at full
                                             opacity, others faded. The x and y variables can be changed
                                             without re-clicking Update."), hr(),
                                    layout_column_wrap(
                                          width = 1/2,
                                          selectInput("model_delta_region", "Summarize over",
                                                      choices = c("Modeling domain" = "domain",
                                                                  "Focal site"      = "site",
                                                                  "Study area"      = "polygon"),
                                                      selected = "domain") %>%
                                                tooltip("Region over which deltas are summarized. 'Modeling domain'
                                                        averages over the current bounding box; 'Focal site' samples
                                                        at the clicked site; 'Study area' averages over the uploaded
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
                                    "Connectivity",
                                    helpText("This module contains settings for modeling circuit-based landscape connectivity
                                             between historical (source) and future (target) FIA surfaces.
                                             Click 'enable' to activate options.
                                             To view connectivity results, select 'connectivity' as the raster variable under 'Map layer'.
                                             Julia compiler startup may take a couple minutes the first time you run a connectivity result, but should be fast thereafter."), hr(),
                                    input_switch("cs_enable", "Enable connectivity", FALSE) %>%
                                          tooltip("Adds the connectivity raster and current animation/field vector layers
                                                  to the 'Map layers' menus in the left sidebar."),
                                    conditionalPanel(
                                          condition = "input.cs_enable",
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
      focal_id <- reactive({
            switch(input$track,
                   "Species"      = input$species,
                   "Forest type"  = input$forest_type,
                   "Forest group" = input$forest_group,
                   "All trees"    = "all_trees")
      })

      focal_label <- reactive({
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
            input$stat == "FIA variable" && input$time != "delta"
      })

      # ggplot fill-scale equivalent of the leaflet choose_palette() logic.
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


      ## DATA SETUP -----------------------

      # Padded bbox: includes plots just outside the modeling domain so they
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
      occupied_values <- reactive({
            po <- plot_occ()
            req(isTruthy(input$track))

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

            if(input$vector == "focal site analogs"){
                  pd <- switch(input$analog_direction,
                               "Reverse analogs" = site_analogs()$rev_analogs,
                               "Forward analogs" = site_analogs()$fwd_analogs,
                               "Contemporary analogs" = site_analogs()$bsl_analogs) %>%
                        mutate(x = analog_x, y = analog_y)
                  # Index the per-plot response by analog_index rather than
                  # reading a `response` y-column: pool()$values carries
                  # occ/num/denom (not `response`) for the pooled-ratio
                  # variables, so response_values() is the uniform source.
                  pd$display <- response_values()[pd$analog_index]
                  pd$occupied <- occupied_values()[pd$analog_index]
            } else {
                  # Show the full analog pool, including plots in the padded
                  # ring just outside the modeling bbox — these are still
                  # contributing as analogs to edge cells in the prediction.
                  pd <- plot_occ()
                  pd$display <- response_values()
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
      metric_choices <- function(track) {
            switch(track,
                   "Species"      = c("occurrence", "abundance", "basal area",
                                      "proportion basal area", "regeneration fraction"),
                   "All trees"    = c("occurrence", "abundance", "basal area",
                                      "regeneration fraction"),
                   "Forest type"  = c("occurrence"),
                   "Forest group" = c("occurrence"))
      }
      stage_choices <- function(track, metric) {
            if (track %in% c("Forest type", "Forest group")) return(character(0))
            switch(metric,
                   "occurrence"            = c("seedling", "sapling", "adult", "any"),
                   "abundance"             = c("seedling", "sapling", "adult"),
                   "regeneration fraction" = c("seedling", "sapling"),
                   character(0))  # basal area / proportion basal area: no stage
      }
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
                   "proportion basal area" = "proportion basal area",
                   "regeneration fraction" = switch(stage,
                                                    "seedling" = "seedling fraction",
                                                    "sapling"  = "sapling fraction"))
      }

      ## Repopulate metric choices when the track changes.
      observeEvent(input$track, {
            mc <- metric_choices(input$track)
            updateSelectInput(session, "metric", choices = mc, selected = mc[1])
      }, ignoreInit = TRUE)

      ## Repopulate stage choices when track or metric changes (staged metrics only).
      observeEvent(list(input$track, input$metric), {
            req(input$metric)
            sc <- stage_choices(input$track, input$metric)
            if (length(sc) > 0) {
                  updateSelectInput(session, "stage", choices = sc, selected = sc[1])
            }
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
      ## be exactly that set. Preserve the current pick when it's still available,
      ## otherwise fall back to the first (v1) / second (v2) selected variable.
      observeEvent(input$clim_vars, {
            cv <- input$clim_vars
            req(length(cv) >= 1)
            keep <- function(cur, default) if (isTruthy(cur) && cur %in% cv) cur else default
            updateSelectInput(session, "cr_v1", choices = cv,
                              selected = keep(input$cr_v1, cv[1]))
            updateSelectInput(session, "cr_v2", choices = cv,
                              selected = keep(input$cr_v2, if (length(cv) >= 2) cv[2] else cv[1]))
      })

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
            updateSelectInput(session, "bbox_x", selected = c(bb$xmin, bb$xmax))
            updateSelectInput(session, "bbox_y", selected = c(bb$ymin, bb$ymax))
      })

      ## Snap domain to current map view
      observeEvent(input$window_domain, {
            bb <- input$map_bounds
            updateSelectInput(session, "bbox_x", selected = c(bb$west, bb$east))
            updateSelectInput(session, "bbox_y", selected = c(bb$south, bb$north))
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
      observeEvent(input$cs_enable, {
            on <- isTRUE(input$cs_enable)

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

      # Swap eval_stat choices to the metric set cv_performance() produces
      # for the current response type. Binary (presence-type) responses get
      # ranking/threshold metrics; continuous responses get error metrics.
      observe({
            if(isTruthy(input$fia_var) && grepl("presence", input$fia_var)){
                  updateSelectInput(session, "eval_stat",
                                    choices = c("AUC", "TSS", "Brier"),
                                    selected = "AUC")
            }else{
                  updateSelectInput(session, "eval_stat",
                                    choices = c("R-squared", "RMSE"),
                                    selected = "R-squared")
            }
      })

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
                        color = "#e05c00",
                        weight = 2,
                        dashArray = "5, 5",
                        fill = FALSE,
                        group = "user_polygon"
                  ) %>%
                  addPolygons(
                        data = poly,
                        color = "darkred",
                        weight = 2,
                        fillOpacity = 0.1,
                        fillColor = "darkred",
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

            if (region == "site") {
                  shiny::validate(need(!is.null(site$x) && !is.null(site$y),
                                       "Click the map to select a focal site, or choose a different summary region."))
            } else if (region == "polygon") {
                  shiny::validate(need(isTruthy(polygon_upload_active()),
                                       "Upload a study-area polygon, or choose a different summary region."))
            }

            # Extract ALL climate variables (not just the current axis pair),
            # so the user can re-pick chart axes without triggering a recompute
            # -- the plot just reselects columns from the cached result. One
            # extra cost is reading 6 raster layers per file instead of 2,
            # which is minor next to the per-file I/O.
            all_vars <- clim_vars()
            hst <- rast("data/clim_historical.tif")[[all_vars]]

            # Region geometry resolved once, outside the per-file loop.
            bb_ext <- ext(c(input$bbox_x, input$bbox_y))
            poly_vect <- if (region == "polygon") {
                  terra::vect(st_transform(user_polygon(), crs(hst)))
            } else NULL

            # Summarize a (future - historical) delta SpatRaster to a named
            # numeric vector, one value per climate variable.
            summarize_delta <- function(delta) {
                  if (region == "site") {
                        v <- terra::extract(delta, cbind(site$x, site$y),
                                            method = "bilinear")
                        stats::setNames(as.numeric(v[1, ]), names(delta))
                  } else if (region == "polygon") {
                        v <- terra::extract(delta, poly_vect, fun = "mean",
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
            bindCache(input$model_delta_region,
                      input$bbox_x, input$bbox_y,
                      site$x, site$y,
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

            list(index_hst = build_analog_index(cbind(xy, clim_hst), coord_type = "lonlat",
                                                geog_res_adj = 4),
                 index_fut = build_analog_index(cbind(xy, clim_fut), coord_type = "lonlat",
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

            aim_query <- function(x){
                  suppressWarnings(analogs::analog_impact(
                        x = x,
                        pool = pool()$index_hst,
                        env =  analogs::kernel(max = input$theta_clim * 3, weight = "gaussian", theta = input$theta_clim),
                        geog = analogs::kernel(max = input$theta_geog * 3, weight = "gaussian", theta = input$theta_geog),
                        stat = c("weighted_mean", weight_stats),
                        y = pool()$values[, y_cols, drop = FALSE]))
            }

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

            aim_bsl <- name_outcomes(aim_query(clim()$hst))
            aim_fut <- name_outcomes(aim_query(clim()$fut))

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
                                   opacity = 1, fillOpacity = .5, color = "#4c32a8",
                                   group = "focal_site") %>%
                        addPolygons(lng = r1[,1], lat = r1[,2],
                                    color = "#4c32a8", weight = 2, fill = FALSE,
                                    dashArray = "12, 5",
                                    group = "focal_site") %>%
                        addPolygons(lng = r3[,1], lat = r3[,2],
                                    color = "#4c32a8", weight = 2, fill = FALSE,
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
      display_layer <- reactive({
            req(focal_id(), aim())  # only run once we have data

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
                  r <- switch(input$time,
                              "delta" = aim()$fut$response - aim()$bsl$response,
                              "historical" = aim()$bsl$response,
                              "future" = aim()$fut$response)
            }

            r
      })

      raster_title <- reactive({
            ## legend title ##
            if(input$stat == "ESS"){
                  title <- input$stat
            }else{
                  pre <- paste0(focal_label(), "<br/>")
                  title <- paste0(pre, input$time, "<br/>", tolower(sub(" ", "<br/>", input$fia_var)))
                  if(input$fia_var == "total basal area"){
                        title <- paste0(title, "<br/>(ft2 BA / acre)")
                  }
                  if(input$fia_var == "number of trees"){
                        title <- paste0(title, "<br/>(trees / acre)")
                  }
                  if(input$fia_var == "number of saplings"){
                        title <- paste0(title, "<br/>(saplings / acre)")
                  }
                  if(input$fia_var == "number of seedlings"){
                        title <- paste0(title, "<br/>(seedlings / acre)")
                  }
                  if(input$stat == "connectivity"){
                        title <- paste0(title, " --<br/>", input$cs_stat)
                  }
                  if(input$stat == "focal site similarity"){
                        title <- paste0(gsub(" analogs", "", input$analog_direction), "<br>climate distance<br>from focal site")
                  }
            }

            title
      })

      observe({
            r <- display_layer()
            req(r)

            trans <- switch(input$rast_transform,
                            "identity" = identity,
                            "log10" = log10,
                            "sqrt" = sqrt,
                            "uniform" = function(x) scales::rescale(rank(x, na.last = "keep")))
            r <- setValues(r, trans(r[]))
            r[!is.finite(r[])] <- NA

            proxy <- leafletProxy("map") %>%
                  clearImages() %>%
                  removeControl("raster_legend")

            ## palettes ##
            proj_method <- "bilinear"
            mm <- minmax(r)
            dom <- if(input$rast_transform == "identity") c(0, mm[2]) else mm

            choose_palette <- function(x) if(input$rast_palette == "default") x else input$rast_palette

            if(input$time == "delta" && input$stat == "FIA variable"){
                  # units are deltas
                  r <- display_layer()
                  # Fill non-analog NAs -- except for the regeneration fractions,
                  # whose NAs are a deliberate mask (species not found within
                  # range). Filling those would read as "no change" on the
                  # diverging scale rather than "no information".
                  if(!grepl("fraction", input$fia_var)) r[is.na(r[])] <- 0
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
                  if(grepl("fraction", input$fia_var)){
                        # Regeneration fractions diverge around 0.5 (relict <->
                        # colonizing), so they default to RdBu. Fixing the domain to
                        # the variable's full [0, 1] range keeps 0.5 on the palette
                        # midpoint rather than floating with the observed max.
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
                  proxy <- proxy %>%
                        addLegend(
                              position = "topright",
                              pal      = rast_pal,
                              values   = dom,
                              title    = HTML(raster_title()),
                              layerId  = "raster_legend",
                              opacity  = input$rast_opacity
                        )
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
                  # 0.5 stays on the diverging palette's midpoint.
                  dom <- if(response_is_binary()) 0:1
                  else if(grepl("fraction", input$fia_var)) c(0, 1)
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
                  } else if(grepl("fraction", input$fia_var)){
                        "RdBu"
                  } else {
                        c("gray90", green)
                  }
                  plot_pal <- colorNumeric(palette = marker_palette,
                                           domain = dom)
                  breaks <- pretty(seq(dom[1], dom[2], length.out = 100), n = 5)

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
                  proxy <- proxy %>%
                        addControl(
                              circleLegend(
                                    title  = HTML(paste("FIA plots: historical", input$fia_var, sep = "<br/>")),
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

            # See note in the FIA-plot marker layer: pool()$values has no
            # `response` column for the pooled-ratio variables, so take the
            # per-plot value from response_values() via analog_index.
            d$value <- response_values()[d$analog_index]
            d <- arrange(d, value) %>% mutate(panel = "Climate space")


            # search perimeters
            perims <- rbind(clim_circle(c(s[[v1]], s[[v2]]), input$theta_clim,
                                        c(scl1$mean, scl2$mean), c(scl1$sd, scl2$sd)) %>%
                                  mutate(r = "b"),
                            clim_circle(c(s[[v1]], s[[v2]]), input$theta_clim * 3,
                                        c(scl1$mean, scl2$mean), c(scl1$sd, scl2$sd)) %>%
                                  mutate(r = "a"))

            dom <- if(!response_is_binary()) seq_domain(d$value) else 0:1

            p <- ggplot() +
                  geom_path(data = perims,
                            aes(v1, v2, group = r, linetype = r),
                            color = "#4c32a8", linewidth = .5, alpha = .5) +
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
                             color = "#4c32a8", fill = "#4c32a8", size = 4, alpha = .5) +
                  marker_fill_scale(limits = dom) +
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
                                 color = "#4c32a8", alpha = .5)
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
            d$value <- response_values()[d$analog_index]
            d <- arrange(d, value) %>% mutate(panel = "Distance space")

            dom <- if(!response_is_binary()) seq_domain(d$value) else 0:1

            ymax <- if(input$filter_clim) input$theta_clim * 3.3 else NA
            ggplot(d, aes(geog_dist, env_dist, fill = value)) +
                  geom_vline(xintercept = input$theta_geog, linetype = "longdash", color = "#4c32a8", alpha = .5) +
                  geom_vline(xintercept = input$theta_geog * 3, color = "#4c32a8", alpha = .5) +
                  geom_hline(yintercept = input$theta_clim, linetype = "longdash", color = "#4c32a8", alpha = .5) +
                  geom_hline(yintercept = input$theta_clim * 3, color = "#4c32a8", alpha = .5) +
                  scale_x_continuous(expand = c(0, 0), limits = c(0, input$theta_geog * 3.3)) +
                  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +

                  geom_point(data = filter(d, samp_method_cd == 2), color = "gray70", size = 1.5) +
                  geom_point(data = filter(d, samp_method_cd == 1), shape = 21, color = "black", size = 3) +

                  geom_point(data = tibble::tibble(geog_dist = 0, env_dist = 0,
                                                   shape_label = paste0("selected site,\n",
                                                                        eras()[2], " climate")),
                             aes(x = geog_dist, y = env_dist, shape = shape_label),
                             color = "#4c32a8", fill = "#4c32a8",
                             size = 10, alpha = .5, inherit.aes = FALSE) +

                  marker_fill_scale(limits = dom) +
                  scale_shape_manual(values = 22) +
                  facet_wrap(~ panel) +
                  theme(legend.position = "top") +
                  strip_style +
                  labs(y = "Climate difference from focal site (z-score)",
                       x = "Geographic distance to focal site (km)",
                       shape = NULL)
      })

      ### forest composition --------------------------------
      output$comp_plot <- renderPlot({
            # Pool analog rows from the three analog directions with their weights.
            d <- bind_rows(site_analogs()$rev_analogs %>% mutate(type = "reverse"),
                           site_analogs()$fwd_analogs %>% mutate(type = "forward"),
                           site_analogs()$bsl_analogs %>% mutate(type = "contemporary")) %>%
                  filter(env_dist <= input$theta_clim * 3)

            # Plot-level analog set: one row per (analog plot × analog direction),
            # joined to plot_occ for plot-level attributes (forest_type, forest_group).
            plot_analogs <- plot_occ()[d$analog_index, ] %>%
                  select(x, y, forest_type, forest_group) %>%
                  bind_cols(select(d, env_dist, geog_dist, type)) %>%
                  mutate(weight = dnorm(env_dist, sd = input$theta_clim) *
                               dnorm(geog_dist, sd = input$theta_geog))

            # SPECIES MODE: per-species decomposition using tree-level data.
            # Shows the selected response variable (presence / BA / BA proportion)
            # broken down by species across analog plots.
            if (input$track == "Species") {

                  analogs <- plot_analogs %>%
                        select(x, y, weight, type) %>%
                        left_join(fia_plots_species_all(), by = join_by(x, y),
                                  relationship = "many-to-many")

                  stat <- weighted.mean

                  a <- analogs %>%
                        select(x, y, weight, species, type) %>%
                        tidyr::expand(tidyr::nesting(x, y, weight, type), species) %>%
                        left_join(analogs, by = join_by(x, y, weight, species, type)) %>%
                        mutate(basal_area = ifelse(is.na(basal_area), 0, basal_area),
                               n_trees = ifelse(is.na(n_trees), 0, n_trees),
                               n_saplings = ifelse(is.na(n_saplings), 0, n_saplings),
                               n_seedlings = ifelse(is.na(n_seedlings), 0, n_seedlings)) %>%
                        group_by(x, y, type) %>%
                        mutate(baplot = sum(basal_area),
                               baprop = basal_area / sum(basal_area),
                               baprop = ifelse(is.finite(baprop), baprop, 0)) %>%
                        group_by(species, type) %>%
                        summarize(
                              pres = stat((n_trees + n_saplings) > 0 | n_seedlings > 0, weight), # any life stage
                              pres_seed = stat(n_seedlings > 0, weight),
                              pres_sap = stat(n_saplings > 0, weight),
                              pres_tree = stat(n_trees > 0, weight),
                              batot = stat(basal_area, weight),
                              baprop = stat(baprop, weight),
                              nttot = stat(n_trees, weight),
                              satot = stat(n_saplings, weight),
                              nstot = stat(n_seedlings, weight),
                              .groups = "drop"
                        ) %>%
                        filter(!is.na(species)) %>%
                        # Per-species regeneration fractions, on the same
                        # occurrence basis as the map (see ratio_pieces): the
                        # weighted presence terms above are exactly the occurrence
                        # weights, so the fraction is one over their sum.
                        mutate(seed_frac = ifelse((pres_seed + pres_tree) > 0,
                                                  pres_seed / (pres_seed + pres_tree), NA_real_),
                               sap_frac  = ifelse((pres_sap + pres_tree) > 0,
                                                  pres_sap / (pres_sap + pres_tree), NA_real_))

                  v <- switch(input$fia_var,
                              "presence probability"  = "pres",
                              "seedling presence"      = "pres_seed",
                              "sapling presence"       = "pres_sap",
                              "tree presence"          = "pres_tree",
                              "proportion basal area" = "baprop",
                              "total basal area"      = "batot",
                              "number of trees"       = "nttot",
                              "number of saplings"    = "satot",
                              "number of seedlings"   = "nstot",
                              "seedling fraction"      = "seed_frac",
                              "sapling fraction"       = "sap_frac")
                  a$value <- a[[v]]
                  a$category <- a$species
                  lab <- paste("historical", input$fia_var)

            } else {
                  # TYPE / GROUP / ALL-TREES MODES: weighted prevalence of the
                  # relevant categorical unit (forest type, forest group, or
                  # species) across the analog plots, computed per analog
                  # direction.
                  if (input$track == "Forest type") {
                        a <- plot_analogs %>%
                              group_by(type, category = forest_type) %>%
                              summarize(value = sum(weight), .groups = "drop_last") %>%
                              mutate(value = value / sum(value)) %>%   # normalize within analog direction
                              ungroup() %>%
                              mutate(category = ft_name(category))
                  } else if (input$track == "Forest group") {
                        a <- plot_analogs %>%
                              group_by(type, category = forest_group) %>%
                              summarize(value = sum(weight), .groups = "drop_last") %>%
                              mutate(value = value / sum(value)) %>%
                              ungroup() %>%
                              mutate(category = fg_name(category))
                  } else {
                        # All trees: per-species prevalence (presence-weighted by analog weight)
                        analogs <- plot_analogs %>%
                              select(x, y, weight, type) %>%
                              left_join(fia_plots_species_all(),
                                        by = join_by(x, y),
                                        relationship = "many-to-many") %>%
                              filter(!is.na(species))
                        # Total analog weight per direction (denominator)
                        totals <- plot_analogs %>%
                              group_by(type) %>%
                              summarize(total_weight = sum(weight), .groups = "drop")
                        # Weight of analog plots where each species is present
                        a <- analogs %>%
                              distinct(x, y, weight, type, species) %>%   # one row per (plot × species × direction)
                              group_by(type, category = species) %>%
                              summarize(value = sum(weight), .groups = "drop") %>%
                              left_join(totals, by = "type") %>%
                              mutate(value = value / total_weight) %>%
                              select(-total_weight)
                  }

                  # Expand to every (category × analog direction) combination,
                  # filling categories absent from a direction with 0. Without
                  # this, missing combinations render as blank tiles rather
                  # than explicit zeros.
                  a <- a %>%
                        tidyr::complete(category, type, fill = list(value = 0))

                  # Drop categories with zero weighted prevalence everywhere
                  a <- a %>%
                        group_by(category) %>%
                        filter(any(value > 0)) %>%
                        ungroup()

                  lab <- "historical weighted prevalence"
            }

            req(nrow(a) > 0)

            a <- a %>%
                  arrange(category) %>%
                  mutate(svalue = value / max(value, na.rm = TRUE))

            b <- a %>%
                  select(type, category, value) %>%
                  mutate(value = value / max(value),
                         type = factor(type, levels = c("reverse", "contemporary", "forward")),
                         itype = as.integer(type) - 0.4,
                         itype2 = itype + 0.8,
                         icat = as.integer(factor(category))) %>%
                  group_by(type) %>%
                  mutate(pad = .05,
                         ymax = cumsum(value) + cumsum(pad) - pad,
                         ymin = lag(ymax) + pad,
                         ymin = ifelse(is.na(ymin), 0, ymin),
                         y = ymin/2 + ymax/2) %>%
                  tidyr::gather(rep, itype, itype, itype2)

            labs <- b %>%
                  group_by(category) %>%
                  filter(value == max(value)) %>%
                  summarize(itype = mean(itype),
                            y = mean(y))

            ggplot() +
                  geom_ribbon(data = b,
                              aes(itype, ymin = ymin, ymax = ymax,
                                  group = category, fill = category, color = category),
                              alpha = .7, linewidth = 0.05) +
                  geom_rect(data = b %>% tidyr::spread(rep, itype),
                            aes(xmin = itype, xmax = itype2, ymin = ymin, ymax = ymax,
                                group = category, fill = category, color = category)) +
                  geom_text(data = labs,
                            aes(itype, y, label = category),
                            hjust = .5, vjust = .5, fontface = "bold.italic") +
                  scale_x_continuous(breaks = 1:3,
                                     labels = c("reverse\nanalogs", "contemporary\nanalogs", "forward\nanalogs")) +
                  theme(legend.position = "none",
                        panel.grid = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank()) +
                  labs(y = "weighted composition of analog plots")
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

            domain_vals <- na.omit(values(display_layer(), mat = FALSE))

            # Polygon-dependent values: NULL when no polygon is uploaded. The
            # corresponding curves (study area, buffer) are then simply absent
            # from the plot. Check the upload input directly to avoid
            # triggering user_polygon()'s internal req() short-circuit.
            has_polygon <- isTruthy(polygon_upload_active())
            poly_vals   <- if (has_polygon) polygon_raster_values() else NULL
            buffer_vals <- if (has_polygon) polygon_buffered_raster_values() else NULL

            study_color   <- "darkred"   # matches the map polygon
            buffer_color  <- "#e05c00"   # matches the buffer outline on the map
            domain_color  <- "gray60"
            focal_color   <- "#4c32a8"   # matches the focal site marker elsewhere

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
            x_lab <- gsub("<br/>", " ", raster_title())
            if (input$density_transform == "uniform") {
                  ecdf_ref <- ecdf(domain_vals)
                  if (has_polygon) {
                        poly_vals   <- ecdf_ref(poly_vals)
                        buffer_vals <- ecdf_ref(buffer_vals)
                  }
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
                  tibble::tibble(value = domain_vals, source = "Modeling domain")
            )
            if (has_polygon) {
                  parts <- c(parts, list(
                        tibble::tibble(value = buffer_vals, source = "Buffered area"),
                        tibble::tibble(value = poly_vals,   source = "Study area")
                  ))
            }
            d <- bind_rows(parts) %>%
                  mutate(source = factor(source,
                                         levels = c("Modeling domain",
                                                    "Buffered area",
                                                    "Study area",
                                                    "Focal site")))

            pal <- c("Modeling domain" = domain_color,
                     "Buffered area"   = buffer_color,
                     "Study area"      = study_color,
                     "Focal site"      = focal_color)
            ltypes <- c("Modeling domain" = "solid",
                        "Buffered area"   = "dashed",
                        "Study area"      = "solid",
                        "Focal site"      = "solid")

            panel_levels <- c("Relative density", "ECDF")
            factor_panel <- function(df, lvl) {
                  df$panel <- factor(lvl, levels = panel_levels)
                  df
            }

            # Optional focal-site reference line. Skipped (NULL → no-op when
            # added to ggplot) if no focal site is selected or the cell value
            # isn't finite. Built as a data-driven geom_vline so the color and
            # linetype scales pick it up for the legend.
            focal_layer <- if (!is.na(focal_val)) {
                  geom_vline(data = tibble::tibble(xint = focal_val,
                                                   source = factor("Focal site",
                                                                   levels = c("Modeling domain",
                                                                              "Buffered area",
                                                                              "Study area",
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
            legend_breaks <- c("Modeling domain", "Buffered area", "Study area")
            if (!is.na(focal_val)) legend_breaks <- c(legend_breaks, "Focal site")

            p <- ggplot(d, aes(x = value, color = source,
                               fill = source, linetype = source)) +
                  geom_density(aes(group = source, y = after_stat(scaled)),
                               data = ~ factor_panel(., "Relative density"),
                               alpha = 0.25, linewidth = 0.6,
                               bounds = d_bounds) +
                  stat_ecdf(aes(group = source),
                            data = ~ factor_panel(., "ECDF"),
                            geom = "step", linewidth = 0.7, pad = FALSE) +
                  focal_layer +
                  facet_wrap(~ panel, nrow = 1) +
                  scale_y_continuous(
                        breaks = c(0, 0.25, 0.5, 0.75, 1),
                        expand = expansion(0, 0)
                  ) +
                  scale_fill_manual(values = c("Modeling domain" = domain_color,
                                               "Buffered area"   = buffer_color,
                                               "Study area"      = study_color,
                                               "Focal site"      = NA),
                                    breaks = legend_breaks,
                                    guide = "none") +
                  scale_color_manual(values = pal,
                                     breaks = legend_breaks) +
                  scale_linetype_manual(values = ltypes,
                                        breaks = legend_breaks) +
                  theme(legend.position = "top") +
                  strip_style +
                  labs(x = x_lab,
                       y = NULL,
                       fill = NULL, color = NULL, linetype = NULL)

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
                   "polygon"  = "Upload a study area polygon under 'Study area' in the left sidebar to evaluate over it.",
                   "buffered" = "Upload a study area polygon under 'Study area' in the left sidebar to evaluate over its buffered neighborhood.",
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
            # derive the response from that same object rather than calling
            # response_values() (which re-derives from plot_occ() and could,
            # in principle, see a different value).
            po <- plot_occ()
            xy <- as.data.frame(po[, c("x", "y")])
            clim_hst <- terra::extract(clim_r, xy, ID = FALSE, method = "bilinear")
            d <- cbind(xy, clim_hst)
            values <- response_values()

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
                  cv <- suppressWarnings(analogs::analog_cv(
                        fun       = analogs::analog_impact,
                        pool      = d,
                        y         = values,
                        cv_method = "loo",
                        env = analogs::kernel(max = theta_clim * 3, weight = "gaussian", theta = theta_clim),
                        geog = analogs::kernel(max = theta_geog * 3, weight = "gaussian", theta = theta_geog),
                        stat      = "weighted_mean"))

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
                                 "No results for the selected evaluation scope."))

            updateSliderInput(session, "theta_geog", value = opt$theta_geog)
            updateSliderInput(session, "theta_clim", value = opt$theta_clim)
            updateSelectizeInput(session, "clim_vars", selected = unlist(strsplit(opt$clim_vars, "_")))
      })

      output$eval_heatmap <- renderPlot({
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
                  else sprintf("No evaluation results for '%s'. Either the last run predates this selection, or fewer than %d plots fall inside it.",
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
                                 sprintf("'%s' is not defined for this response within '%s'.",
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
                       fill = paste0(input$eval_stat, " (", input$fia_var, ")"))

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
            e <- eval_scatter_data()
            s <- eval_scatter_scope()
            req(e, s)

            d <- tibble::tibble(obs      = e$cv$obs,
                                pred     = e$cv$obs - e$cv$residual,
                                in_scope = s$keep) %>%
                  filter(!is.na(obs), !is.na(pred))

            focus <- d[d$in_scope, , drop = FALSE]

            shiny::validate(need(
                  nrow(focus) >= eval_min_n,
                  sprintf("Only %d plot(s) fall within '%s' \u2014 too few to read a fit from.",
                          nrow(focus), s$label)))

            # The domain cloud stays visible underneath. Seeing that a unit's
            # plots sit in the tail of the domain's error distribution says more
            # than seeing those plots alone.
            p <- ggplot(mapping = aes(obs, pred))
            if (!all(d$in_scope)) {
                  p <- p + geom_point(data = d[!d$in_scope, , drop = FALSE],
                                      color = "gray80", size = 1)
            }
            p +
                  geom_point(data = focus) +
                  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "dodgerblue") +
                  geom_smooth(data = focus, method = lm, se = F, color = "red", formula = 'y ~ x') +
                  annotate(geom = "point", color = "red", size = 4,
                           x = mean(focus$obs), y = mean(focus$pred)) +
                  scale_x_sqrt() +
                  scale_y_sqrt() +
                  coord_fixed() +
                  theme_bw() +
                  labs(x = paste0("observed ", input$fia_var),
                       y = paste0("predicted ", input$fia_var),
                       subtitle = if (all(d$in_scope)) s$label
                       else sprintf("%s \u2014 %d of %d plots",
                                    s$label, nrow(focus), nrow(d)))
      })

      # Held-out accuracy for the current settings, over the selected scope.
      # Reports the whole metric set cv_performance() returns rather than only
      # input$eval_stat: for a single configuration there is no cost to showing
      # all of them, and they answer different questions (bias and R-squared
      # disagree in useful ways). Deliberately NOT faded by gg_stale -- this
      # section is computed from current settings, so unlike the heatmap it is
      # never stale, and that contrast is the point.
      output$eval_scatter_stats <- renderUI({
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

            # cr_v1/cr_v2 track input$clim_vars via an observer; guard against the
            # brief flush after a climate-variable change where the previous pick
            # isn't yet in the layer set.
            req(input$cr_v1 %in% input$clim_vars, input$cr_v2 %in% input$clim_vars)

            v1 <- input$cr_v1
            v2 <- input$cr_v2
            n_groups <- input$cr_n


            # Climate for the x-axis and grouping variable. Historical by default;
            # 'future' lets you pair a future prediction with the future climate
            # that produced it, rather than reading it against current climate.
            use_fut <- isTRUE(input$cr_clim == "future")

            if(input$cr_data == "modeled grid cells"){
                  cr <- if (use_fut) clim()$fut else clim()$hst
                  d <- rast(c(y = display_layer(),
                              clim1 = cr[[v1]],
                              clim2 = cr[[v2]])) %>%
                        as.data.frame()
            }

            if(input$cr_data == "FIA plots"){
                  # response_values() is the per-plot value for every variable
                  # (including the pooled-ratio variables, whose pool()$values
                  # carries occ/num/denom rather than a `response` column). It is
                  # the raw historical plot observation.
                  cp <- if (use_fut) pool()$clim_fut else pool()$clim_hst
                  d <- data.frame(y = response_values(),
                                  clim1 = cp[[v1]],
                                  clim2 = cp[[v2]])
            }

            # Climate rasters are stored standardized, so back-convert to native
            # units (z * sd + mean) before grouping and plotting -- this puts both
            # the x-axis and the legend's interval labels in mm / degrees C. Falls
            # back to leaving a variable standardized if it has no scale entry.
            back_transform <- function(x, v) {
                  s <- clim_scales[[v]]
                  if (is.null(s) || !is.finite(s$sd)) return(x)
                  x * s$sd + s$mean
            }
            d$clim1 <- back_transform(d$clim1, v1)
            d$clim2 <- back_transform(d$clim2, v2)

            # Group clim2 into equal-COUNT (quantile) bins rather than equal-width
            # bins: skewed climate variables leave equal-width tail bins too sparse
            # (or empty) to fit a smooth, which drops groups from the plot. Fall
            # back to equal-width if the data can't support n distinct quantiles.
            if(n_groups > 1) d$group <- tryCatch(cut_number(d$clim2, n_groups),
                                                 error = function(e) cut(d$clim2, n_groups))
            if(n_groups == 1) d$group <- "All values"

            # Smooth on a link scale so fitted curves respect each variable's
            # natural bounds: logit for [0, 1] variables, log for non-negative
            # variables, identity for differences. Uses a GLM with a natural-spline
            # basis; quasi-families accept the continuous grid-cell values without
            # integer/proportion warnings. The FIA-plots dataset is the raw
            # historical per-plot variable, so its family depends on the variable
            # alone; only the modeled-grid dataset honors the raster stat and time.
            bounded01 <- grepl("presence|fraction", input$fia_var) ||
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
            fits <- unlist(lapply(split(dd, dd$group), function(g) {
                  tryCatch(fitted(glm(y ~ splines::ns(clim1, 4), family = fam, data = g)),
                           error = function(e) numeric(0))
            }))
            ylim <- if (length(fits)) range(fits, na.rm = TRUE) else NULL

            p_main <- ggplot(d, aes(clim1, y, color = group, fill = group)) +
                  geom_smooth(method = "glm",
                              formula = y ~ splines::ns(x, 4),
                              method.args = list(family = fam)) +
                  scale_fill_viridis_d() +
                  scale_color_viridis_d() +
                  coord_cartesian(ylim = ylim) +
                  # Label the climate axes with their timeframe (so historical vs
                  # future is never ambiguous) and their native units.
                  labs(x = paste0(input$cr_clim, " ", v1, " (", unit_lab(v1), ")"),
                       fill = paste0(input$cr_clim, "\n", v2, " (", unit_lab(v2), ")"),
                       color = paste0(input$cr_clim, "\n", v2, " (", unit_lab(v2), ")"),
                       # FIA-plots data is raw historical observation, so label it
                       # by the variable rather than the (stat/time-aware) raster title.
                       y = if (input$cr_data == "FIA plots") paste("observed", input$fia_var)
                       else gsub("<br/>", " ", raster_title()))

            # Marginal density of the x variable, split by group. This shows each
            # clim2 bin's coverage along clim1, which is what explains sparse-tail
            # extrapolation in the smooths below. It carries its own y-scale, so it
            # doesn't interact with the fitted-line ylim on the main panel.
            p_dens <- ggplot(d, aes(clim1, fill = group, color = group)) +
                  geom_density(alpha = .3, linewidth = .3, position = "stack") +
                  scale_fill_viridis_d() +
                  scale_color_viridis_d() +
                  theme_void() +
                  theme(legend.position = "none")

            # patchwork aligns the two panels on a shared x-axis (both are built
            # from the same clim1), and collects the single legend from p_main.
            p_dens / p_main +
                  plot_layout(heights = c(1, 8), guides = "collect")
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
