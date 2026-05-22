library(jsonlite)
library(shiny)
library(leaflet)
library(ggplot2)
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

cs_setup()

clim_vars <- function() c("tmean", "ppt", "aet", "cwd", "tmincm", "tmaxwm")
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

.fia_global <- read_csv("data/fia.csv", show_col_types = FALSE)
spp <- sort(unique(.fia_global$species))
forest_types_in_data <- sort(unique(na.omit(.fia_global$forest_type)))
# Initial map bounds: bounding box of all FIA plots
fia_xy_bounds <- list(
      xmin = min(.fia_global$x, na.rm = TRUE),
      xmax = max(.fia_global$x, na.rm = TRUE),
      ymin = min(.fia_global$y, na.rm = TRUE),
      ymax = max(.fia_global$y, na.rm = TRUE)
)
rm(.fia_global)

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

theme_set(theme_bw(base_size = 16))


# HELPER FUNCTIONS ------------------------------------

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

auc <- function(obs, pred){
      if(!all(obs %in% 0:1)) return(NA)
      as.vector(pROC::roc(obs, pred, quiet = TRUE)$auc)
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

      theme = bs_theme(version = 5, bootswatch = "sandstone"),

      fillable = TRUE,   # lets the page body fill the browser height

      title = "TreeShift",


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
                        p("This tool uses 'analog impact models' to interpolate FIA variables across
                          a spatial grid, based on geographic proximity and climatic similarity to
                          FIA plots, optionally accounting for pedicted future climate change.
                          Use it to assess potential shifts in species distributions and other forest characteristics,
                          model connectivity between current and potential future distributions,
                          and explore properties of plots that are analogs to your focal sites.",
                          style = "font-size: 14px; color: #555;"),
                        br(),
                        p("NOTE: This tool is a prototype and is still under development.
                          Please contact matthew.kling@umontana.edu with questions or bug reports.",
                          style = "font-size: 14px; color: #555;")
                  ),

                  accordion_panel(
                        "FIA variable",
                        selectInput("track", "Mode",
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
                        selectInput("fia_var", "Variable",
                                    c("presence probability", "proportion basal area", "total basal area")) %>%
                              tooltip("Choose which focal variable gets interpolated and predicted.")
                  ),

                  accordion_panel(
                        "Spatial grid",
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
                        sliderInput("clim_agg", "Raster grain size", 1, 10, value = 10, step = 1, post = " km") %>%
                              tooltip("Select the resolution of the modeling grid; larger values are faster but less precise.")
                  ),

                  accordion_panel(
                        "Climate data",
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
                        selectInput("scenario_model", "Climate model",
                                    choices = model_choices,
                                    selected = default_model) %>%
                              tooltip("Global climate model (GCM).")
                  ),

                  accordion_panel(
                        "Analog settings",
                        helpText("Analogs are weighted using Gaussian kernels based on climatic and geographic distances from a site.
                                 Bandwidths specify the standard deviations of these kernels. Plots farther than 3 stdevs
                                 are considered non-analog and excluded."), hr(),
                        sliderInput("theta_geog", "Geographic bandwidth (km)", 5, 100, value = 30),
                        sliderInput("theta_clim", "Climate bandwidth (z)", .05, .5, value = .1) %>%
                              tooltip("Climate bandwidth units are in standard deviations relative to climate variation across all FIA plots in the contiguous US"),
                        selectInput("summary_stat", "Summary statistic", c("weighted mean", "weighted sum")) %>%
                              tooltip(HTML("How to summarize the focal variable across analog FIA plots.<br><br>
                                        <dl style='margin:0;'>
                                            <dt>Weighted mean</dt>
                                            <dd>Expected conditions given climate and location. Can be high even when analogs are few and distant.</dd>
                                            <dt>Weighted sum</dt>
                                            <dd>Propagule pressure potential: suitability weighted by accessibility to currently-occupied analog plots.</dd>
                                        </dl>"))
                  ),

                  accordion_panel(
                        "Connectivity",
                        helpText("Settings for circuit-based connectivity modeling between historical (source)
                                 and future (target) FIA surfaces. To view connectivity results, select
                                 'connectivity' as the raster variable under 'Map layer'."), hr(),
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
                              tooltip("Relative strength of loss conductance, compared to source, target, and horizontal components.")
                  ),

                  accordion_panel(
                        "Map layer",
                        helpText("Choose which model results are drawn on the map."), hr(),
                        selectInput("stat", "Raster variable", c("FIA variable", "effective sample size" = "ESS",
                                                                 "focal site similarity", "connectivity")) %>%
                              tooltip(HTML("<dl style='margin:0;'>
                                              <dt>FIA variable</dt> <dd>Interpolated values for the plot variable selected above</dd>
                                              <dt>Effective sample size</dt> <dd>Effective number of FIA plots within climatic and geographic range, accounting for weights</dd>
                                              <dt>Focal site similarity</dt> <dd>Climatic distance to the site selected on map, based on the selected analog direction</dd>
                                              <dt>Connectivity</dt> <dd>Circuit-based connectivity outputs between historical and future surfaces; choose the specific connectivity variable below</dd>
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
                              tooltip("Choose whether to display historical or future model predictions, or their difference.
                                      Only relevant for `FIA variable` and `ESS` raster variables."),
                        hr(),
                        selectInput("vector", "Vector layer",
                                    c("occupied FIA plots", "all FIA plots", "focal site analogs", "current animation", "current field", "none")) %>%
                              tooltip(HTML("Select which point or line layer is displayed:<br><br>
                                          <dl style='margin:0;'>
                                              <dt>Occupied FIA plots</dt> <dd>Plots where the selected species was recorded</dd>
                                              <dt>All FIA plots</dt> <dd>All plots including non-forested plots not surveyed</dd>
                                              <dt>Focal site analogs</dt> <dd>Plots within climatic and geographic range of the selected site</dd>
                                              <dt>Current animation</dt> <dd>Animation with particle movement showing direction and magnitude of mean current flow</dd>
                                              <dt>Current field</dt> <dd>Vector field showing showing direction and magnitude of mean current flow</dd>
                                          </dl>"))
                  ),

                  accordion_panel(
                        "Cartography",
                        sliderInput("vector_scale", "Point/line size", -1, 1, 0, step = .1) %>%
                              tooltip("Use this slider to adjust the size of FIA plot markers, or the length of current lines"),
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
                        helpText("Download the raster data displayed on the map as a .tif, or the current model settings as a .txt for reproducibility."), hr(),
                        textInput("download_filename", "Filename (without extension)", "fia_aim_result"),
                        downloadButton("download_raster", "Download raster"),
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

                              accordion_panel(
                                    "Focal site",
                                    helpText(
                                          "This module lets you explore climate analogs for individual sites of interest.
                                          It does not influence the model results shown on the map except for the 'Focal site similarity' raster variable.
                                          Click the map to select a site; the plots below display properties of FIA plots that are both geographically close and climatically similar to your site.
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
                                          div(style = "display: flex; gap: 0.5rem;",
                                              actionButton("highlight_site", "Highlight on map",
                                                           style = "flex: 1;"),
                                              actionButton("site_clear", "Clear") %>%
                                                    tooltip("Remove the focal site from the map and plots.")),
                                          input_switch("filter_clim", "Exclude non-analog plots", TRUE) %>%
                                                tooltip("Whether to include sites that are more than 3 climate bandwidths away from the focal site.
                                                        Sites more than 3 geographic bandwidths away are never shown.")
                                    ),

                                    layout_column_wrap(
                                          width = 1/2,
                                          plotOutput("clim_plot"),
                                          plotOutput("dist_plot")
                                    ),
                                    plotOutput("comp_plot")
                              ),

                              accordion_panel(
                                    "Study area",
                                    helpText("Upload a polygon layer to display on the map and compare the distribution of raster values inside your study area against the broader modeling context.
                                             The plot shows three curves: the modeling domain (all cells in the current bbox), your study area (cells inside the polygon), and the buffered area
                                             (cells within 3× the geographic bandwidth of the polygon — the predictive neighborhood used for analog matching)."), hr(),

                                    layout_column_wrap(
                                          width = NULL,
                                          style = htmltools::css(grid_template_columns = "1fr auto"),
                                          fileInput("polygon_upload", NULL,
                                                    multiple = TRUE,
                                                    accept = c(".shp", ".dbf", ".shx", ".prj", ".kml", ".kmz",
                                                               ".gpkg", ".geojson"),
                                                    placeholder = "Shapefile, KML, GPKG, or GeoJSON") %>%
                                                tooltip("Upload a polygon layer. For shapefiles, select all component files (.shp, .dbf, .shx, .prj)."),
                                          div(style = "padding-top: 0.5rem;",
                                              actionButton("polygon_clear", "Clear") %>%
                                                    tooltip("Remove the polygon from the map and plot."))
                                    ),

                                    plotOutput("polygon_density"),
                                    selectInput("density_transform", "X-axis scale",
                                                choices = c("identity", "uniform", "sqrt", "log10"),
                                                selected = "identity") %>%
                                          tooltip("Optionally transform the x-axis for skewed variables. 'uniform' rank-transforms the values across the modeling domain, spreading skewed distributions evenly.")
                              ),

                              accordion_panel(
                                    "Model evaluation",
                                    helpText("This module uses cross-validation to assess the accuracy of contemporary AIM predictions for the selected species and FIA variable.
                                             The first plot shows predicted versus observed values for the currently-selected model parameters."), br(),
                                    plotOutput("eval_scatter"), hr(),

                                    helpText("The tool below helps you explore how different parameter values affect predictive accuracy.
                                             For bandwidth parameters, it tests five values bracketing the value selected under 'Analog bandwidths'.
                                             For climate variables, it tests all combinations of the variables selected under 'Climate data', using the selected grain size."), br(),

                                    layout_column_wrap(
                                          width = 1/2,
                                          checkboxGroupInput("eval_params", NULL,
                                                             choiceNames = c("Test geographic bandwidths", "Test climate bandwidths", "Test climate variables"),
                                                             choiceValues = c("geog", "clim", "vars")),
                                          selectInput("eval_stat", "Evaluation statistic",  c("R-squared", "RMSE", "delta AIC", "AUC")),
                                          sliderInput("eval_n_vars", "Max n climate vars", 1, 6, value = c(1, 2), step = 1),
                                          actionButton("optimal_params", "Apply optimal parameters"),
                                          actionButton("run_eval", "Run evaluations")
                                    ),

                                    plotOutput("eval_heatmap", height = "600px")
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

      # TRUE if the active response is a 0/1 indicator (presence variables);
      # FALSE if continuous (basal area variables).
      response_is_binary <- reactive({
            isTruthy(input$fia_var) && input$fia_var == "presence probability"
      })

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
            read_csv("data/fia.csv", show_col_types = FALSE) %>%
                  select(x, y) %>%
                  distinct() %>%
                  filter(between(x, pb$x[1], pb$x[2]),
                         between(y, pb$y[1], pb$y[2])) %>%
                  mutate(index = 1:nrow(.))
      })

      fia_plots_species_all <- reactive({
            read_csv("data/fia.csv", show_col_types = FALSE) %>%
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
                        pres = as.integer(any(focal)), # focal species presence/absence
                        batot = sum(basal_area[focal]), # total focal species basal area
                        baprop = weighted.mean(focal, basal_area), #  proportion ba
                        baplot = sum(basal_area), # total plot ba, all species
                        forest_type = dplyr::first(forest_type), # plot-level FORTYPCD
                        .groups = "drop") %>%
                  left_join(ft_to_fg, by = "forest_type") %>%
                  mutate(baprop = ifelse(baplot == 0, 0, baprop))
      })

      # Single active response vector, aligned with plot_occ() rows.
      # Computed from the (track, fia_var, focal_id) tuple.
      response_values <- reactive({
            po <- plot_occ()
            req(isTruthy(input$track), isTruthy(input$fia_var))

            v <- if (input$track == "Species") {
                  switch(input$fia_var,
                         "presence probability"  = po$pres,
                         "proportion basal area" = po$baprop,
                         "total basal area"      = po$batot)
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
                         "presence probability" = as.integer(po$baplot > 0),
                         "total basal area"     = po$baplot)
            }

            # Guard against transient NULL during track/fia_var transitions
            req(!is.null(v))
            v
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

      ## Swap FIA variable choices when track changes
      observeEvent(input$track, {
            choices <- switch(input$track,
                              "Species"      = c("presence probability",
                                                 "proportion basal area",
                                                 "total basal area"),
                              "Forest type"  = c("presence probability"),
                              "Forest group" = c("presence probability"),
                              "All trees"    = c("presence probability",
                                                 "total basal area"))
            updateSelectInput(session, "fia_var",
                              choices = choices,
                              selected = choices[1])
      }, ignoreInit = TRUE)

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

      # Only show AUC option for binary (presence-type) variables
      observe({
            if(isTruthy(input$fia_var) && input$fia_var == "presence probability"){
                  updateSelectInput(session, "eval_stat",
                                    choices = c("R-squared", "RMSE", "delta AIC", "AUC"),
                                    selected = "AUC")
            }else{
                  updateSelectInput(session, "eval_stat",
                                    choices = c("R-squared", "RMSE", "delta AIC"))
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


      ## ANALOG IMPACT MODEL -----------------------------------------

      # Raw (uncropped, unaggregated) historical and future climate rasters,
      # restricted to the selected climate variables. The future raster depends
      # on the scenario selection; the historical raster is shared across
      # scenarios. Used as the source for both the modeling grid (via clim())
      # and for plot-level climate extraction (in pool(), site_analogs(), and
      # eval_var_set()).
      # Active scenario filename, reconstructed from the three scenario inputs.
      # If the corresponding file doesn't exist on disk (partial factorial),
      # surfaces a user-visible message via validate() rather than silently
      # failing — downstream outputs (the map, the focal-site plots, etc.)
      # will display this message instead of erroring.
      scenario_file <- reactive({
            req(isTruthy(input$scenario_ssp),
                isTruthy(input$scenario_period),
                isTruthy(input$scenario_model))
            f <- paste0(input$scenario_ssp, "_",
                        input$scenario_period, "_",
                        input$scenario_model, ".tif")
            path <- file.path("data/clim_future", f)
            shiny::validate(need(file.exists(path),
                                 "Selected scenario combination is not available. Try a different emissions pathway, era, or climate model."))
            path
      })

      clim_pool <- reactive({
            # Resolve scenario_file() at the top level so a req() short-circuit
            # (e.g., user selected a (ssp, period, model) combination whose
            # file doesn't exist) propagates as a clean reactive cancellation
            # rather than an rast() argument-evaluation error.
            sf <- scenario_file()
            r_hst <- rast("data/clim_historical.tif")[[input$clim_vars]]
            r_fut <- rast(sf)[[input$clim_vars]]
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


      pool <- reactive({
            req(plot_occ())
            xy <- plot_occ()[, c("x", "y")]
            # Extract climate from the raw (uncropped, unaggregated) rasters at
            # each plot's actual location. This gives the climate the plot is
            # actually experiencing rather than the spatial average over a
            # coarsened cell, and means pool() doesn't depend on the modeling
            # grid (input$clim_agg) or the bbox crop.
            cp <- clim_pool()
            clim_hst <- extract(cp$hst, xy, ID = FALSE, method = "bilinear")
            clim_fut <- extract(cp$fut, xy, ID = FALSE, method = "bilinear")
            list(index_hst = build_analog_index(cbind(xy, clim_hst), coord_type = "lonlat", index_res = 8),
                 index_fut = build_analog_index(cbind(xy, clim_fut), coord_type = "lonlat", index_res = 8),
                 clim_hst = clim_hst,
                 clim_fut = clim_fut,
                 # Single active response column for AIM smoothing,
                 # plus samp_method_cd carried through for plot styling.
                 values = tibble::tibble(response = response_values(),
                                         samp_method_cd = plot_occ()$samp_method_cd))
      })

      ### full raster -----------------------
      aim <- reactive({

            stat <- gsub(" ", "_", input$summary_stat)

            aim_query <- function(x){
                  suppressWarnings(analogs::analog_impact(
                        x = x,
                        pool = pool()$index_hst,
                        max_clim = input$theta_clim * 3,
                        max_geog = input$theta_geog * 3,
                        stat = c("count", "sum_weights", stat, "ess"),
                        kernel = "gaussian_joint",
                        theta = c(input$theta_clim, input$theta_geog),
                        y = pool()$values))
            }

            aim_bsl <- aim_query(clim()$hst)
            aim_fut <- aim_query(clim()$fut)
            names(aim_bsl) <- names(aim_fut) <- gsub("weighted_mean_|weighted_sum_", "", names(aim_fut))

            list(fut = aim_fut,
                 bsl = aim_bsl)
      })

      ### focal site ------------------------
      site_analogs <- reactive({
            req(site$x, site$y)

            xy <- t(as.matrix(c(x = site$x, y = site$y), nrow = 1))
            # Extract from raw rasters to match pool() (consistent scaling).
            cp <- clim_pool()
            hst <- cbind(xy, extract(cp$hst, xy, method = "bilinear"))
            fut <- cbind(xy, extract(cp$fut, xy, method = "bilinear"))

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
                        max_geog = input$theta_geog * 3,
                        max_clim = max_clim,
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

            # Ground nodes
            grounds <- (electrodes()$destinations + leakage()) %>% setNames("grounds")

            # Run circuitscape
            cs <- cs_advanced(resistance()$resistance, electrodes()$sources, grounds,
                              resistance_is = "resistances", ground_is = "conductances")

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
                  title <- paste0(pre, input$time, "<br/>", input$summary_stat, "<br/>", tolower(sub(" ", "<br/>", input$fia_var)))
                  if(input$fia_var == "total basal area" &&
                     input$summary_stat == "weighted mean"){
                        title <- paste0(title, "<br/>(m2 BA / plot)")
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
                  r[is.na(r[])] <- 0 # fill non-analog NAs
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
                  rast_pal <- colorNumeric(
                        palette = c("gray90", green) %>% choose_palette(),
                        domain = dom,
                        na.color = "transparent"
                  )
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

                  if(input$vector == "focal site analogs"){
                        pd <- switch(input$analog_direction,
                                     "Reverse analogs" = site_analogs()$rev_analogs,
                                     "Forward analogs" = site_analogs()$fwd_analogs,
                                     "Contemporary analogs" = site_analogs()$bsl_analogs) %>%
                              mutate(x = analog_x, y = analog_y)
                        pd$display <- pd$response
                  } else {
                        # Show the full analog pool, including plots in the padded
                        # ring just outside the modeling bbox — these are still
                        # contributing as analogs to edge cells in the prediction.
                        pd <- plot_occ()
                        pd$display <- response_values()
                  }

                  if(input$vector == "occupied FIA plots") pd <- filter(pd, display > 0)

                  dom <- if(!response_is_binary()) c(0, max(pd$display, na.rm = TRUE)) else 0:1

                  # Markers follow the raster palette only when the raster is
                  # itself showing FIA response on a sequential scale (see
                  # marker_follows_raster_palette); otherwise the canonical
                  # gray90 -> green ramp. colorNumeric accepts either a color
                  # vector or a named-palette string.
                  marker_palette <- if(input$rast_palette != "default" &&
                                       marker_follows_raster_palette()){
                        input$rast_palette
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
                              data = filter(pd, samp_method_cd == 1),
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

            d$value <- d$response
            d <- arrange(d, value) %>% mutate(panel = "Climate space")


            # search perimeters
            perims <- rbind(clim_circle(c(s[[v1]], s[[v2]]), input$theta_clim,
                                        c(scl1$mean, scl2$mean), c(scl1$sd, scl2$sd)) %>%
                                  mutate(r = "b"),
                            clim_circle(c(s[[v1]], s[[v2]]), input$theta_clim * 3,
                                        c(scl1$mean, scl2$mean), c(scl1$sd, scl2$sd)) %>%
                                  mutate(r = "a"))

            dom <- if(!response_is_binary()) c(0, max(d$value, na.rm = TRUE)) else 0:1

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
                  theme(legend.position = "top",
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black", color = "black")) +
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

            if(nrow(filter(d, clim_dist <= input$theta_clim * 3)) == 0){
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
            d$value <- d$response
            d <- arrange(d, value) %>% mutate(panel = "Distance space")

            dom <- if(!response_is_binary()) c(0, max(d$value, na.rm = TRUE)) else 0:1

            ymax <- if(input$filter_clim) input$theta_clim * 3.3 else NA
            ggplot(d, aes(geog_dist, clim_dist, fill = value)) +
                  geom_vline(xintercept = input$theta_geog, linetype = "longdash", color = "#4c32a8", alpha = .5) +
                  geom_vline(xintercept = input$theta_geog * 3, color = "#4c32a8", alpha = .5) +
                  geom_hline(yintercept = input$theta_clim, linetype = "longdash", color = "#4c32a8", alpha = .5) +
                  geom_hline(yintercept = input$theta_clim * 3, color = "#4c32a8", alpha = .5) +
                  scale_x_continuous(expand = c(0, 0), limits = c(0, input$theta_geog * 3.3)) +
                  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +

                  geom_point(data = filter(d, samp_method_cd == 2), color = "gray70", size = 1.5) +
                  geom_point(data = filter(d, samp_method_cd == 1), shape = 21, color = "black", size = 3) +

                  geom_point(data = tibble::tibble(geog_dist = 0, clim_dist = 0,
                                                   shape_label = paste0("selected site,\n",
                                                                        eras()[2], " climate")),
                             aes(x = geog_dist, y = clim_dist, shape = shape_label),
                             color = "#4c32a8", fill = "#4c32a8",
                             size = 10, alpha = .5, inherit.aes = FALSE) +

                  marker_fill_scale(limits = dom) +
                  scale_shape_manual(values = 22) +
                  facet_wrap(~ panel) +
                  theme(legend.position = "top",
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black", color = "black")) +
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
                  filter(clim_dist <= input$theta_clim * 3)

            # Plot-level analog set: one row per (analog plot × analog direction),
            # joined to plot_occ for plot-level attributes (forest_type, forest_group).
            plot_analogs <- plot_occ()[d$analog_index, ] %>%
                  select(x, y, forest_type, forest_group) %>%
                  bind_cols(select(d, clim_dist, geog_dist, type)) %>%
                  mutate(weight = dnorm(clim_dist, sd = input$theta_clim) *
                               dnorm(geog_dist, sd = input$theta_geog))

            # SPECIES MODE: per-species decomposition using tree-level data.
            # Shows the selected response variable (presence / BA / BA proportion)
            # broken down by species across analog plots.
            if (input$track == "Species") {

                  analogs <- plot_analogs %>%
                        select(x, y, weight, type) %>%
                        left_join(fia_plots_species_all(), by = join_by(x, y),
                                  relationship = "many-to-many")

                  stat <- switch(input$summary_stat,
                                 "weighted mean" = weighted.mean,
                                 "weighted sum" = function(x, w) sum(x * w, na.rm = T))

                  a <- analogs %>%
                        select(x, y, weight, species, type) %>%
                        tidyr::expand(tidyr::nesting(x, y, weight, type), species) %>%
                        left_join(analogs, by = join_by(x, y, weight, species, type)) %>%
                        mutate(basal_area = ifelse(is.na(basal_area), 0, basal_area)) %>%
                        group_by(x, y, type) %>%
                        mutate(baplot = sum(basal_area),
                               baprop = basal_area / sum(basal_area),
                               baprop = ifelse(is.finite(baprop), baprop, 0)) %>%
                        group_by(species, type) %>%
                        summarize(
                              pres = stat(basal_area > 0, weight),
                              batot = stat(basal_area, weight),
                              baprop = stat(baprop, weight),
                              .groups = "drop"
                        ) %>%
                        filter(!is.na(species))

                  v <- switch(input$fia_var,
                              "presence probability"  = "pres",
                              "proportion basal area" = "baprop",
                              "total basal area"      = "batot")
                  a$value <- a[[v]]
                  a$category <- a$species
                  lab <- paste("historical", input$summary_stat, input$fia_var)

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
                  theme(legend.position = "top",
                        strip.text = element_text(color = "white"),
                        strip.background = element_rect(fill = "black", color = "black")) +
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
                                   scatter = NULL)


      eval_var_set <- function(clim_vars, vary_bandwidths = TRUE){
            # Extract historical climate from the raw raster at each pool plot,
            # matching pool().
            clim_r <- rast("data/clim_historical.tif")[[clim_vars]]

            po <- plot_occ()
            xy <- po[, c("x", "y")]
            clim_hst <- extract(clim_r, xy, ID = FALSE, method = "bilinear")
            d <- cbind(xy, clim_hst)
            values <- response_values()

            set.seed(123)
            train <- sample(c(TRUE, FALSE), nrow(xy), replace = T)

            aicc <- function(observed, predicted, k) {
                  # note: k should be n_clim_vars + n_bandwidth_params
                  n <- length(observed)
                  sse <- sum((observed - predicted)^2)
                  aic <- n * log(sse / n) + 2 * k # Basic AIC
                  aicc <- aic + (2 * k^2 + 2 * k) / (n - k - 1) # Corrected AICc
                  return(aicc)
            }

            eval <- function(theta_clim, theta_geog, x, index, values_train, values_test, summarize = TRUE){
                  aim <- suppressWarnings(analogs::analog_impact(
                        x = x,
                        pool = index,
                        max_clim = theta_clim * 3,
                        max_geog = theta_geog * 3,
                        stat = "weighted_mean",
                        kernel = "gaussian_joint",
                        theta = c(theta_clim, theta_geog),
                        y = values_train)) %>%
                        rename(pred = weighted_mean) %>%
                        mutate(obs = values_test) %>%
                        na.omit()
                  if(summarize){
                        aim <- aim  %>%
                              summarize(rsq = cor(pred, obs)^2,
                                        rmse = sqrt(mean((obs - pred)^2)),
                                        daicc = aicc(obs, pred, length(clim_vars) + 2),
                                        auc = auc(obs, pred)) %>%
                              mutate(theta_clim = theta_clim,
                                     theta_geog = theta_geog)
                  }
                  return(aim)
            }

            index <- build_analog_index(d[train,], coord_type = "lonlat", index_res = 8)

            if(vary_bandwidths){
                  mult <- c(.25, .5, 1, 2, 4)
                  clim_mult <- if("clim" %in% input$eval_params) mult else 1
                  geog_mult <- if("geog" %in% input$eval_params) mult else 1
                  params <- tidyr::expand_grid(theta_clim = input$theta_clim * clim_mult,
                                               theta_geog = input$theta_geog * geog_mult)

                  e <- purrr::pmap(params, eval,
                                   x = d[!train,], index = index,
                                   values_train = values[train], values_test = values[!train]) %>%
                        purrr::list_rbind() %>%
                        mutate(daicc = daicc - min(daicc)) %>% # convert AIC to delta AIC
                        tidyr::gather(stat, value, rsq, rmse, daicc, auc) %>%
                        mutate(clim_vars = paste(clim_vars, collapse = "_"))
            }else{
                  e <- eval(input$theta_clim, input$theta_geog,
                            x = d[!train,], index = index,
                            values_train = values[train], values_test = values[!train],
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
      })

      observeEvent(input$optimal_params, {
            req(model_eval$comparison)

            metric <- switch(input$eval_stat,
                             "RMSE" = "rmse",
                             "R-squared" = "rsq",
                             "delta AIC" = "daicc",
                             "AUC" = "auc")

            opt <- model_eval$comparison %>% ungroup() %>%
                  filter(stat == metric) %>%
                  mutate(value = ifelse(stat %in% c("rsq", "auc"), value, -value)) %>%
                  filter(value == max(value))

            updateSliderInput(session, "theta_geog", value = opt$theta_geog)
            updateSliderInput(session, "theta_clim", value = opt$theta_clim)
            updateSelectizeInput(session, "clim_vars", selected = unlist(strsplit(opt$clim_vars, "_")))
      })

      output$eval_heatmap <- renderPlot({
            req(model_eval$comparison)

            metric <- switch(input$eval_stat,
                             "RMSE" = "rmse",
                             "R-squared" = "rsq",
                             "delta AIC" = "daicc",
                             "AUC" = "auc")

            opt <- if(metric %in% c("rsq", "auc")) 1 else -1

            e <- model_eval$comparison %>% ungroup() %>%
                  filter(stat == metric) %>%
                  mutate(label = signif(value, 2),
                         optimal = value * opt == max(value * opt),
                         text_color = case_when(optimal ~ "red",
                                                scales::rescale(value) < .5 ~ "white",
                                                TRUE ~ "black"),
                         clim_vars = gsub("_", " + ", clim_vars),
                         n_vars = stringr::str_count(clim_vars, "\\+") + 1)

            ggplot(e, aes(theta_geog, theta_clim, fill = value, label = label)) +
                  facet_wrap(~clim_vars) +
                  geom_tile() +
                  geom_text(aes(color = I(text_color), fontface = "bold")) +
                  scale_x_log10(breaks = unique(e$theta_geog)) +
                  scale_y_log10(breaks = unique(e$theta_clim)) +
                  scale_fill_viridis_c() +
                  guides(fill = guide_colorbar(barwidth = 12)) +
                  theme(legend.position = "bottom") +
                  labs(y = "climate bandwidth (z)",
                       x = "geographic bandwidth (km)",
                       fill = paste0(input$eval_stat, " (", input$fia_var, ")"))
      })


      eval_scatter_data <- reactive({
            eval_var_set(input$clim_vars, vary_bandwidths = FALSE)
      })

      output$eval_scatter <- renderPlot({
            e <- eval_scatter_data()
            req(e)

            ggplot(e, aes(obs, pred)) +
                  geom_point() +
                  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "dodgerblue") +
                  geom_smooth(method = lm, se = F, color = "red", formula = 'y ~ x') +
                  annotate(geom = "point", color = "red", size = 4,
                           x = mean(e$obs), y = mean(e$pred)) +
                  scale_x_sqrt() +
                  scale_y_sqrt() +
                  coord_fixed() +
                  theme_bw() +
                  labs(x = paste0("observed ", input$fia_var),
                       y = paste0("predicted ", input$fia_var))
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
