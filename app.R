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

# Helper functions and constants live in helpers/, the UI's parts in ui/, and the
# server's in server/. They are sourced explicitly into this file's environment.
# (Not from an R/ folder: Shiny sources that into a parent environment, where
# the functions could not see the data loaded below.) helpers/ and ui/ load in
# alphabetical order, so a value computed at top level in one of them may only
# use definitions from its own file.
for (.src in list.files(c("helpers", "ui"), pattern = "[.]R$", full.names = TRUE)) {
      source(.src, local = TRUE, encoding = "UTF-8", keep.source = TRUE)
}


# Fetch climate data if absent. Runs once at app process startup, before any
# rast("data/...") call below. On the production shiny-server deployment,
# prefer running this once as an explicit provisioning step (from an R
# console) so app startup stays fast; this in-app call then just confirms
# presence and returns immediately.
ensure_climate_data()

# The focal-polygon library (see helpers/boundaries.R). Small enough to hold in
# memory: one row per unit, EPSG:4326.
ensure_boundaries()
boundary_units <- st_read("data/boundaries.gpkg", layer = "units", quiet = TRUE)

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

theme_set(theme_bw(base_size = 16))

# House style for facet strips: white text on a black background. Defined
# once here and added to every faceted plot via `+ strip_style`. It is a
# partial theme() (only the two strip properties), so each plot keeps its
# own legend / axis theme settings and simply layers this on top.
strip_style <- theme(strip.text = element_text(color = "white"),
                     strip.background = element_rect(fill = "black", color = "black"))


# The valid choices that links are checked against (helpers/url_state.R)
url_ctx <- url_context()

# App and data versions, for the settings record (helpers/provenance.R)
app_info <- app_versions()


# UI ---------------------------------------------------
# A function of the request, so the page is built from the settings in the
# link it was opened with (helpers/url_state.R).
ui <- function(request) {
      s <- url_decode(request$QUERY_STRING, url_ctx)$state
      page_navbar(
            id = "main_nav",

            header = tagList(ui_head(),
                             add_busy_spinner(spin = "cube-grid", position = "bottom-right")),

            ## Monochrome header: very dark gray. The selected tab is white;
            ## the title and other tabs are gray, turning white on hover.
            ## (Flatly's own links turn green on hover and when selected.)
            theme = bs_theme(version = 5, bootswatch = "flatly",
                             "navbar-dark-color" = "rgba(255, 255, 255, 0.55)",
                             "navbar-dark-brand-color" = "rgba(255, 255, 255, 0.55)",
                             "navbar-dark-hover-color" = "#ffffff",
                             "navbar-dark-active-color" = "#ffffff"),
            navbar_options = navbar_options(bg = "#212529", theme = "dark"),

            fillable = "map",   # lets the map tab fill the browser height; About scrolls

            title = "TreeShift [beta]",

            nav_panel(
                  "Map",
                  value = "map",
                  layout_sidebar(
                        fillable = TRUE,
                        padding = 0,
                        border = FALSE,
                        border_radius = FALSE,

                        ## LEFT SIDEBAR ------------------------------------

                        sidebar = ui_left_sidebar(s),

                        card(
                              class = "h-100 w-100 rounded-0 overflow-hidden border-0",
                              layout_sidebar(
                                    fillable = TRUE,
                                    padding = 0,   # the map fills its area, edge to edge


                                    ## RIGHT SIDEBAR ---------------------------------------

                                    sidebar = ui_right_sidebar(s),

                                    ## MAP ---------------------------------------
                                    div(class = "h-100 position-relative",
                                        leafletOutput("map", height = "100%"),
                                        uiOutput("map_pending"))
                              )
                        )
                  )
            ),

            nav_panel("About", value = "about", ui_about())
      )
}



# SERVER ---------------------------------------------------

# The server is split by feature across server/*.R. Each file is sourced into
# this function's environment, in the order below, so the files share reactives
# exactly as if they were written inline. The order matters: observers run in
# the order they are created, and some files use objects defined in earlier ones.
server_files <- c("shared", "url_state", "model_settings", "data", "controls", "climate", "climate_models",
                  "aim", "connectivity", "map_base", "community", "map_layers",
                  "site_plots", "polygon_values", "evaluation", "climate_response", "downloads")

server <- function(input, output, session) {
      for (.src in file.path("server", paste0(server_files, ".R"))) {
            source(.src, local = TRUE, encoding = "UTF-8", keep.source = TRUE)
      }
}

shinyApp(ui, server)
