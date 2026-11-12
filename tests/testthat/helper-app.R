# Shared by the tests. The app directory is two levels up from tests/testthat.
app_dir <- normalizePath(file.path(testthat::test_path(), "..", ".."))

# Evaluate app.R's global setup (everything before the UI) into the global
# environment, for tests of the helpers that need the app's data and choices.
load_app_globals <- function() {
      if (exists("url_ctx", envir = globalenv())) return(invisible())
      lines <- readLines(file.path(app_dir, "app.R"))
      end <- grep("^# UI ----", lines) - 1
      withr::with_dir(app_dir, suppressMessages(
            eval(parse(text = lines[seq_len(end)]), envir = globalenv())))
      invisible()
}

# A fresh session of a running app, optionally opened on a link's query
open_session <- function(url) {
      app <- shinytest2::AppDriver$new(url, width = 1400, height = 900,
                                       load_timeout = 600000, timeout = 600000)
      settle(app)
      app
}

# Wait until the app has been idle long enough for the debounced address-bar
# update to land
settle <- function(app) app$wait_for_idle(1000, timeout = 600000)

# Set inputs one at a time, letting each one's dependent menus update before
# the next (track -> stage -> metric, collection -> unit)
set_in_order <- function(app, ...) {
      vals <- list(...)
      for (id in names(vals)) {
            do.call(app$set_inputs, c(stats::setNames(vals[id], id), list(wait_ = FALSE)))
            settle(app)
      }
}

# Press Run model, applying the model settings set so far
apply_settings <- function(app) {
      app$click("update_model")
      settle(app)
}

# Browser Back or Forward. These reload the page when the settings differ,
# and shinytest2 can't wait on a reloaded page, so this waits for a new
# session (the marker set on the old page is gone) to connect and go idle.
# Returns whether the page reloaded.
history_nav <- function(app, how = "back") {
      app$run_js(sprintf("window.__before_nav = true; history.%s();", how))
      reloaded <- FALSE
      for (i in 1:300) {
            Sys.sleep(1)
            state <- tryCatch(app$get_js(
                  "window.__before_nav ? 'same' : (window.Shiny && Shiny.shinyapp && Shiny.shinyapp.isConnected() && !document.documentElement.classList.contains('shiny-busy')) ? 'ready' : 'loading'"),
                  error = function(e) "loading")
            if (state == "ready") { reloaded <- TRUE; break }
            if (state == "same" && i >= 5) break
      }
      Sys.sleep(2)
      reloaded
}

# An input's value, read in the page (works after a reload, unlike get_value)
page_input <- function(app, id) app$get_js(sprintf("Shiny.shinyapp.$inputValues['%s']", id))

click_map <- function(app, lon, lat) {
      app$run_js(sprintf("Shiny.setInputValue('map_click', {lng: %s, lat: %s, '.nonce': Math.random()});", lon, lat))
      settle(app)
}

address_bar <- function(app) app$get_js("window.location.search")

notices <- function(app) {
      app$get_js("Array.from(document.querySelectorAll('.shiny-notification-content')).map(n => n.innerText).join('\\n')")
}

# The displayed raster, downloaded. Download links are only filled in once the
# Export panel has been opened. The link's address is resolved by the browser:
# AppDriver$get_download() appends it to the session URL, which breaks when
# that URL has a query string (it fetches the page instead).
download_raster <- function(app) {
      app$run_js("var b = Array.from(document.querySelectorAll('.accordion-button')).find(b => b.textContent.trim() === 'Export'); if (b.classList.contains('collapsed')) b.click();")
      ## The link's address can arrive after the app goes idle
      app$wait_for_js("document.getElementById('download_raster').href.indexOf('/download/') >= 0",
                      timeout = 600000)
      href <- app$get_js("document.getElementById('download_raster').href")
      f <- tempfile(fileext = ".tif")
      utils::download.file(href, f, quiet = TRUE, mode = "wb")
      stopifnot(identical(readBin(f, "raw", 2), charToRaw("II")))   # a GeoTIFF, not a page
      f
}

raster_md5 <- function(app) {
      f <- download_raster(app)
      on.exit(unlink(f))
      unname(tools::md5sum(f))
}

# The displayed raster's values
raster_values <- function(app) {
      f <- download_raster(app)
      on.exit(unlink(f))
      terra::values(terra::rast(f), mat = FALSE)
}

# The input each link key sets
key_inputs <- c(e = "track", es = "entity_scope", sp = "species", ft = "forest_type",
                fg = "forest_group", st = "stage", m = "metric", t = "time",
                cm = "comm_measure", n = "comm_top_n", sf = "support_floor", g = "clim_agg",
                bx = "bbox_x", by = "bbox_y", poly = "polygon_unit", cv = "clim_vars",
                ssp = "scenario_ssp", era = "scenario_period", gcm = "scenario_model",
                tg = "theta_geog", tc = "theta_clim", cs = "cs_enable", csrc = "cs_sources",
                ctgt = "cs_targets", ccl = "cs_conductance_layer", cll = "cs_loss_layer",
                ccs = "cs_conductance_scale", cls = "cs_loss_scale", r = "stat",
                rev = "evidence_stat", rcs = "cs_stat", pts = "vector", an = "analog_direction")

# The inputs a link sets: those its keys name, plus the variable they resolve
# to and the map layers, which links leave out when at their defaults
inputs_set_by <- function(query) {
      keys <- sub("=.*$", "", strsplit(sub("^\\?", "", query), "&", fixed = TRUE)[[1]])
      unique(c(unname(key_inputs[intersect(keys, names(key_inputs))]),
               "fia_var", "stat", "vector", "analog_direction"))
}
