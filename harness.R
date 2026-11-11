# ---- treeshift_harness.R ------------------------------------------------
#
# Scripted control of the TreeShift Shiny app via shinytest2::AppDriver.
#
#   app <- ts_start("path/or/url")
#   ts_render(app, track = "Species", species = "Pinus ponderosa",
#             metric = "occurrence", stage = "any", focal_domain = TRUE,
#             stat = "FIA variable", time = "historical",
#             download = "raster", outdir = "out", tag = "pipo_hist")
#
# Every app control is a named argument of ts_render(); buttons are logical
# flags. Arguments may be given in any order -- the harness applies them in
# the order the app's observers require. ts_batch() loops ts_render() over a
# data frame whose columns are ts_render() argument names.
#
# Requires no changes to app.R.

library(shinytest2)


# ==========================================================================
#  App schema
# ==========================================================================

# Flush order. Controls sharing a phase are set together; lower phases first.
# The ordering matters because several observers repopulate downstream
# controls and reset them to their first choice:
#
#   track  -> updateSelectInput("metric", ..., selected = mc[1])
#   metric -> updateSelectInput("stage",  ..., selected = sc[1])
#
# so `stage` set alongside `metric` would be silently overwritten. The snap
# buttons need track and species, and must precede any explicit bbox.
#
# `fia_var` is deliberately absent: the app derives it from
# (track, metric, stage) and it must not be set directly.
TS_PHASES <- c(
      track                = 1,

      species              = 2,
      forest_type          = 2,
      forest_group         = 2,

      focal_domain         = 3,   # button: snap to focal entity's range

      metric               = 4,

      stage                = 5,

      bbox_x               = 6,
      bbox_y               = 6,
      clim_agg             = 6,
      clim_vars            = 6,
      scenario_ssp         = 6,
      scenario_period      = 6,
      scenario_model       = 6,
      theta_geog           = 6,
      theta_clim           = 6,
      cs_sources           = 6,
      cs_targets           = 6,
      cs_conductance_layer = 6,
      cs_loss_layer        = 6,
      cs_conductance_scale = 6,
      cs_loss_scale        = 6,
      stat                 = 6,
      time                 = 6,
      vector               = 6,
      vector_scale         = 6,
      rast_transform       = 6,
      rast_palette         = 6,

      cs_stat              = 7,   # conditionalPanel on stat == "connectivity"

      window_domain        = 8    # button: snap to current map view
)

TS_BUTTONS <- c("focal_domain", "window_domain")

TS_DOWNLOADS    <- c(raster = "download_raster",
                     points = "download_points",
                     metadata = "download_metadata")
TS_DOWNLOAD_EXT <- c(raster = "tif", points = "csv", metadata = "txt")

# Accordion panel -> a selector present only once that panel's outputs have
# been resumed. The app's accordion is open = "Info", multiple = FALSE, and
# no output sets suspendWhenHidden = FALSE, so these outputs never receive a
# value until their panel has been visible once.
TS_RESUME <- list(
      "Spatial grid" = "#focal_domain",
      "Export"       = "a#download_raster:not([href=''])"
)


# ==========================================================================
#  Session
# ==========================================================================

ts_start <- function(app = ".",
                     load_timeout = 15 * 60 * 1000,
                     timeout = 60 * 60 * 1000,
                     resume = TRUE,
                     ...) {
      drv <- AppDriver$new(app,
                           name = "treeshift-batch",
                           load_timeout = load_timeout,
                           timeout = timeout,
                           ...)
      if (resume) ts_resume(drv)
      drv
}

ts_js_str <- function(x) paste0("'", gsub("'", "\\\\'", x), "'")

ts_exists <- function(app, selector) {
      isTRUE(app$get_js(sprintf("!!document.querySelector(%s)", ts_js_str(selector))))
}

# Open each accordion panel once so Shiny resumes its outputs. Resumption is
# sticky -- a download link keeps its href and a rendered uiOutput keeps its
# markup after the panel collapses again -- so it does not matter that
# accordion(multiple = FALSE) closes each panel as the next one opens.
#
# One panel at a time: Bootstrap drops a collapse toggle issued while another
# collapse is mid-transition. Only the panels in TS_RESUME are touched;
# opening the right sidebar's analysis panels would resume their plot outputs
# and start real computation.
ts_resume <- function(app, panels = names(TS_RESUME),
                      timeout = 30, interval = 0.25) {
      for (p in panels) {
            sel <- TS_RESUME[[p]]
            if (ts_exists(app, sel)) next

            found <- app$get_js(sprintf("
                  (() => {
                        const b = Array.from(document.querySelectorAll('.accordion-button'))
                              .find(x => x.textContent.trim().startsWith(%s));
                        if (!b) return false;
                        if (b.classList.contains('collapsed')) b.click();
                        return true;
                  })()
            ", ts_js_str(p)))

            if (!isTRUE(found)) {
                  stop("no accordion panel whose header starts with '", p, "'",
                       call. = FALSE)
            }

            deadline <- Sys.time() + timeout
            while (!ts_exists(app, sel) && Sys.time() < deadline) Sys.sleep(interval)
            if (!ts_exists(app, sel)) {
                  stop("opened panel '", p, "' but '", sel, "' never appeared",
                       call. = FALSE)
            }
      }
      invisible(app)
}


# ==========================================================================
#  Discovery
# ==========================================================================

ts_input_ids <- function(app) sort(names(app$get_values()$input))

ts_choices <- function(app, input_id) {
      unlist(app$get_js(sprintf(
            "Array.from(document.querySelectorAll('#%s option')).map(e => e.value)",
            input_id)))
}

# What the app currently thinks it is modelling, including the hidden
# canonical variable.
ts_state <- function(app) {
      v <- app$get_values()$input
      v[intersect(c("track", "species", "metric", "stage", "fia_var",
                    "stat", "time", "vector", "theta_clim", "theta_geog",
                    "bbox_x", "bbox_y", "clim_agg", "clim_vars",
                    "scenario_ssp", "scenario_period", "scenario_model"),
                  names(v))]
}


# ==========================================================================
#  Internals
# ==========================================================================

ts_idle <- function(app, duration = 2000, timeout = NULL) {
      if (is.null(timeout)) app$wait_for_idle(duration = duration)
      else app$wait_for_idle(duration = duration, timeout = timeout)
      invisible(app)
}

# Named values -> ordered list of stages, each a named list of inputs or a
# click. Buttons flagged FALSE are dropped.
ts_build_stages <- function(vals) {
      vals <- vals[!vapply(vals, is.null, logical(1))]

      is_btn <- names(vals) %in% TS_BUTTONS
      keep <- rep(TRUE, length(vals))
      keep[is_btn] <- vapply(vals[is_btn], isTRUE, logical(1))
      vals <- vals[keep]
      if (!length(vals)) return(list())

      phase <- TS_PHASES[names(vals)]
      stages <- list()

      for (p in sort(unique(phase))) {
            this <- vals[phase == p]
            btn <- names(this)[names(this) %in% TS_BUTTONS]
            inp <- this[!(names(this) %in% TS_BUTTONS)]
            if (length(inp)) stages <- c(stages, list(inp))
            for (b in btn) {
                  stages <- c(stages, list(structure(list(input = b),
                                                     class = "ts_click")))
            }
      }
      stages
}

ts_apply <- function(app, stages, idle_ms = 2000, timeout = NULL) {
      for (stage in stages) {
            if (inherits(stage, "ts_click")) {
                  if (!ts_exists(app, paste0("#", stage$input))) {
                        stop("button '", stage$input, "' is not in the DOM; ",
                             "its output is still suspended (see ts_resume)",
                             call. = FALSE)
                  }
                  app$click(input = stage$input)
            } else {
                  do.call(app$set_inputs, c(stage, list(wait_ = FALSE)))
            }
            ts_idle(app, duration = idle_ms, timeout = timeout)
      }
      invisible(app)
}

# Read inputs back and compare against what was asked for. Silent divergence
# is the dominant failure mode: a select given a value that isn't among its
# options resolves to "" without error, an observer can reset a control, and
# shiny::validate() cancels downstream reactives without raising one. All
# produce plausible output for the wrong parameters.
ts_verify <- function(app, wanted, skip = character(0)) {
      got <- app$get_values()$input
      bad <- character(0)

      for (nm in setdiff(names(wanted), c(TS_BUTTONS, skip))) {
            w <- wanted[[nm]]
            if (is.null(w)) next
            g <- got[[nm]]

            ok <- if (is.null(g)) FALSE
            else if (is.numeric(w) && is.numeric(g))
                  length(w) == length(g) && all(abs(w - g) < 1e-8)
            else identical(sort(as.character(w)), sort(as.character(g)))

            if (!ok) {
                  bad <- c(bad, sprintf("%s: asked %s, got %s", nm,
                                        paste(w, collapse = "/"),
                                        if (is.null(g)) "NULL" else paste(g, collapse = "/")))
            }
      }

      if (length(bad)) {
            stop("inputs did not take:\n  ", paste(bad, collapse = "\n  "),
                 "\n  (check ts_choices() for the legal values)", call. = FALSE)
      }
      invisible(TRUE)
}

# A Shiny error paints red text in the UI without failing the driver.
# validate() messages carry the same class, so an unavailable scenario
# combination is caught here too.
ts_errors <- function(app) {
      html <- tryCatch(app$get_html(".shiny-output-error", outer_html = TRUE),
                       error = function(e) NULL)
      if (length(html) == 0) return(character(0))
      trimws(gsub("\\s+", " ", gsub("<[^>]*>", " ", html)))
}


# ==========================================================================
#  ts_render
# ==========================================================================

# `app` is a live AppDriver, or a path/URL (a driver is then started and
# stopped around the call -- convenient for one-offs, wasteful in a loop).
#
# `download` names any of "raster", "points", "metadata". Files are written
# to file.path(outdir, paste0(tag, "_", name, ".", ext)). With download = NULL
# the app is simply left in the requested state.
ts_render <- function(app,

                      # FIA variable
                      track = NULL, species = NULL,
                      forest_type = NULL, forest_group = NULL,
                      metric = NULL, stage = NULL,

                      # spatial grid
                      focal_domain = FALSE, window_domain = FALSE,
                      bbox_x = NULL, bbox_y = NULL, clim_agg = NULL,

                      # climate
                      clim_vars = NULL, scenario_ssp = NULL,
                      scenario_period = NULL, scenario_model = NULL,

                      # analog settings
                      theta_geog = NULL, theta_clim = NULL,

                      # connectivity
                      cs_sources = NULL, cs_targets = NULL,
                      cs_conductance_layer = NULL, cs_loss_layer = NULL,
                      cs_conductance_scale = NULL, cs_loss_scale = NULL,

                      # map layer / cartography
                      stat = NULL, cs_stat = NULL, time = NULL, vector = NULL,
                      vector_scale = NULL, rast_transform = NULL,
                      rast_palette = NULL,

                      # export
                      download = NULL, outdir = ".", tag = "treeshift_result",

                      verify = TRUE, idle_ms = 2000, timeout = NULL) {

      if (is.character(app)) {
            drv <- ts_start(app)
            on.exit(try(drv$stop(), silent = TRUE), add = TRUE)
            app <- drv
      }

      stages <- ts_build_stages(mget(names(TS_PHASES)))

      # A click stage must exist for exactly the buttons flagged TRUE.
      clicked <- vapply(Filter(function(s) inherits(s, "ts_click"), stages),
                        function(s) s$input, character(1))
      asked <- TS_BUTTONS[vapply(TS_BUTTONS, function(b) isTRUE(get(b)), logical(1))]
      if (!identical(sort(clicked), sort(asked))) {
            stop("stage construction is wrong: would click [",
                 paste(clicked, collapse = ", "), "] but was asked for [",
                 paste(asked, collapse = ", "), "]", call. = FALSE)
      }

      # The export filename rides along in the first input stage; it appears
      # inside the metadata file.
      nm <- list(download_filename = tag)
      i <- which(!vapply(stages, inherits, logical(1), "ts_click"))
      if (length(i)) stages[[i[1]]] <- c(nm, stages[[i[1]]])
      else stages <- c(list(nm), stages)

      ts_apply(app, stages, idle_ms = idle_ms, timeout = timeout)

      if (verify) {
            # A snap button rewrites the bbox by design, so don't hold the
            # sliders to values the caller didn't set.
            skip <- if (isTRUE(focal_domain) || isTRUE(window_domain))
                  c("bbox_x", "bbox_y")[c(is.null(bbox_x), is.null(bbox_y))]
            else character(0)
            ts_verify(app, mget(names(TS_PHASES)), skip = skip)
      }

      errs <- ts_errors(app)
      if (length(errs)) {
            stop("app reported output error(s): ", paste(errs, collapse = " | "),
                 call. = FALSE)
      }

      if (is.null(download)) return(invisible(character(0)))

      bad <- setdiff(download, names(TS_DOWNLOADS))
      if (length(bad)) {
            stop("unknown download(s): ", paste(bad, collapse = ", "),
                 "; available: ", paste(names(TS_DOWNLOADS), collapse = ", "),
                 call. = FALSE)
      }

      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

      # Attempt all of them: if the raster handler fails, the metadata file is
      # often what explains why.
      res <- lapply(download, function(d) {
            dest <- file.path(outdir, sprintf("%s_%s.%s", tag, d, TS_DOWNLOAD_EXT[[d]]))
            tryCatch({
                  app$get_download(TS_DOWNLOADS[[d]], filename = dest)
                  list(path = dest, error = NA_character_)
            }, error = function(e) list(path = NA_character_,
                                        error = conditionMessage(e)))
      })

      failed <- Filter(function(x) !is.na(x$error), res)
      if (length(failed)) {
            stop("download failed: ",
                 paste(vapply(failed, function(x) x$error, character(1)),
                       collapse = "; "), call. = FALSE)
      }

      stats::setNames(vapply(res, function(x) x$path, character(1)), download)
}


# ==========================================================================
#  ts_batch
# ==========================================================================

ts_default_tag <- function(row) {
      v <- vapply(row, function(x) paste(as.character(x), collapse = "-"), character(1))
      paste(gsub("^-|-$", "", gsub("[^A-Za-z0-9.]+", "-", v)), collapse = "_")
}

ts_append_manifest <- function(path, row) {
      ex <- file.exists(path)
      utils::write.table(row, path, sep = ",", row.names = FALSE,
                         col.names = !ex, append = ex, qmethod = "double")
}

# `grid` columns and `constants` names must be ts_render() argument names.
# Constants hold anything fixed across the batch, so invariants don't inflate
# the cartesian product.
ts_batch <- function(grid,
                     app = ".",
                     outdir = "batch_output",
                     download = c("raster", "metadata"),
                     constants = list(),
                     tag_fun = ts_default_tag,
                     restart_every = 10,
                     skip_existing = TRUE,
                     idle_ms = 2000,
                     load_timeout = 15 * 60 * 1000,
                     timeout = 60 * 60 * 1000,
                     manifest = file.path(outdir, "manifest.csv"),
                     logdir = file.path(outdir, "logs")) {

      bad <- setdiff(c(names(grid), names(constants)), names(TS_PHASES))
      if (length(bad)) {
            stop("not app controls: ", paste(bad, collapse = ", "),
                 if ("fia_var" %in% bad)
                       "\n  (fia_var is derived from track/metric/stage; set those instead)",
                 call. = FALSE)
      }

      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

      driver <- NULL
      since_restart <- Inf
      on.exit(if (!is.null(driver)) try(driver$stop(), silent = TRUE), add = TRUE)

      for (i in seq_len(nrow(grid))) {

            row <- lapply(as.list(grid[i, , drop = FALSE]),
                          function(x) if (is.factor(x)) as.character(x) else x)
            tag <- tag_fun(row)

            expected <- file.path(outdir, sprintf("%s_%s.%s", tag, download,
                                                  TS_DOWNLOAD_EXT[download]))
            if (skip_existing && all(file.exists(expected))) {
                  message(sprintf("[%d/%d] %s - skipped (exists)", i, nrow(grid), tag))
                  next
            }

            if (since_restart >= restart_every) {
                  if (!is.null(driver)) try(driver$stop(), silent = TRUE)
                  message("starting app driver...")
                  driver <- ts_start(app, load_timeout = load_timeout,
                                     timeout = timeout, resume = TRUE)
                  since_restart <- 0
            }

            message(sprintf("[%d/%d] %s", i, nrow(grid), tag))
            t0 <- Sys.time()

            res <- tryCatch({
                  paths <- do.call(ts_render,
                                   c(list(app = driver, download = download,
                                          outdir = outdir, tag = tag,
                                          idle_ms = idle_ms),
                                     constants, row))
                  list(status = "ok", message = "",
                       files = paste(basename(paths), collapse = ";"))
            }, error = function(e) {
                  try(writeLines(utils::capture.output(print(driver$get_logs())),
                                 file.path(logdir, paste0(tag, ".log"))),
                      silent = TRUE)
                  list(status = "error", message = conditionMessage(e), files = "")
            })

            if (res$status == "error") {
                  since_restart <- Inf   # session may be wedged; start fresh
                  warning(sprintf("%s failed: %s", tag, res$message), call. = FALSE)
            } else {
                  since_restart <- since_restart + 1
            }

            ts_append_manifest(manifest,
                               cbind(as.data.frame(row, stringsAsFactors = FALSE),
                                     data.frame(tag = tag,
                                                status = res$status,
                                                elapsed_s = round(as.numeric(
                                                      difftime(Sys.time(), t0, units = "secs")), 1),
                                                files = res$files,
                                                message = res$message,
                                                stringsAsFactors = FALSE)))
      }

      invisible(utils::read.csv(manifest, stringsAsFactors = FALSE))
}
