# The settings record: a text file for citing results in reports. It pairs
# the link that reopens them (helpers/url_state.R) with what a link can't
# carry -- the settings in plain language, and the versions of the app and
# data that produced them.

# Versions, looked up once at startup (app.R)
app_versions <- function(dir = ".") {
      list(app = app_commit(dir),
           climate = formals(ensure_climate_data)$tag,
           boundaries = tryCatch(sf::st_read(file.path(dir, "data/boundaries.gpkg"), layer = "about",
                                             quiet = TRUE)$version,
                                 error = function(e) "unknown"),
           fia = substr(unname(tools::md5sum(file.path(dir, "data/fia.csv"))), 1, 12))
}

# The app's git commit. With git available, it also says whether the checkout
# has uncommitted changes to tracked files. Without git (a server may not have
# it, or may refuse a checkout owned by another user), the commit is read from
# .git directly and changes can't be checked, which the record says.
app_commit <- function(dir = ".") {
      git <- Sys.which("git")
      run <- function(...) tryCatch(suppressWarnings(system2(git, c("-C", dir, ...), stdout = TRUE, stderr = FALSE)),
                                    error = function(e) character(0))
      if (nzchar(git)) {
            sha <- run("rev-parse", "--short", "HEAD")
            if (length(sha) == 1 && is.null(attr(sha, "status")) && grepl("^[0-9a-f]{7,}$", sha)) {
                  changes <- run("status", "--porcelain", "--untracked-files=no")
                  return(paste0(sha, if (length(changes)) " (with uncommitted changes)"))
            }
      }
      head <- file.path(dir, ".git", "HEAD")
      if (!file.exists(head)) return("unknown")
      ref <- readLines(head, 1, warn = FALSE)
      sha <- if (startsWith(ref, "ref: ")) {
            ref <- sub("^ref: ", "", ref)
            loose <- file.path(dir, ".git", ref)
            if (file.exists(loose)) {
                  readLines(loose, 1, warn = FALSE)
            } else {
                  packed <- file.path(dir, ".git", "packed-refs")
                  hit <- if (file.exists(packed)) grep(paste0(" ", ref, "$"), readLines(packed, warn = FALSE), value = TRUE)
                  if (length(hit)) sub(" .*$", "", hit[1]) else NA_character_
            }
      } else ref
      if (is.na(sha) || !grepl("^[0-9a-f]{40}$", sha)) return("unknown")
      paste0(substr(sha, 1, 7), " (uncommitted changes not checked)")
}

# State (helpers/url_state.R) -> the record's lines. Lists only the settings
# that apply, labeled as the app's menus label them.
provenance_text <- function(s, ctx, link, versions, created = Sys.time()) {
      label <- function(choices, value) {
            l <- names(choices)[match(value, choices)]
            ifelse(is.na(l) | !nzchar(l), value, l)
      }
      num <- function(x) url_num(x)
      line <- function(name, value) sprintf("  %-24s %s", paste0(name, ":"), value)
      multi <- url_multi(s)

      entity <- switch(s$track,
                       "Species" = if (multi) "all species" else s$species,
                       "Forest type" = if (multi) "all forest types" else label(ctx$forest_types, s$forest_type),
                       "Forest group" = if (multi) "all forest groups" else label(ctx$forest_groups, s$forest_group),
                       "All trees" = "all trees")
      timeframe <- switch(s$time, historical = "historical", future = "future",
                          delta = "change from historical to future")
      floor_spec <- support_floor_spec(multi, s$fia_var, s$track)

      c("TreeShift settings record",
        paste("Created:", format(created, "%Y-%m-%d %H:%M %Z", tz = "UTC")),
        "",
        "Link to these results:",
        paste0("  ", link),
        "",
        "Model",
        line("Entity", paste0(label(track_choices, s$track), ": ", entity)),
        if (multi) line("Community measure", s$comm_measure),
        if (multi && s$comm_measure %in% categorical_measures) line("Categories shown", num(s$comm_top_n)),
        if (url_has_stage(s)) line("Life stage", label(stage_choices, s$stage)),
        line("Variable", s$metric),
        line("Timeframe", timeframe),
        if (url_has_floor(s)) line(floor_spec$label, num(s$support_floor)),
        "",
        "Study area",
        line("Grain size", paste(num(s$clim_agg), "km")),
        line("Longitude range", paste(num(s$bbox_x[1]), "to", num(s$bbox_x[2]))),
        line("Latitude range", paste(num(s$bbox_y[1]), "to", num(s$bbox_y[2]))),
        line("Focal polygon", if (nzchar(s$polygon_unit))
              sprintf("%s (%s)", ctx$unit_names[[s$polygon_unit]], s$polygon_unit) else "none"),
        line("Focal site", if (is.null(s$site)) "none" else
              sprintf("longitude %s, latitude %s", num(s$site[1]), num(s$site[2]))),
        "",
        "Climate",
        line("Climate variables", paste(toupper(s$clim_vars), collapse = ", ")),
        line("Emissions pathway", label(ctx$ssps, s$scenario_ssp)),
        line("Future time period", label(ctx$periods, s$scenario_period)),
        line("Climate models", paste(label(ctx$models, s$scenario_model), collapse = ", ")),
        "",
        "Interpolation",
        line("Geographic bandwidth", paste(num(s$theta_geog), "km")),
        line("Climate bandwidth", paste(num(s$theta_clim), "(z)")),
        if (url_cs_on(s)) c(
              "",
              "Connectivity",
              line("Sources", s$cs_sources),
              line("Targets", s$cs_targets),
              line("Horizontal conductance", s$cs_conductance_layer),
              line("Loss", s$cs_loss_layer),
              line("log10 conductance scale", num(s$cs_conductance_scale)),
              line("log10 loss scale", num(s$cs_loss_scale))),
        "",
        "Map layers",
        line("Raster variable", switch(s$stat,
              connectivity = paste("connectivity:", s$cs_stat),
              evidence = label(evidence_choices, s$evidence_stat),
              label(stat_choices_cs, s$stat))),
        line("Vector variable", s$vector),
        line("Analog direction", s$analog_direction),
        "",
        "Versions",
        line("TreeShift", paste("commit", versions$app)),
        line("Climate data", paste("release", versions$climate)),
        line("Boundary library", versions$boundaries),
        line("FIA plot data", paste("data/fia.csv, checksum", versions$fia)),
        "",
        "The link reopens these settings in TreeShift. Results can differ if the app or its data have changed since.")
}
