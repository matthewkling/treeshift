# Result cache -------------------------------------------------------------------
#
# Computed model surfaces (single-entity and community) and connectivity solves,
# kept on disk so that settings
# already run -- the defaults, a shared link, a step back -- reload instead of
# recomputing. Shared by every session and R process using the same directory.
#
# A cached result is reused only when its key matches, and the key covers
# everything the result is computed from:
#   - the settings the computation reads (ts_surfaces_spec(), ts_connectivity_spec());
#   - the content of the data files it reads;
#   - the versions of the packages it calls, and the installed files of the
#     non-CRAN ones (analogs, circuitscaper), which change without a version bump;
#   - for connectivity, the Julia version;
#   - the content of the code files that compute it.
# Anything that can't be fingerprinted disables caching for that result (the key
# is NULL) rather than risking a stale one. Any failure to read or write the cache
# falls back to computing.
#
# The default settings are put in the cache ahead of time by
# scripts/precompute_defaults.R.

#' The result cache
#'
#' @param dir Cache directory: `TREESHIFT_CACHE_DIR` if set, else `cache/` in the
#'   working directory (the app directory, for the app).
#' @return A `cachem` disk cache: least-recently-used entries are evicted past
#'   `max_size`, and entries expire after `max_age`.
ts_cache <- function(dir = Sys.getenv("TREESHIFT_CACHE_DIR", "cache"),
                     max_size = 5 * 1024^3, max_age = 30 * 24 * 3600) {
      cachem::cache_disk(dir, max_size = max_size, max_age = max_age, evict = "lru")
}

#' A cached result, or a freshly computed one
#'
#' @param cache `ts_cache()`, or NULL for no caching.
#' @param key A key from `ts_surfaces_key()` or `ts_connectivity_key()`, or NULL
#'   when the result can't be fingerprinted.
#' @param compute A function of no arguments that computes the result: a
#'   SpatRaster, or a list of them.
#' @param report Optional function, called with where the result came from:
#'   "hit" (the cache), "miss" (computed, then cached) or "off" (computed, no
#'   cache or key).
ts_cached <- function(cache, key, compute, report = function(source) NULL) {
      if (is.null(cache) || is.null(key)) {
            report("off")
            return(compute())
      }
      hit <- tryCatch(cache$get(key), error = function(e) cachem::key_missing())
      if (!cachem::is.key_missing(hit)) {
            value <- tryCatch(ts_unpack(hit), error = function(e) NULL)
            if (!is.null(value)) {
                  report("hit")
                  return(value)
            }
      }
      report("miss")
      value <- compute()
      tryCatch(cache$set(key, ts_pack(value)),
               error = function(e) warning("Couldn't save to the result cache: ",
                                           conditionMessage(e), call. = FALSE))
      value
}

# SpatRasters hold a pointer to C++ memory, so they are stored wrapped. wrap()
# stores the values themselves, not a reference to the raster's source file
# (checked: a wrapped file-backed raster still unwraps after its file is deleted).
# But it records the extent as text, to 15 significant digits, which moves the
# grid by ~1e-14 degrees; the exact extent is stored alongside and restored.
ts_pack <- function(x) {
      if (!inherits(x, "SpatRaster")) return(lapply(x, ts_pack))
      structure(list(raster = terra::wrap(x), extent = as.vector(terra::ext(x))),
                class = "ts_packed_raster")
}
ts_unpack <- function(x) {
      if (!inherits(x, "ts_packed_raster")) return(lapply(x, ts_unpack))
      r <- terra::unwrap(x$raster)
      terra::ext(r) <- terra::ext(x$extent)
      r
}


# Settings -----------------------------------------------------------------------

#' The settings a run's surfaces are computed from
#'
#' For a single entity, the entity is read only on its own track (ts_plots(),
#' ts_response()): the species on the Species track, the forest type or group
#' on theirs. A community run models every entity the track and variable give
#' (ts_multi_response()), so no single entity enters. The life stage and metric
#' enter only through `fia_var`. Values are normalized so a setting reached by
#' different routes (a slider, a link, a script) gives the same key; menu
#' orders are kept as given, since they set layer order.
#'
#' @param settings `ts_settings()`, or the app's model settings.
#' @param community Whether the surfaces are a community run's per-entity
#'   stack. Stated by the caller, not read from the settings, so a key can't
#'   pair one kind of result with the other.
ts_surfaces_spec <- function(settings, community = FALSE) {
      s <- settings
      list(track = as.character(s$track),
           entity = if (community) "community" else switch(s$track,
                           "Species"      = as.character(s$species),
                           "Forest type"  = as.integer(s$forest_type),
                           "Forest group" = as.integer(s$forest_group),
                           NULL),
           fia_var = as.character(s$fia_var),
           clim_vars = as.character(s$clim_vars),
           scenario = c(as.character(s$scenario_ssp), as.character(s$scenario_period),
                        as.character(s$scenario_model)),
           bbox = as.numeric(c(s$bbox_x, s$bbox_y)),
           clim_agg = as.numeric(s$clim_agg),
           theta = as.numeric(c(s$theta_geog, s$theta_clim)))
}

# The settings ts_surfaces_spec() reads, for callers that pass only these
surfaces_setting_ids <- c("track", "species", "forest_type", "forest_group", "fia_var",
                          "clim_vars", "scenario_ssp", "scenario_period", "scenario_model",
                          "bbox_x", "bbox_y", "clim_agg", "theta_geog", "theta_clim")

#' The settings a connectivity map is computed from
#'
#' The surfaces it is built on, the dispersal reach and land use weight
#' settings, and which map:
#' connectivity within the historical or future habitat, or the range shift
#' between them (server/connectivity.R).
#'
#' @param era "hist", "fut" or "shift".
ts_connectivity_spec <- function(settings, era) {
      list(surfaces = ts_surfaces_spec(settings),
           reach = as.character(settings$cs_reach),
           landuse = as.numeric(settings$cs_landuse),
           era = era)
}

connectivity_setting_ids <- c("cs_reach", "cs_landuse")


# Fingerprints -------------------------------------------------------------------

# Content hash of files; NULL if any is missing
file_hashes <- function(paths) {
      if (!all(file.exists(paths))) return(NULL)
      vapply(paths, function(p) digest::digest(file = p, algo = "xxhash64"), character(1),
             USE.NAMES = FALSE)
}

# Package versions, plus the installed files of the given packages
package_fingerprint <- function(packages, by_files = character(0)) {
      list(r = R.version.string,
           versions = vapply(packages, function(p) as.character(utils::packageVersion(p)),
                             character(1)),
           files = lapply(stats::setNames(nm = by_files), function(p) {
                 dir <- system.file(package = p)
                 file_hashes(sort(list.files(dir, recursive = TRUE, full.names = TRUE)))
           }))
}

# The code that computes the surfaces: the core (R/) and the server glue that
# assembles its inputs (clim(), pool() and aim()). Paths relative to the app
# directory, which scripts and the app both run from.
surfaces_code_files <- function() {
      c(sort(list.files("R", pattern = "^[^_].*[.]R$", full.names = TRUE)),
        "server/aim.R", "server/climate.R", "server/climate_models.R", "server/data.R")
}
community_code_files <- function() c(surfaces_code_files(), "server/community.R")
connectivity_code_files <- function() c(surfaces_code_files(), "server/connectivity.R")

# Packages whose code computes the surfaces
surfaces_packages <- c("analogs", "terra", "dplyr", "tidyr", "purrr", "readr", "tibble")

#' Cache key for a run's surfaces
#'
#' @param settings `ts_settings()`, or the app's model settings.
#' @param data `ts_data()`.
#' @param community See `ts_surfaces_spec()`.
#' @return A key string, or NULL if a file it depends on is missing.
ts_surfaces_key <- function(settings, data, community = FALSE) {
      inputs <- c(file.path(data$dir, c("fia.csv", "forest_type_codes.csv", "clim_scales.rds")),
                  data$clim_historical, data$conus,
                  ts_scenario_path(data, settings$scenario_ssp, settings$scenario_period,
                                    settings$scenario_model))
      code <- if (community) community_code_files() else surfaces_code_files()
      parts <- list(spec = ts_surfaces_spec(settings, community),
                    data = file_hashes(inputs),
                    code = file_hashes(code),
                    packages = package_fingerprint(surfaces_packages, by_files = "analogs"))
      if (is.null(parts$data) || is.null(parts$code)) return(NULL)
      paste0(if (community) "community-" else "surfaces-", rlang::hash(parts))
}

# The Julia version, read without starting Julia; NULL if it can't be found.
# The dipole methods use only Julia's standard library, which comes with it, and
# their own Julia code, which the package fingerprint covers.
julia_fingerprint <- function() {
      v <- tryCatch(circuitscaper:::julia_version(circuitscaper:::find_julia_bin()),
                    error = function(e) NULL)
      if (!length(v)) NULL else as.character(v)
}

#' Cache key for a connectivity map
#'
#' @param era See `ts_connectivity_spec()`.
#' @return A key string, or NULL if the surfaces, the human modification layer
#'   or the solver can't be fingerprinted.
ts_connectivity_key <- function(settings, data, era) {
      surfaces <- ts_surfaces_key(settings, data)
      julia <- julia_fingerprint()
      parts <- list(surfaces = surfaces,
                    spec = ts_connectivity_spec(settings, era),
                    ghm = file_hashes(file.path(data$dir, "ghm.tif")),
                    code = file_hashes(connectivity_code_files()),
                    packages = package_fingerprint(c("circuitscaper", "JuliaConnectoR", "terra"),
                                                   by_files = "circuitscaper"),
                    julia = julia)
      if (is.null(surfaces) || is.null(parts$ghm) || is.null(parts$code) || is.null(julia))
            return(NULL)
      paste0("connectivity-", rlang::hash(parts))
}
