# Input data ---------------------------------------------------------------------
#
# Everything the model reads from data/, loaded once. The app builds one of these
# at startup (app.R); scripts and tests build their own. The climate rasters are
# large, so only their paths are held here; they are read per run.

#' Load the model's input data
#'
#' @param dir Directory holding fia.csv, forest_type_codes.csv,
#'   clim_historical.tif, clim_scales.rds and clim_future/.
#' @return A `ts_data` list: `fia` (one row per plot and species, with
#'   `forest_type` 0 for non-forest plots), `forest_types` (FORTYPCD lookup,
#'   including the non-forest code 0), `scenarios` (one row per future climate
#'   file), `clim_scales`, and the climate file paths.
ts_data <- function(dir = "data") {
      # Absolute, so the climate paths held below still resolve from another
      # working directory
      dir <- normalizePath(dir, mustWork = TRUE)
      fia <- readr::read_csv(file.path(dir, "fia.csv"), show_col_types = FALSE)
      # Non-forest plots have no FORTYPCD. Recoded to the sentinel 0, which the
      # lookup below names, so filters and joins treat them like any other type.
      fia$forest_type <- ifelse(is.na(fia$forest_type), 0L, fia$forest_type)

      forest_types <- readr::read_csv(file.path(dir, "forest_type_codes.csv"),
                                      show_col_types = FALSE) |>
            dplyr::bind_rows(tibble::tibble(forest_type = 0L,
                                            forest_type_name = "Non-forest",
                                            forest_group = 0L,
                                            forest_group_name = "Non-forest",
                                            forest_group_abbr = "nonfor"))

      # Future climate files are named <ssp>_<period>_<model>.tif
      files <- list.files(file.path(dir, "clim_future"), pattern = "\\.tif$")
      parts <- strsplit(tools::file_path_sans_ext(files), "_", fixed = TRUE)
      scenarios <- tibble::tibble(file   = files,
                                  ssp    = vapply(parts, `[`, character(1), 1),
                                  period = vapply(parts, `[`, character(1), 2),
                                  model  = vapply(parts, `[`, character(1), 3))

      structure(list(dir = dir,
                     fia = fia,
                     forest_types = forest_types,
                     scenarios = scenarios,
                     clim_scales = readRDS(file.path(dir, "clim_scales.rds")),
                     clim_historical = file.path(dir, "clim_historical.tif"),
                     clim_future_dir = file.path(dir, "clim_future")),
                class = "ts_data")
}

#' Forest type and group names from their FIA codes (NA for unknown codes)
ts_forest_type_name <- function(data, code) {
      ft <- data$forest_types
      ft$forest_type_name[match(code, ft$forest_type)]
}

ts_forest_group_name <- function(data, code) {
      g <- dplyr::distinct(data$forest_types, forest_group, forest_group_name)
      g$forest_group_name[match(code, g$forest_group)]
}
