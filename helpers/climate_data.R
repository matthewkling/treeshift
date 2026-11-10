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
