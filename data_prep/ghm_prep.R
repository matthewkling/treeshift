# Global Human Modification (gHM) v3, onto the app's climate grid.
#
#   Rscript data_prep/ghm_prep.R <path to HMv20240801_2022s_AA_300.tif>
#
# Zenodo serves the 9.3 GB file slowly, so it can instead be read remotely,
# fetching only the contiguous US. Zenodo refuses GDAL's default user agent,
# hence the override, and drops some tile requests, hence the retries:
#
#   GDAL_HTTP_USERAGENT=curl/8.7.1 GDAL_DISABLE_READDIR_ON_OPEN=EMPTY_DIR \
#   GDAL_HTTP_MAX_RETRY=10 GDAL_HTTP_RETRY_DELAY=5 \
#   CPL_VSIL_CURL_ALLOWED_EXTENSIONS=.tif Rscript data_prep/ghm_prep.R \
#   /vsicurl/https://zenodo.org/records/14502573/files/HMv20240801_2022s_AA_300.tif
#
# The window is fetched in strips of one degree of latitude, each saved to
# data/ghm_strips/ as it completes and retried on its own if a read fails, so
# an interrupted run picks up where it stopped when rerun. Delete the strips
# once data/ghm.tif is built.
#
# Run from the project root. The source is the 2022 static snapshot of all
# threats combined (0 = unmodified, 1 = fully modified; 300 m, EPSG:4326),
# from Theobald et al. (2024), Global human modification datasets of
# terrestrial ecosystems for 2022, doi:10.5281/zenodo.14502573, CC BY 4.0.
# The connectivity module uses it to lower conductance where land is modified
# (server/connectivity.R).
#
# Output: data/ghm.tif, the mean gHM in each cell of data/clim_historical.tif,
# NA outside the climate data.

library(terra)

src <- commandArgs(trailingOnly = TRUE)[1]
# A local copy, or a GDAL remote path (/vsicurl/...), which reads only the
# blocks the crop needs
if (is.na(src) || (!startsWith(src, "/vsi") && !file.exists(src)))
      stop("Give the path to HMv20240801_2022s_AA_300.tif.")

template <- rast("data/clim_historical.tif")[[1]]
window <- extend(ext(template), 0.05)
ghm_src <- rast(src)

strip_dir <- "data/ghm_strips"
dir.create(strip_dir, showWarnings = FALSE)
edges <- seq(floor(ymin(window)), ceiling(ymax(window)), by = 1)

# One strip, written under a temporary name and renamed when complete, so a
# failed read never leaves a strip that looks finished
fetch_strip <- function(y0, attempts = 5) {
      out <- file.path(strip_dir, sprintf("ghm_%+03d.tif", y0))
      if (file.exists(out)) return(out)
      e <- ext(xmin(window), xmax(window), max(y0, ymin(window)), min(y0 + 1, ymax(window)))
      tmp <- paste0(out, ".part.tif")
      for (i in seq_len(attempts)) {
            ok <- tryCatch({
                  crop(ghm_src, e, filename = tmp, overwrite = TRUE)
                  TRUE
            }, error = function(err) {
                  message(sprintf("Strip %d, attempt %d failed: %s", y0, i, conditionMessage(err)))
                  FALSE
            })
            if (ok) {
                  file.rename(tmp, out)
                  message(sprintf("Strip %d done", y0))
                  return(out)
            }
            Sys.sleep(30 * i)
      }
      stop(sprintf("Strip %d failed %d times; rerun to try again.", y0, attempts))
}

strips <- vapply(head(edges, -1), fetch_strip, character(1))

ghm <- vrt(strips) |>
      resample(template, method = "average", threads = TRUE) |>
      mask(template) |>
      setNames("ghm")

writeRaster(ghm, "data/ghm.tif", datatype = "FLT4S",
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"), overwrite = TRUE)
