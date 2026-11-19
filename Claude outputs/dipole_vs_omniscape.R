# Dipole method vs Omniscape across kernel radii.
#
#   Rscript sandbox/dipole_vs_omniscape.R
#
# Run from the project root. Needs sandbox/dipole.R and
# sandbox/connectivity_bench/suitability.tif (from sandbox/export_surface.R).
#
# For each radius, times and compares:
#   omniscape (practical)  block size ~ radius / 10 (odd), as usually run
#   omniscape (exact)      block size 1, at the radii in `exact_radii` only;
#                          at 50 km this is the long job (perhaps 15-30 min)
#   dipole                 32 draws of 1e5 dipoles, with flow potential and
#                          normalization
#
# Every job writes its GeoTIFF and a row in timings.csv as it finishes. If the
# script is interrupted, rerunning it skips finished jobs and carries on.
# Set RUN_OMNISCAPE=false to run only the dipole jobs.
#
# Sources are suitability with values below 0.01 set to zero for both methods
# (Omniscape's source_threshold), and resistance is 1 / conductance as in
# server/connectivity.R.

library(terra)
source("sandbox/dipole.R")

radii <- c(25, 50, 100, 200)
exact_radii <- c(25, 50)
threads <- 4L
source_threshold <- 0.01
solver <- "cg+amg"      # faster than cholmod for Omniscape's larger windows here
n_draws <- 32
n_dipoles <- 1e5
run_omniscape <- !identical(tolower(Sys.getenv("RUN_OMNISCAPE", "true")), "false")

out_dir <- "sandbox/connectivity_bench/radius"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
timings_file <- file.path(out_dir, "timings.csv")
set.seed(1)


# Inputs ----------------------------------------------------------------------

suit <- rast("sandbox/connectivity_bench/suitability.tif")
conductance <- suit / global(suit, "mean", na.rm = TRUE)[[1]] + 1e-3
resistance <- 1 / conductance
sources <- suit
sources[sources < source_threshold] <- 0
px_km <- sqrt(global(cellSize(suit, unit = "km"), "mean", na.rm = TRUE)[[1]])

# Omniscape's usual rule of thumb: an odd block size of about radius / 10
practical_block <- function(radius_px) {
      b <- max(1L, as.integer(floor(radius_px / 10)))
      if (b %% 2 == 0) b - 1L else b
}

jobs <- list()
for (r in radii) {
      rpx <- as.integer(round(r / px_km))
      blocks <- unique(c(practical_block(rpx), if (r %in% exact_radii) 1L))
      if (run_omniscape) for (b in blocks)
            jobs[[length(jobs) + 1]] <- list(method = "omniscape", radius = r, radius_px = rpx, block = b,
                                             file = sprintf("omniscape_r%d_b%d.tif", r, b))
      jobs[[length(jobs) + 1]] <- list(method = "dipole", radius = r, radius_px = rpx, block = NA,
                                       file = sprintf("dipole_r%d.tif", r))
}

done <- if (file.exists(timings_file)) read.csv(timings_file) else NULL
is_done <- function(job) !is.null(done) && job$file %in% done$file && file.exists(file.path(out_dir, job$file))
record <- function(job, seconds, detail) {
      row <- data.frame(method = job$method, radius_km = job$radius, radius_px = job$radius_px,
                        block = job$block, seconds = round(seconds, 1), detail = detail, file = job$file)
      write.table(row, timings_file, sep = ",", row.names = FALSE,
                  col.names = !file.exists(timings_file), append = file.exists(timings_file))
      message(sprintf("%-10s r = %3d km  block %-3s %8.1f s   %s", job$method, job$radius,
                      ifelse(is.na(job$block), "-", job$block), seconds, detail))
}


# Run -------------------------------------------------------------------------

todo <- Filter(Negate(is_done), jobs)
message(sprintf("%d jobs, %d left to run", length(jobs), length(todo)))

if (any(vapply(todo, function(j) j$method == "omniscape", logical(1)))) {
      library(circuitscaper)
      cs_setup(threads = threads)
      # Warm up Julia's compilation on a small crop, so it isn't in the timings
      mid <- c(mean(as.vector(ext(suit))[1:2]), mean(as.vector(ext(suit))[3:4]))
      e <- ext(mid[1] - 0.75, mid[1] + 0.75, mid[2] - 0.75, mid[2] + 0.75)
      invisible(os_run(crop(resistance, e), radius = 3L, source_strength = crop(sources, e),
                       block_size = 1L, source_threshold = source_threshold,
                       resistance_is = "resistances", parallelize = threads > 1,
                       julia_threads = threads, solver = solver))
}

if (any(vapply(todo, function(j) j$method == "dipole", logical(1)))) {
      t0 <- Sys.time()
      g <- dp_graph(conductance)
      g0 <- dp_graph(mask(setValues(conductance, 1), conductance))
      t_graphs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      message(sprintf("Dipole graphs and factorizations: %.1f s (counted in each dipole job)", t_graphs))
}

for (job in todo) {
      t0 <- Sys.time()
      if (job$method == "omniscape") {
            out <- os_run(resistance, radius = job$radius_px, source_strength = sources,
                          block_size = job$block, source_threshold = source_threshold,
                          resistance_is = "resistances", parallelize = threads > 1,
                          julia_threads = threads, solver = solver)
            secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
            detail <- sprintf("%d threads, %s", threads, solver)
      } else {
            out <- dp_connectivity(sources, conductance, job$radius, n_draws, n_dipoles, g = g, g0 = g0)
            secs <- as.numeric(difftime(Sys.time(), t0, units = "secs")) + t_graphs
            h <- attr(out, "halves")
            detail <- sprintf("%d draws x %g dipoles; halves rel. diff: current %.2f, normalized %.2f",
                              n_draws, n_dipoles, h["current", 2], h["normalized", 2])
      }
      writeRaster(out, file.path(out_dir, job$file), overwrite = TRUE)
      record(job, secs, detail)
}


# Compare ---------------------------------------------------------------------

timings <- read.csv(timings_file)
v <- function(x) values(x, mat = FALSE)
hab <- v(suit) >= source_threshold & !is.na(v(suit))
rho <- function(x, y) round(cor(v(x)[hab], v(y)[hab], method = "spearman", use = "complete.obs"), 3)
layers <- c(current = "cumulative_current", flow_potential = "flow_potential", normalized = "normalized_current")
load <- function(f) { p <- file.path(out_dir, f); if (file.exists(p)) rast(p) }

cmp <- list()
for (r in radii) {
      rows <- timings[timings$radius_km == r, ]
      dp <- load(sprintf("dipole_r%d.tif", r))
      os_rows <- rows[rows$method == "omniscape", ]
      for (k in seq_len(nrow(os_rows))) {
            os <- load(os_rows$file[k])
            if (is.null(os) || is.null(dp)) next
            cmp[[length(cmp) + 1]] <- data.frame(
                  radius_km = r, comparison = sprintf("dipole vs omniscape block %d", os_rows$block[k]),
                  current = rho(dp$current, os[[layers["current"]]]),
                  flow_potential = rho(dp$flow_potential, os[[layers["flow_potential"]]]),
                  normalized = rho(dp$normalized, os[[layers["normalized"]]]),
                  speedup = round(os_rows$seconds[k] / rows$seconds[rows$method == "dipole"], 1))
      }
      # Omniscape's own block approximation, where the exact run exists
      if (nrow(os_rows) == 2) {
            a <- load(os_rows$file[os_rows$block > 1]); b <- load(os_rows$file[os_rows$block == 1])
            cmp[[length(cmp) + 1]] <- data.frame(
                  radius_km = r, comparison = sprintf("omniscape block %d vs block 1", max(os_rows$block)),
                  current = rho(a[[layers["current"]]], b[[layers["current"]]]),
                  flow_potential = rho(a[[layers["flow_potential"]]], b[[layers["flow_potential"]]]),
                  normalized = rho(a[[layers["normalized"]]], b[[layers["normalized"]]]),
                  speedup = round(os_rows$seconds[os_rows$block == 1] / max(os_rows$seconds[os_rows$block > 1]), 1))
      }
}
message("\nTimings (s):")
print(timings[, c("method", "radius_km", "block", "seconds")], row.names = FALSE)
if (length(cmp)) {
      cmp <- do.call(rbind, cmp)
      write.csv(cmp, file.path(out_dir, "comparison.csv"), row.names = FALSE)
      message("\nSpearman rho within habitat. speedup: omniscape time / dipole time, or for the")
      message("'block b vs block 1' rows, block-1 time / block-b time.")
      print(cmp, row.names = FALSE)
}

# Normalized maps, rank-stretched: practical Omniscape vs dipole, by radius
rank_stretch <- function(r) { x <- v(r); ok <- is.finite(x); x[ok] <- rank(x[ok]) / sum(ok); x[!ok] <- NA; setValues(r, x) }
panels <- list()
for (r in radii) {
      os_rows <- timings[timings$radius_km == r & timings$method == "omniscape", ]
      os <- if (nrow(os_rows)) load(os_rows$file[which.max(os_rows$block)])
      dp <- load(sprintf("dipole_r%d.tif", r))
      if (!is.null(os)) panels[[length(panels) + 1]] <- list(os[[layers["normalized"]]],
            sprintf("Omniscape r = %d km, block %d: normalized", r, max(os_rows$block)))
      if (!is.null(dp)) panels[[length(panels) + 1]] <- list(dp$normalized, sprintf("Dipole r = %d km: normalized", r))
}
if (length(panels)) {
      nc <- if (run_omniscape) 2 else 1
      png(file.path(out_dir, "normalized_by_radius.png"), width = 1000 * nc, height = 900 * ceiling(length(panels) / nc), res = 150)
      par(mfrow = c(ceiling(length(panels) / nc), nc), mar = c(1, 1, 3, 1))
      for (p in panels) plot(rank_stretch(p[[1]]), col = hcl.colors(100, "viridis"), legend = FALSE, axes = FALSE, main = p[[2]])
      dev.off()
}
message("Wrote ", out_dir)
