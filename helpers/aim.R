# Analog impact surfaces --------------------------------------------------------
#
# AIM surfaces are computed on plot and cell coordinates in CONUS Albers
# equal-area (EPSG:5070), in km so the geographic bandwidth keeps its units. The
# FFT method needs planar coordinates, and planar distance also makes the exact
# method about 1.5x faster than great-circle distance; over CONUS the two differ
# by about 1%. Results are written back onto the lon/lat modeling grid.
aim_crs <- "EPSG:5070"

aim_xy <- function(xy) {
      p <- terra::project(as.matrix(xy[, c("x", "y")]), "EPSG:4326", aim_crs) / 1000
      colnames(p) <- c("x", "y")
      p
}

aim_kernels <- function(theta_clim, theta_geog) {
      list(env  = analogs::kernel(max = theta_clim * 3, weight = "gaussian", theta = theta_clim),
           geog = analogs::kernel(max = theta_geog * 3, weight = "gaussian", theta = theta_geog))
}

# FFT runs are tiled geographically, so each tile's grid spans only the climates
# found there. ~300 km tiles were fastest in benchmarks.
aim_tile_km <- 300

# Which method a run uses: "exact" visits every cell-plot pair, "fft" (binned
# convolution) costs about the same whatever the grain, and its error is small
# (99.9th percentile ~0.007 on a 0-1 occurrence scale; identical NA cells). Exact
# wins at coarse grain or narrow bandwidths, FFT at fine grain and wide ones.
#
# The choice must be a deterministic function of the settings, so that the same
# settings always give the same map -- hence estimated work, not timings. Exact
# work is the number of cell-plot pairs within the kernel limits, from analog
# counts at a fixed sample of cells. FFT work is the number of grid nodes, from
# each tile's extent in geography and climate (nodes are a quarter-bandwidth
# apart, plus a 3-bandwidth halo), times the number of kernel sums. The
# coefficients were fitted on runs at 1-10 km over 6-20 degree domains with
# 1-80 outcome columns; only their ratio matters, so machine speed cancels.
# One method serves every timeframe of a run, so differences between timeframes
# never mix methods.
aim_method <- function(pts_list, index, theta_clim, theta_geog, n_y, n_sums) {
      env_vars <- setdiff(names(pts_list[[1]]), c("x", "y"))
      if (length(env_vars) > 3) return("exact")   # FFT supports up to 3

      est <- vapply(pts_list, function(pts) {
            n <- nrow(pts)
            probe <- pts[unique(round(seq(1, n, length.out = min(n, 2000)))), , drop = FALSE]
            count <- suppressWarnings(analogs::analog_search(
                  x = probe, pool = index, select = "all", stat = "count",
                  env = analogs::kernel(max = theta_clim * 3),
                  geog = analogs::kernel(max = theta_geog * 3)))$count
            pairs <- n * mean(count, na.rm = TRUE)

            spans <- pts %>%
                  group_by(tx = floor(x / aim_tile_km), ty = floor(y / aim_tile_km)) %>%
                  summarise(across(everything(), ~ max(.x) - min(.x)), .groups = "drop")
            nodes_along <- function(span, theta) (span + 6 * theta) / (0.25 * theta) + 1
            nodes <- sum(nodes_along(spans$x, theta_geog) * nodes_along(spans$y, theta_geog) *
                         Reduce(`*`, lapply(spans[env_vars], nodes_along, theta = theta_clim)))

            c(exact = 5.0e-6 * n + 4.5e-9 * pairs * (1 + 0.20 * (n_y - 1)),
              fft   = 1.8e-6 * n + 2.1e-9 * nodes * (1 + 0.42 * n_sums))
      }, numeric(2))
      if (sum(est["fft", ]) < sum(est["exact", ])) "fft" else "exact"
}

# AIM surfaces for several focal rasters (e.g. historical and future climate)
# against one pool, on each raster's own grid. `pool` is pool()'s list, of which
# this uses the AIM-coordinate plot table (`aim_pts`) and its index
# (`aim_index`). Returns a list of SpatRasters named like `focals`, with the
# layers analog_impact() would give for a raster.
aim_surfaces <- function(focals, pool, y, theta_clim, theta_geog, stat) {
      k <- aim_kernels(theta_clim, theta_geog)
      cells <- lapply(focals, function(r) as.data.frame(r, xy = TRUE, cells = TRUE, na.rm = TRUE))
      pts <- lapply(cells, function(fd) data.frame(aim_xy(fd), fd[, setdiff(names(fd), c("cell", "x", "y")), drop = FALSE]))

      run <- function(method) lapply(pts, function(p) {
            if (method == "exact") {
                  suppressWarnings(analogs::analog_impact(
                        x = p, pool = pool$aim_index, env = k$env, geog = k$geog,
                        stat = stat, y = y, coord_type = "projected"))
            } else {
                  span <- vapply(p[, c("x", "y")], function(v) diff(range(v)), numeric(1))
                  suppressWarnings(analogs::tiled_analog_search(
                        x = p, pool = pool$aim_pts, y = y,
                        n_tiles = max(1, round(prod(span) / aim_tile_km^2)),
                        fun = analogs::analog_impact, env = k$env, geog = k$geog,
                        stat = stat, method = "fft", coord_type = "projected",
                        progress = FALSE))
            }
      })

      # Kernel sums FFT computes: one per outcome column (two with
      # weighted_sd), the weights, and the squared weights if a stat needs them
      n_sums <- ncol(y) * (1 + ("weighted_sd" %in% stat)) + 1 +
            any(c("ess", "effective_weight") %in% stat)
      method <- aim_method(pts, pool$aim_index, theta_clim, theta_geog,
                           n_y = ncol(y), n_sums = n_sums)
      res <- if (method == "fft") tryCatch(run("fft"), error = function(e) {
            # A tile whose FFT grid exceeds the size limit (narrow bandwidths
            # over climatically varied ground) sends the whole run to exact,
            # rather than mixing methods between tiles. Any other error is real.
            if (!grepl("max_cells", conditionMessage(e))) stop(e)
            message("AIM: FFT grid too large for a tile; using the exact method")
            NULL
      })
      if (is.null(res)) res <- run("exact")

      Map(function(r, fd, out) {
            cols <- setdiff(names(out), c("index", "x", "y"))
            s <- terra::rast(r, nlyrs = length(cols), names = cols, vals = NA)
            s[fd$cell] <- as.matrix(out[, cols, drop = FALSE])
            s
      }, focals, cells, res)
}
