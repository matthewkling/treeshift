# Real Omniscape on the 64 x 64 crop used by the exact factorial decomposition,
# to validate its R reimplementation of Omniscape.
#
#   Rscript sandbox/dipole_factorial_omniscape.R
#
# Run from the project root. Crop: rows 305-368, columns 337-400 of
# sandbox/connectivity_bench/suitability.tif (2 km P. edulis). Conductance is
# suitability / full-surface mean + 1e-3, sources suitability with values
# below 0.01 set to 0, radius 6 and 12 pixels, block size 1 (exact).
# Writes omniscape_R6.tif and omniscape_R12.tif to
# sandbox/connectivity_bench/factorial/.
library(terra); library(circuitscaper)
out_dir <- "sandbox/connectivity_bench/factorial"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
suit_full <- rast("sandbox/connectivity_bench/suitability.tif")
mean_full <- global(suit_full, "mean", na.rm = TRUE)[[1]]
suit <- suit_full[305:368, 337:400, drop = FALSE]
conductance <- suit / mean_full + 1e-3
sources <- suit; sources[sources < 0.01] <- 0
writeRaster(suit, file.path(out_dir, "crop_suitability.tif"), overwrite = TRUE)
for (R in c(6L, 12L)) {
      t0 <- Sys.time()
      out <- os_run(1 / conductance, radius = R, source_strength = sources, block_size = 1L,
                    source_threshold = 0.01, resistance_is = "resistances", solver = "cholmod")
      writeRaster(out, file.path(out_dir, sprintf("omniscape_R%d.tif", R)), overwrite = TRUE)
      message(sprintf("R = %d: %.0f s, layers %s", R, as.numeric(difftime(Sys.time(), t0, units = "secs")),
                      paste(names(out), collapse = ", ")))
}
