# Choose which units go in the focal-polygon library, and write them to the
# manifest data_prep/boundaries_units.csv that boundaries_prep.R builds from.
#
# The manifest is committed, and this script is rerun only to change the
# selection policy, not on every build. That keeps unit ids stable: a park
# re-surveyed at 99 km2 does not silently drop out of the library and break the
# links that cite it.
#
# Policy (lower 48 only -- there is no climate data elsewhere):
#   usfs     every National Forest System administrative unit
#   usfs_rd  every ranger district
#   nps      National Park Service units of at least NPS_MIN_KM2. Smaller units
#            (historic sites, memorials, battlefields) are a fraction of one map
#            cell at the default 10 km grain and typically contain no FIA plots.
# and, for all sources, only units with at least one 1 km climate cell inside:
# a unit that is essentially open water (Dry Tortugas) has nothing to summarize.
#
# Run from the project root:  Rscript data_prep/boundaries_select.R

library(sf)

NPS_MIN_KM2 <- 100
MANIFEST    <- "data_prep/boundaries_units.csv"
source("data_prep/boundaries_sources.R")

# Unit centroids and areas from a coarsely generalized copy of each layer: enough
# to apply the policy, not the geometry the library is built from.
coarse <- function(src) {
      url <- paste0(src$url, "/query?where=1%3D1&outFields=", src$key,
                    "&outSR=4326&maxAllowableOffset=0.002&f=geojson")
      tmp <- tempfile(fileext = ".geojson")
      on.exit(unlink(tmp))
      utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
      # Per feature, in the equal-area projection, then totaled by unit (a
      # unit can span several features); no geometry merging needed
      x <- st_collection_extract(st_make_valid(st_transform(st_read(tmp, quiet = TRUE), 5070)),
                                 "POLYGON", warn = FALSE)
      # Compare coordinates with the bounds directly: as a spherical polygon,
      # a lon/lat box has great-circle edges, and its southern edge bows far
      # enough north to drop south Florida and south Texas
      p <- st_coordinates(st_transform(st_point_on_surface(x), 4326))
      inside <- p[, "X"] >= lower48["xmin"] & p[, "X"] <= lower48["xmax"] &
            p[, "Y"] >= lower48["ymin"] & p[, "Y"] <= lower48["ymax"]
      # 1 km climate cells whose centers fall inside
      cells <- terra::extract(clim, terra::vect(st_transform(x, 4326)), ID = TRUE)
      cells <- tabulate(cells$ID[!is.na(cells[, 2])], nrow(x))
      f <- data.frame(unit = x[[src$key]], km2 = as.numeric(st_area(x)) / 1e6,
                      inside = inside, cells = cells)
      data.frame(unit = unique(f$unit),
                 lower48 = tapply(f$inside, f$unit, any)[unique(f$unit)],
                 km2 = tapply(f$km2, f$unit, sum)[unique(f$unit)],
                 cells = tapply(f$cells, f$unit, sum)[unique(f$unit)])
}

# The climate rasters' extent, and the 1 km climate grid itself
lower48 <- c(xmin = -125, xmax = -66, ymin = 24, ymax = 50)
clim <- terra::rast("data/clim_historical.tif")[[1]]

picked <- do.call(rbind, lapply(names(boundary_sources), function(s) {
      u <- coarse(boundary_sources[[s]])
      keep <- u$lower48 & u$cells > 0 & (s != "nps" | u$km2 >= NPS_MIN_KM2)
      dropped_water <- u$unit[u$lower48 & u$cells == 0 & (s != "nps" | u$km2 >= NPS_MIN_KM2)]
      message(sprintf("%-8s %4d units, %4d in the lower 48, %4d selected%s",
                      s, nrow(u), sum(u$lower48), sum(keep),
                      if (length(dropped_water)) paste0(" (no climate cells: ", paste(dropped_water, collapse = ", "), ")") else ""))
      data.frame(source = s, unit = sort(u$unit[keep]))
}))

write.csv(picked, MANIFEST, row.names = FALSE, quote = FALSE)
message(sprintf("Wrote %s: %d units", MANIFEST, nrow(picked)))
