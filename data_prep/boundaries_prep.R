# Build the focal-polygon library: data/boundaries.gpkg
#
# The app offers a fixed library of management-unit boundaries in place of
# user uploads. This script builds that library from the agencies' own
# services, at build time only -- the app never contacts them. Run from the
# project root:
#
#   Rscript data_prep/boundaries_prep.R
#
# Which units go in is listed in data_prep/boundaries_units.csv (source, unit),
# written by boundaries_select.R; where each source comes from is in
# boundaries_sources.R.
#
# Output: data/boundaries.gpkg, with two layers
#   units  one row per unit, EPSG:4326:
#          id (e.g. "nps:YELL", stable, for citable links), name, agency,
#          collection (the selector group), unit_type, state, source
#   about  one row: library version, build date, and sources
#
# Processing: units are read with read_source() (boundaries_sources.R); from
# ArcGIS layers they are fetched in batches, lightly generalized by the server
# (FETCH_OFFSET, a fifth of the simplification tolerance) to keep transfers
# small. Each unit is then merged to a single (multi)polygon, repaired, and
# simplified to SIMPLIFY_M in an equal-area projection. FIA plot coordinates
# are perturbed by up to about 1 km and the app's finest grain is 1 km, so
# detail below ~100 m carries no information here and only costs size.
#
# Distribution: the .gpkg is not committed; it is a piggyback release asset,
# fetched by ensure_boundaries() (helpers/boundaries.R). After a rebuild, bump
# VERSION here and in ensure_boundaries(), then upload:
#
#   piggyback::pb_release_create("matthewkling/treeshift", VERSION)
#   piggyback::pb_upload("data/boundaries.gpkg", "matthewkling/treeshift", tag = VERSION)

library(sf)
source("data_prep/boundaries_sources.R")

VERSION      <- "boundaries-v1.2"
SIMPLIFY_M   <- 100
FETCH_OFFSET <- 0.0002   # degrees, ~20 m
BATCH        <- 40
OUT          <- "data/boundaries.gpkg"

# Planar (GEOS) repair, as one MULTIPOLYGON. Repair can return a collection
# holding stray lines or points; only the polygons are kept.
repair <- function(geom) {
      old <- sf_use_s2(FALSE)
      on.exit(sf_use_s2(old))
      geom <- st_make_valid(geom)
      if (any(st_geometry_type(geom) == "GEOMETRYCOLLECTION")) {
            # Planar on purpose, so sf's note saying so is dropped
            geom <- suppressMessages(st_union(st_collection_extract(geom, "POLYGON")))
      }
      st_cast(geom, "MULTIPOLYGON")
}

# One unit's features -> one simplified (multi)polygon, with diagnostics
build_unit <- function(raw) {
      geom <- raw |>
            st_transform(5070) |>
            st_make_valid() |>
            st_union() |>
            st_simplify(preserveTopology = TRUE, dTolerance = SIMPLIFY_M) |>
            repair()
      area_raw <- as.numeric(sum(st_area(st_make_valid(st_transform(raw, 5070)))))
      area_new <- as.numeric(st_area(geom))
      # Reprojecting can introduce tiny self-intersections, so repair again
      # in the library's own coordinates
      list(geom = repair(st_transform(geom, 4326)),
           simplify_change = area_new / area_raw - 1,
           acres_check = if ("gis_acres" %in% names(raw))
                 area_new / 4046.856 / sum(raw$gis_acres) - 1 else NA_real_,
           vertices = nrow(st_coordinates(geom)))
}

manifest <- read.csv("data_prep/boundaries_units.csv", colClasses = "character")
stopifnot(!anyDuplicated(manifest), all(manifest$source %in% names(boundary_sources)))

built <- list()
for (s in unique(manifest$source)) {
      src <- boundary_sources[[s]]
      wanted <- manifest$unit[manifest$source == s]
      raw <- read_source(src, wanted, offset = FETCH_OFFSET, batch = BATCH)
      missing <- setdiff(wanted, raw[[src$key]])
      if (length(missing)) stop(s, ": no features for ", paste(missing, collapse = ", "))

      # Downloads above run one at a time (network calls crash in forked
      # processes on macOS); building units, the slow part, runs in parallel
      res <- parallel::mclapply(wanted, function(u) {
            r <- raw[raw[[src$key]] == u, ]
            b <- build_unit(r)
            attrs <- src$standardize(st_drop_geometry(r))[1, ]
            list(
                  sf = st_sf(id = paste0(s, ":", u), name = attrs$name,
                             agency = src$agency, collection = src$collection,
                             unit_type = attrs$unit_type, state = attrs$state,
                             source = source_location(src), geometry = b$geom),
                  diag = data.frame(source = s, unit = u, name = attrs$name,
                                    simplify_change = b$simplify_change,
                                    acres_check = b$acres_check, vertices = b$vertices))
      }, mc.cores = parallel::detectCores() - 2)
      failed <- !vapply(res, function(r) is.list(r) && !inherits(r, "try-error"), NA)
      if (any(failed)) stop(s, ": building failed for ", paste(wanted[failed], collapse = ", "), ":\n",
                            paste(unique(vapply(res[failed], paste, "", collapse = " ")), collapse = "\n"))
      built <- c(built, res)
      message(sprintf("%-8s %4d units fetched and built", s, length(wanted)))
}

units <- do.call(rbind, lapply(built, `[[`, "sf"))
# Collections in the order the sources are defined, which is the order the
# app's selector lists them in; units alphabetical within each
collections <- vapply(boundary_sources, `[[`, "", "collection")
units <- units[order(match(units$collection, collections), units$name), ]
diag <- do.call(rbind, lapply(built, `[[`, "diag"))

# Valid both as spherical (s2, sf's default for lon/lat) and planar (GEOS, as
# terra uses them in the app) geometries
stopifnot(!anyDuplicated(units$id), !any(st_is_empty(units)), all(st_is_valid(units)),
          all(st_geometry_type(units) == "MULTIPOLYGON"))
sf_use_s2(FALSE)
stopifnot(all(st_is_valid(units)))
sf_use_s2(TRUE)
message(sprintf("simplification area change: median %+.3f%%, largest %+.2f%% (%s)",
                100 * median(diag$simplify_change),
                100 * diag$simplify_change[which.max(abs(diag$simplify_change))],
                diag$name[which.max(abs(diag$simplify_change))]))
chk <- diag[!is.na(diag$acres_check), ]
message(sprintf("area vs Forest Service gis_acres: median %+.2f%%, largest %+.2f%% (%s)",
                100 * median(chk$acres_check),
                100 * chk$acres_check[which.max(abs(chk$acres_check))],
                chk$name[which.max(abs(chk$acres_check))]))
message(sprintf("vertices per unit: median %d, max %d (%s)", as.integer(median(diag$vertices)),
                max(diag$vertices), diag$name[which.max(diag$vertices)]))

about <- data.frame(
      version = VERSION,
      built = format(Sys.Date()),
      simplify_m = SIMPLIFY_M,
      sources = paste(sprintf("%s: %s (retrieved %s)", names(boundary_sources),
                              vapply(boundary_sources, source_location, ""), Sys.Date()),
                      collapse = "; ")
)

if (file.exists(OUT)) invisible(file.remove(OUT))
st_write(units, OUT, layer = "units", quiet = TRUE)
st_write(about, OUT, layer = "about", quiet = TRUE)
message(sprintf("Wrote %s: %d units, %.1f MB", OUT, nrow(units), file.size(OUT) / 1e6))
