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
# Processing: units are fetched in batches, lightly generalized by the server
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

VERSION      <- "boundaries-v1"
SIMPLIFY_M   <- 100
FETCH_OFFSET <- 0.0002   # degrees, ~20 m
BATCH        <- 40
OUT          <- "data/boundaries.gpkg"

# Features for a batch of units from an ArcGIS REST layer, as sf in EPSG:4326.
# Downloaded to a file first: the Forest Service server resets GDAL's own HTTP
# requests.
fetch_batch <- function(src, units) {
      where <- sprintf("%s IN (%s)", src$key, paste0("'", units, "'", collapse = ","))
      url <- paste0(src$url, "/query?",
                    "where=", utils::URLencode(where, reserved = TRUE),
                    "&outFields=", paste(src$fields, collapse = ","),
                    "&outSR=4326&maxAllowableOffset=", FETCH_OFFSET, "&f=geojson")
      tmp <- tempfile(fileext = ".geojson")
      on.exit(unlink(tmp))
      utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
      st_read(tmp, quiet = TRUE)
}

# One unit's features -> one simplified (multi)polygon, with diagnostics
build_unit <- function(raw) {
      geom <- raw |>
            st_transform(5070) |>
            st_make_valid() |>
            st_union() |>
            st_simplify(preserveTopology = TRUE, dTolerance = SIMPLIFY_M) |>
            st_make_valid()
      # Repair can return a collection holding stray lines or points; keep
      # only the polygons
      if (any(st_geometry_type(geom) == "GEOMETRYCOLLECTION")) {
            geom <- st_union(st_collection_extract(geom, "POLYGON"))
      }
      geom <- st_cast(geom, "MULTIPOLYGON")
      area_raw <- as.numeric(sum(st_area(st_make_valid(st_transform(raw, 5070)))))
      area_new <- as.numeric(st_area(geom))
      list(geom = st_transform(geom, 4326),
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
      batches <- split(wanted, ceiling(seq_along(wanted) / BATCH))
      raw <- do.call(rbind, lapply(batches, fetch_batch, src = src))
      missing <- setdiff(wanted, raw[[src$key]])
      if (length(missing)) stop(s, ": no features for ", paste(missing, collapse = ", "))

      for (u in wanted) {
            r <- raw[raw[[src$key]] == u, ]
            b <- build_unit(r)
            attrs <- src$standardize(st_drop_geometry(r))[1, ]
            built[[length(built) + 1]] <- list(
                  sf = st_sf(id = paste0(s, ":", u), name = attrs$name,
                             agency = src$agency, collection = src$collection,
                             unit_type = attrs$unit_type, state = attrs$state,
                             source = src$url, geometry = b$geom),
                  diag = data.frame(source = s, unit = u, name = attrs$name,
                                    simplify_change = b$simplify_change,
                                    acres_check = b$acres_check, vertices = b$vertices))
      }
      message(sprintf("%-8s %4d units fetched and built", s, length(wanted)))
}

units <- do.call(rbind, lapply(built, `[[`, "sf"))
# Collections in the order the sources are defined, which is the order the
# app's selector lists them in; units alphabetical within each
collections <- vapply(boundary_sources, `[[`, "", "collection")
units <- units[order(match(units$collection, collections), units$name), ]
diag <- do.call(rbind, lapply(built, `[[`, "diag"))

stopifnot(!anyDuplicated(units$id), all(st_is_valid(units)), !any(st_is_empty(units)))
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
                              vapply(boundary_sources, `[[`, "", "url"), Sys.Date()),
                      collapse = "; ")
)

if (file.exists(OUT)) invisible(file.remove(OUT))
st_write(units, OUT, layer = "units", quiet = TRUE)
st_write(about, OUT, layer = "about", quiet = TRUE)
message(sprintf("Wrote %s: %d units, %.1f MB", OUT, nrow(units), file.size(OUT) / 1e6))
