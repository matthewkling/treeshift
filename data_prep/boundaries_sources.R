# Where each kind of focal-polygon library unit comes from, and how to read it.
# Shared by boundaries_select.R and boundaries_prep.R.
#
# Each source is either an ArcGIS REST layer (`url` ending in the layer number),
# queried for just the units wanted, or a zipped shapefile (`zip`, `layer`),
# downloaded whole. It names the field holding the publisher's stable unit code
# (the library id is "<source>:<code>", e.g. "nps:YELL"), the fields to read,
# and how to map them onto the library's attributes (`numeric_key` when that
# field is numeric, so queries don't quote it). `collection` is the group
# a unit is listed under in the app's selector; sources are listed here in the
# order the selector lists them.
#
# `standardize` gets all of a unit's features (a unit can span several) and the
# first row of its result is kept, so per-unit summaries belong there.

census_year <- 2025

boundary_sources <- list(
      nps = list(
            agency = "National Park Service",
            collection = "National Park Service units",
            url = paste0("https://services1.arcgis.com/fBc8EJBxQRMcHlei/arcgis/rest/services/",
                         "NPS_Land_Resources_Division_Boundary_and_Tract_Data_Service/FeatureServer/2"),
            key = "UNIT_CODE",
            fields = c("UNIT_CODE", "UNIT_NAME", "UNIT_TYPE", "STATE"),
            standardize = function(x) data.frame(
                  unit = x$UNIT_CODE, name = x$UNIT_NAME,
                  unit_type = x$UNIT_TYPE, state = x$STATE)
      ),
      usfs = list(
            agency = "USDA Forest Service",
            collection = "National Forests and Grasslands",
            url = "https://apps.fs.usda.gov/arcx/rest/services/EDW/EDW_ForestSystemBoundaries_01/MapServer/0",
            key = "forestorgcode",
            fields = c("forestorgcode", "forestname", "gis_acres"),
            standardize = function(x) data.frame(
                  unit = x$forestorgcode, name = x$forestname,
                  unit_type = "National Forest System unit", state = NA_character_)
      ),
      usfs_rd = list(
            agency = "USDA Forest Service",
            collection = "Forest Service ranger districts",
            url = "https://apps.fs.usda.gov/arcx/rest/services/EDW/EDW_RangerDistricts_01/MapServer/0",
            key = "districtorgcode",
            fields = c("districtorgcode", "districtname", "forestname"),
            # A few district names recur across forests ("Dillon Ranger
            # District"), so the forest is part of the name
            standardize = function(x) data.frame(
                  unit = x$districtorgcode, name = paste0(x$districtname, ", ", x$forestname),
                  unit_type = "Ranger district", state = NA_character_)
      ),
      blm = list(
            agency = "Bureau of Land Management",
            collection = "BLM field offices",
            url = "https://gis.blm.gov/arcgis/rest/services/admin_boundaries/BLM_Natl_AdminUnit/MapServer/3",
            key = "ADM_UNIT_CD",
            fields = c("ADM_UNIT_CD", "ADMU_NAME", "ADMIN_ST"),
            standardize = function(x) data.frame(
                  unit = x$ADM_UNIT_CD,
                  name = paste0(blm_office_name(x$ADMU_NAME), ", ", x$ADMIN_ST),
                  unit_type = "BLM field office", state = x$ADMIN_ST)
      ),
      wilderness = list(
            agency = "Wilderness Connect (University of Montana), for the managing agencies",
            collection = "Wilderness areas",
            url = paste0("https://services1.arcgis.com/ERdCHt0sNM6dENSD/arcgis/rest/services/",
                         "Wilderness_Areas_in_the_United_States/FeatureServer/0"),
            key = "WID",
            numeric_key = TRUE,
            fields = c("WID", "NAME", "STATE", "Agency"),
            # One feature per managing agency within a wilderness; the unit is
            # the whole wilderness. The state goes in the name since a few
            # wilderness names recur.
            standardize = function(x) data.frame(
                  unit = as.character(x$WID), name = paste0(x$NAME, ", ", x$STATE),
                  unit_type = paste0("Wilderness (", paste(sort(unique(x$Agency)), collapse = ", "), ")"),
                  state = x$STATE)
      ),
      state = list(
            agency = "U.S. Census Bureau",
            collection = "States",
            zip = sprintf("https://www2.census.gov/geo/tiger/GENZ%d/shp/cb_%d_us_state_500k.zip",
                          census_year, census_year),
            layer = sprintf("cb_%d_us_state_500k", census_year),
            key = "STUSPS",
            fields = c("STUSPS", "NAME"),
            standardize = function(x) data.frame(
                  unit = x$STUSPS, name = x$NAME, unit_type = "State", state = x$STUSPS)
      ),
      county = list(
            agency = "U.S. Census Bureau",
            collection = "Counties",
            zip = sprintf("https://www2.census.gov/geo/tiger/GENZ%d/shp/cb_%d_us_county_500k.zip",
                          census_year, census_year),
            layer = sprintf("cb_%d_us_county_500k", census_year),
            key = "GEOID",
            fields = c("GEOID", "NAMELSAD", "STATE_NAME", "STUSPS"),
            standardize = function(x) data.frame(
                  unit = x$GEOID, name = paste0(x$NAMELSAD, ", ", x$STATE_NAME),
                  unit_type = "County or equivalent", state = x$STUSPS)
      ),
      eco_l3 = list(
            agency = "U.S. Environmental Protection Agency",
            collection = "EPA Level III ecoregions",
            url = "https://geodata.epa.gov/arcgis/rest/services/ORD/USEPA_Ecoregions_Level_III_and_IV/MapServer/11",
            key = "US_L3CODE",
            fields = c("US_L3CODE", "US_L3NAME"),
            standardize = function(x) data.frame(
                  unit = x$US_L3CODE, name = x$US_L3NAME,
                  unit_type = "Level III ecoregion", state = NA_character_)
      )
)

# BLM office names are upper case in the source ("GRAND JUNCTION FIELD
# OFFICE"), and some are abbreviated or truncated ("... FIELD OFFI", "... FO",
# "... NM") or lack "Field Office" altogether ("MOAB")
blm_office_name <- function(x) {
      x <- stringr::str_to_title(x)
      x <- stringr::str_replace(x, " (Fo|Field Offi|Field Offic)$", " Field Office")
      x <- stringr::str_replace_all(x, c("\\bNm\\b" = "National Monument", "^Nw " = "Northwest ",
                                         " And " = " and ", "D'alene" = "d'Alene"))
      ifelse(grepl("(Field Office|National Monument)$", x), x, paste(x, "Field Office"))
}

# Where a source's features come from, for the library's provenance record
source_location <- function(src) if (is.null(src$zip)) src$url else src$zip

# A source's features as sf in EPSG:4326: all of them, or just `units`.
# `offset` (degrees) is the generalization ArcGIS servers apply before sending;
# file sources are read as published.
read_source <- function(src, units = NULL, offset, batch = 40) {
      x <- if (is.null(src$zip)) {
            if (is.null(units)) {
                  arcgis_query(src, "1=1", offset)
            } else {
                  batches <- split(units, ceiling(seq_along(units) / batch))
                  do.call(rbind, lapply(batches, function(u) {
                        quoted <- if (isTRUE(src$numeric_key)) paste(u, collapse = ",")
                                  else paste0("'", u, "'", collapse = ",")
                        arcgis_query(src, sprintf("%s IN (%s)", src$key, quoted), offset)
                  }))
            }
      } else {
            x <- read_zip(src)
            if (is.null(units)) x else x[x[[src$key]] %in% units, ]
      }
      x[[src$key]] <- as.character(x[[src$key]])
      x
}

# Downloaded to a file first: the Forest Service server resets GDAL's own HTTP
# requests. Each query must fit in one response; a layer with more features
# than the server's record limit would need paging, which none of these do.
arcgis_query <- function(src, where, offset) {
      url <- paste0(src$url, "/query?",
                    "where=", utils::URLencode(where, reserved = TRUE),
                    "&outFields=", paste(src$fields, collapse = ","),
                    "&outSR=4326&maxAllowableOffset=", offset, "&f=geojson")
      tmp <- tempfile(fileext = ".geojson")
      on.exit(unlink(tmp))
      utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
      x <- sf::st_read(tmp, quiet = TRUE)
      count <- jsonlite::fromJSON(paste0(src$url, "/query?where=",
                                         utils::URLencode(where, reserved = TRUE),
                                         "&returnCountOnly=true&f=json"))$count
      if (nrow(x) != count) stop(src$url, ": got ", nrow(x), " of ", count, " features; needs paging")
      x
}

# Zipped shapefiles are downloaded once per R session
read_zip <- function(src) {
      zip <- file.path(tempdir(), basename(src$zip))
      if (!file.exists(zip)) utils::download.file(src$zip, zip, quiet = TRUE, mode = "wb")
      x <- sf::st_read(paste0("/vsizip/", zip), layer = src$layer, quiet = TRUE)
      sf::st_transform(x[, src$fields], 4326)
}
