# Where each kind of focal-polygon library unit comes from: the agencies' own
# ArcGIS REST layers. Shared by boundaries_select.R and boundaries_prep.R.
#
# Each source names its layer, the field holding the agency's stable unit code
# (the library id is "<source>:<code>", e.g. "nps:YELL"), the fields to fetch,
# and how to map them onto the library's attributes. `collection` is the group a
# unit is listed under in the app's selector.

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
      )
)
