# Load TreeShift's modeling core (R/) for a script or test, outside the app:
#
#   source("scripts/load_core.R")
#   data <- ts_data("data")
#   res  <- ts_model(ts_settings(species = "Pinus ponderosa"), data)
#
# Run from the project root. The core is plain files sourced into the global
# environment rather than a package: a DESCRIPTION file at the root would make
# renv treat the project as a package and switch it to a different library.
for (.f in list.files("R", pattern = "^[^_].*[.]R$", full.names = TRUE)) {
      sys.source(.f, envir = globalenv(), keep.source = TRUE)
}
rm(.f)
