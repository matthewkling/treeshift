# The focal-polygon library: management-unit boundaries the user picks from, in
# place of uploads. Built by data_prep/boundaries_prep.R into
# data/boundaries.gpkg and distributed as a piggyback release asset, like the
# climate rasters.
#
# The release tag is the library version. A local copy whose recorded version
# differs (an older build) is replaced, so bumping the tag here is all it takes
# to roll a rebuilt library out to existing installs.
ensure_boundaries <- function(repo = "matthewkling/treeshift",
                              tag  = "boundaries-v1",
                              path = "data/boundaries.gpkg") {
      if (file.exists(path)) {
            version <- tryCatch(sf::st_read(path, layer = "about", quiet = TRUE)$version,
                                error = function(e) NA_character_)
            if (identical(version, tag)) return(invisible(FALSE))
      }

      message("Boundary library missing or out of date; downloading ", tag, " from ", repo, ".")
      if (!requireNamespace("piggyback", quietly = TRUE)) {
            stop("Package 'piggyback' is required to fetch the boundary library. ",
                 "Install it with install.packages('piggyback').")
      }
      piggyback::pb_download(basename(path), repo = repo, tag = tag,
                             dest = dirname(path), overwrite = TRUE)
      if (!file.exists(path)) {
            stop("Boundary library download failed: ", basename(path), " not found on release '",
                 tag, "' of ", repo, ".")
      }
      invisible(TRUE)
}

# Selectize choices for the library: unit names, grouped by collection in the
# order the library lists them (see boundaries_prep.R), with the stable unit id
# as the value.
boundary_choices <- function(units) {
      groups <- factor(units$collection, levels = unique(units$collection))
      lapply(split(units, groups), function(u) stats::setNames(u$id, u$name))
}
