# The focal-polygon library: management-unit boundaries the user picks from, in
# place of uploads. Built by data_prep/boundaries_prep.R into
# data/boundaries.gpkg and distributed as a piggyback release asset, like the
# climate rasters.
#
# The release tag is the library version. A local copy whose recorded version
# differs (an older build) is replaced, so bumping the tag here is all it takes
# to roll a rebuilt library out to existing installs.
ensure_boundaries <- function(repo = "matthewkling/treeshift",
                              tag  = "boundaries-v1.2",
                              path = "data/boundaries.gpkg") {
      version <- function() {
            if (!file.exists(path)) return(NA_character_)
            tryCatch(sf::st_read(path, layer = "about", quiet = TRUE)$version,
                     error = function(e) NA_character_)
      }
      if (identical(version(), tag)) return(invisible(FALSE))

      message("Boundary library missing or out of date; downloading ", tag, " from ", repo, ".")
      if (!requireNamespace("piggyback", quietly = TRUE)) {
            stop("Package 'piggyback' is required to fetch the boundary library. ",
                 "Install it with install.packages('piggyback').")
      }
      # By default piggyback skips the download when the local file is newer
      # than the release asset, which would keep an outdated library copied
      # in after the release was made
      piggyback::pb_download(basename(path), repo = repo, tag = tag,
                             dest = dirname(path), overwrite = TRUE, use_timestamps = FALSE)
      if (!identical(version(), tag)) {
            stop("Boundary library download failed: ", basename(path), " for release '",
                 tag, "' of ", repo, " is missing or has the wrong version.")
      }
      invisible(TRUE)
}

# The selector picks a collection (the kind of unit) first, then a unit within
# it. Collections are in the order the library lists them (see
# boundaries_prep.R).
boundary_collections <- function(units) unique(units$collection)

# The collection a unit id belongs to; the first collection for no unit ("")
boundary_collection_of <- function(units, id) {
      if (nzchar(id) && id %in% units$id) units$collection[units$id == id] else boundary_collections(units)[1]
}

# Unit choices within one collection: names, with the stable unit id as the
# value, after a "None" entry
boundary_choices <- function(units, collection) {
      u <- units[units$collection == collection, ]
      c("None" = "", stats::setNames(u$id, u$name))
}
