# Put the default settings' results in the result cache (R/cache.R), so the app
# opens on them without computing.
#
#   Rscript scripts/precompute_defaults.R
#
# Run from the project root, after each deploy, as the user the app runs as (or
# with TREESHIFT_CACHE_DIR set to the app's cache directory). The defaults come
# from ts_default_settings() (R/settings.R), the same source as the app's, so
# changing them needs nothing here. A change to the defaults, the code, the data
# or the packages changes the cache key; until this is rerun, the app computes
# the new defaults on first load and caches them itself.

source("scripts/load_core.R")

data <- ts_data("data")
settings <- ts_settings()
key <- ts_surfaces_key(settings, data)
if (is.null(key)) stop("Can't fingerprint the default run: an input file is missing.")

cache <- ts_cache()
fresh <- !cache$exists(key)
t <- system.time(ts_cached(cache, key, function() ts_model(settings, data)$surfaces))
message(if (fresh) sprintf("Cached the default surfaces in %.1f s: ", t[["elapsed"]])
        else "The default surfaces were already cached: ", key)
