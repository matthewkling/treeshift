library(raster)
library(hydro)
library(tidyverse)
library(terra)

bb <- extent(-125, -66, 24, 50)

derive_variables <- function(f, bb){{

      # water balance: compute from monthlies
      hst <- f %>%
            stack() %>%
            crop(bb) %>%
            hydro(already_latlong = TRUE, ncores = 8)
      layers <- tolower(names(hst))

      # mean annual temperature
      tmean <- f[grepl("_tas_", f)] %>%
            rast() %>% crop(bb) %>% mean() %>%
            setNames("tmean")

      # tmax of warmest month (bio5)
      tmaxwm <- f[grepl("_tasmax_", f)] %>%
            rast() %>% crop(bb) %>% max() %>%
            setNames("tmaxwm")

      # tmin of coldest month (bio6)
      tmincm <- f[grepl("_tasmin_", f)] %>%
            rast() %>% crop(bb) %>% min() %>%
            setNames("tmincm")

      c(tmean, tmaxwm, tmincm,
        setNames(rast(hst), layers))
}}

postprocess <- function(clim, outfile = "data/clim_quantized.tif", scales = NULL){

      # select focal variables
      vars <- c("tmean", "ppt", "aet", "cwd", "tmincm", "tmaxwm")
      clim <- clim[[vars]]

      # log-transform precip
      clim$ppt <- log10(clim$ppt)

      # compute, apply, and store scales
      if(is.null(scales)){
            scales <- list()
            for(var in vars){
                  x <- values(clim[[var]])

                  m <- mean(x, na.rm = T)
                  s <- sd(x, na.rm = T)

                  clim[[var]] <- (clim[[var]] - m) / s

                  scales[[var]]$mean <- m
                  scales[[var]]$sd <- s
            }
      } else {
            for(var in vars){
                  clim[[var]] <- (clim[[var]] - scales[[var]]$mean) / scales[[var]]$sd
            }
      }

      # mask out cells with 75+ % water
      m <- rast("/Volumes/T7/earthenv/landcover_DISCover_v1.0/consensus_full_class_12.tif") %>%
            crop(clim)
      m[m[] > 75] <- NA
      m[!is.na(m[])] <- 1
      clim <- mask(clim, m)

      # export
      clim <- writeRaster(
            clim,
            outfile,
            datatype = "INT2S",
            scale = 0.01,
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"),
            overwrite = TRUE
      )

      if(is.null(scales)) return(scales)
      return(invisible(outfile))
}


# historical -------------

f <- list.files("/Volumes/T7/CHELSA/v2/raw/", full.names = T)
f <- f[grepl("tas|pr", f)]
clim <- derive_variables(f, bb)
scales <- postprocess(clim, "data/clim_historical.tif")
saveRDS(scales, "data/clim_scales.rds")


# future --------------

for(era in c("2011_2040", "2041_2070", "2071_2100")){
      for(scenario in c("ssp126", "ssp370", "ssp585")){
            for(model in c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll")){
                  message(paste(era, scenario, model))

                  outfile <- paste0("data/clim_future/",
                                    scenario, "_",
                                    sub("_", "-", era), "_",
                                    model, ".tif")
                  if(file.exists(outfile)) next()

                  f <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
                  f <- f[grepl("tas|pr", f) &
                               grepl(era, f) &
                               grepl(scenario, f) &
                               grepl(model, f)]

                  clim <- derive_variables(f, bb)

                  postprocess(clim, outfile, scales = readRDS("data/clim_scales.rds"))
            }
      }
}


# climate deltas --------------

# data/clim_deltas.tif: future minus historical for every future file and
# climate variable, in native units (the stored z-scores times each variable's
# SD; ppt stays log10), masked to the contiguous US and aggregated 10x. One
# layer per <ssp>_<period>_<model>_<var>. The climate storylines rank models
# on it (ts_storyline_model()) and the Climate futures panel plots it.
# Self-contained: needs only data/, after data_prep/conus_prep.R.

scales <- readRDS("data/clim_scales.rds")
vars <- names(scales)
future_files <- list.files("data/clim_future", pattern = "\\.tif$", full.names = TRUE)

delta_files <- parallel::mclapply(future_files, function(f) {
      # terra objects are opened per worker; they don't survive forking
      conus <- rast("data/conus.tif")
      hst <- rast("data/clim_historical.tif")[[vars]]
      fut <- rast(f)[[vars]]
      delta <- (fut - hst) * sapply(vars, function(v) scales[[v]]$sd)
      delta <- aggregate(mask(delta, conus), fact = 10, fun = "mean", na.rm = TRUE)
      names(delta) <- paste(tools::file_path_sans_ext(basename(f)), vars, sep = "_")
      out <- tempfile(fileext = ".tif")
      writeRaster(delta, out, datatype = "FLT4S")
      out
}, mc.cores = 8)

writeRaster(rast(unlist(delta_files)), "data/clim_deltas.tif", datatype = "FLT4S",
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "INTERLEAVE=BAND"),
            overwrite = TRUE)
message("data/clim_deltas.tif: ", round(file.size("data/clim_deltas.tif") / 1024^2, 1), " MB")


# wind -------

wind <- rast("/Volumes/T7/CFSR/multidecadal_roses/wnd10m_annual_p1.tif") %>%
      rotate() %>%
      crop(clim)
writeRaster(
      wind,
      "data/wind_quantized.tif",
      datatype = "INT2S",
      scale = 0.01,
      gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"),
      overwrite = TRUE
)
