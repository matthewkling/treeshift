## ANALOG IMPACT MODEL -----------------------------------------

### full raster -----------------------
# Single-entity surfaces for both timeframes (R/aim.R), from the result cache
# when these settings have been run before (R/cache.R). On a hit, clim() and
# pool() aren't computed at all. A community run's per-entity surfaces are
# multi_aim() (server/community.R).
aim <- reactive({
      settings <- lapply(setNames(nm = surfaces_setting_ids), function(id) model[[id]])
      ts_cached(result_cache, ts_surfaces_key(settings, core_data), function()
            ts_surfaces(clim(), pool(), model$fia_var, model$theta_clim, model$theta_geog))
})

### focal site ------------------------
site_analogs <- reactive({
      req(site$x, site$y)

      xy <- t(as.matrix(c(x = site$x, y = site$y), nrow = 1))
      # Extract from raw rasters to match pool() (consistent scaling).
      cp <- clim_pool()
      hst <- cbind(xy, terra::extract(cp$hst, xy, method = "bilinear"))
      fut <- cbind(xy, terra::extract(cp$fut, xy, method = "bilinear"))

      max_clim <- if(input$filter_clim) model$theta_clim * 3 else NULL

      # Hoist the pool reactive's result so the inner function doesn't
      # shadow it with its `index` argument.
      pool_data <- pool()
      # Climate columns in clim_hst and clim_fut share the same names
      # (e.g., "aet", "cwd") since the post-refactor rasters have no
      # 1/2 suffix. Tag them explicitly so downstream consumers
      # (clim_plot) can pull historical vs. future without ambiguity.
      clim_hst_tagged <- setNames(pool_data$clim_hst,
                                  paste0(names(pool_data$clim_hst), "_hst"))
      clim_fut_tagged <- setNames(pool_data$clim_fut,
                                  paste0(names(pool_data$clim_fut), "_fut"))
      pool_y <- cbind(pool_data$values, clim_hst_tagged, clim_fut_tagged)

      analog_query <- function(x, index){
            suppressWarnings(analogs::analog_search(
                  x = x,
                  pool = index,

                  geog = analogs::kernel(max = model$theta_geog * 3),
                  env  = analogs::kernel(max = max_clim),
                  select = "all",
                  stat = NULL,
                  y = pool_y
            ))
      }

      rev <- analog_query(fut, pool_data$index_hst)
      fwd <- analog_query(hst, pool_data$index_fut)
      bsl <- analog_query(hst, pool_data$index_hst)

      list(hst_site = hst,
           fut_site = fut,
           rev_analogs = rev,
           fwd_analogs = fwd,
           bsl_analogs = bsl
      )
})
