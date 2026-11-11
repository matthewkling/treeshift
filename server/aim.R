## ANALOG IMPACT MODEL -----------------------------------------

### full raster -----------------------
aim <- reactive({

      is_ratio <- !is.null(ratio_pieces())

      # Stats that describe the WEIGHTS rather than an outcome column.
      # analog_impact() returns one layer for each of these regardless of
      # how many outcome columns `y` carries. "count" and "sum_weights"
      # were requested here but read nowhere in the app; they are dropped.
      # (sum_weights becomes useful again if the multi-species module
      # wants a per-species support signal: the weighted count of
      # contributing presence plots is weighted_mean_<sp> * sum_weights.)
      weight_stats <- "ess"

      # Outcome columns for AIM. pool()$values also carries
      # samp_method_cd, which exists for the site_analogs() / FIA-marker
      # path and is not a response -- interpolating it produced an extra
      # unread layer on every run.
      y_cols <- setdiff(names(pool()$values), "samp_method_cd")

      # analog_impact() names outcome layers `weighted_mean_<column>` when
      # `y` has two or more columns, but bare `weighted_mean` when it has
      # exactly one. Stripping the prefix with gsub() therefore yields
      # "weighted_mean" rather than "response" in the single-column case,
      # leaving aim()$bsl$response NULL. That never fired only because
      # samp_method_cd rode along as a spurious second column -- i.e. the
      # rename was load-bearing on an accident. Rename positionally
      # instead: outcome layers are whatever isn't a weight stat, returned
      # in y-column order.
      name_outcomes <- function(r){
            is_outcome <- !names(r) %in% weight_stats
            stopifnot(sum(is_outcome) == length(y_cols))
            names(r)[is_outcome] <- y_cols
            r
      }

      surf <- aim_surfaces(list(bsl = clim()$hst, fut = clim()$fut), pool(),
                           y = pool()$values[, y_cols, drop = FALSE],
                           theta_clim = model$theta_clim, theta_geog = model$theta_geog,
                           stat = c("weighted_mean", weight_stats))
      aim_bsl <- name_outcomes(surf$bsl)
      aim_fut <- name_outcomes(surf$fut)

      # Pooled-ratio surface: ratio-of-means from the weighted_mean layers
      # occ/num/denom, combined into the standard `response` layer so all
      # downstream code is unchanged.
      #
      # num/denom is already conditional on occurrence -- absent plots
      # contribute 0 to both and drop out -- so it reads "the ratio among
      # plots where the focal entity occurs".
      #
      # The regeneration fractions stay conditional and are MASKED (NA)
      # where no occupied plot falls within kernel range. Multiplying by
      # occurrence instead would suppress unsuitable ground by driving the
      # value to 0 -- but 0 already means "occurs, adults only" (relict), so
      # the product would conflate two very different states. Masking keeps
      # them distinct: blank = doesn't occur, 0 = occurs but no regen. This
      # also matches the per-plot markers, which are already NA where the
      # species is absent and 0 where it's adults-only.
      #
      # Proportion basal area keeps the occurrence weighting: there 0 means
      # "the focal species is 0% of a real stand" vs "there is no stand" --
      # both read as "the species isn't here", so collapsing them costs
      # little, and the unconditional surface is better behaved.
      if (is_ratio) {
            conditional <- grepl("fraction", model$fia_var)
            combine_ratio <- function(r){
                  resp <- r[["num"]] / r[["denom"]]
                  if (!conditional) resp <- r[["occ"]] * resp
                  # denom == 0: no occupied plot in range. Undefined for the
                  # fractions; 0 for the occurrence-weighted proportion.
                  resp[!is.finite(resp[])] <- if (conditional) NA else 0
                  names(resp) <- "response"
                  c(r, resp)
            }
            aim_bsl <- combine_ratio(aim_bsl)
            aim_fut <- combine_ratio(aim_fut)
      }

      list(fut = aim_fut,
           bsl = aim_bsl)
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
