## CONNECTIVITY -------------------------------------------------

wind <- reactive({
      rast("data/wind_quantized.tif") %>%
            resample(aim()$bsl[[1]], method = "bilinear", threads = TRUE)
})

### Sources and targets ----------
electrodes <- reactive({
      if(model$cs_sources == "grid cells"){
            sources <- aim()$bsl$response
      }else if(model$cs_sources == "selected site" && !is.null(site$y)){
            xy <- matrix(c(site$x, site$y), nrow = 1)
            sources <- aim()$bsl$response
            sources <- rasterize(xy, sources) # * extract(sources, xy)[1,1]
      }else{
            # fia plots
            v <- plot_occ()[, c("x", "y")]
            v$z <- response_values()
            v <- v[v$z > 0,]
            v <- vect(v, geom = c("x", "y"))
            v <- rasterize(v, aim()$bsl$response, field = "z")
            sources <- v
      }
      sources[is.na(sources[])] <- 0
      sources <- sources / global(sources, "sum")$sum

      if(model$cs_targets == "grid cells"){
            destinations <- aim()$fut$response
      }else if(model$cs_targets == "selected site" && !is.null(site$y)){
            xy <- matrix(c(site$x, site$y), nrow = 1)
            destinations <- aim()$fut$response
            destinations <- rasterize(xy, destinations)
      }else{
            # fia plots
            v <- plot_occ()[, c("x", "y")]
            v$z <- response_values()
            v <- v[v$z > 0,]
            v <- vect(v, geom = c("x", "y"))
            v <- rasterize(v, aim()$fut$response, field = "z")
            destinations <- v
      }
      destinations[is.na(destinations[])] <- 0
      destinations <- destinations / global(destinations, "sum")$sum

      c(sources %>% setNames("sources"),
        destinations %>% setNames("destinations"))
})

suitability <- reactive({
      suit <- sqrt(aim()$bsl$response * aim()$fut$response)
      suit[is.na(suit)] <- 0
      suit
})

### Resistance --------------
resistance <- reactive({

      if(grepl("wind", model$cs_conductance_layer)){
            if(model$cs_sources == "grid cells"){
                  bearing <- gradient_bearing(electrodes()$destinations)
            }else{
                  # convert points to smooth surface for bearing computations
                  bearing <- electrodes()$destinations
                  bearing[bearing == 0] <- NA
                  bearing <- distance(bearing)
                  bearing <- gradient_bearing(-bearing)
            }
      }

      if(model$cs_conductance_layer == "wind"){
            conductance <- aligned_flux(bearing, wind())
            # note: consider adding a minimum conductance to account for local diffusion
      } else if(model$cs_conductance_layer == "suitability") {
            conductance <- suitability()
      } else if(model$cs_conductance_layer == "suitability * wind") {
            conductance <- aligned_flux(bearing, wind()) * suitability()
      } else { # uniform conductance
            # v <- 10^input$cs_conductance
            conductance <- setValues(electrodes()$sources, 1)
      }

      mean_cond <- 10 ^ model$cs_conductance_scale
      conductance <- conductance / global(conductance, "mean", na.rm = T)$mean * mean_cond
      conductance <- setNames(conductance, "conductance")
      resistance <- (1 / (conductance + 1e-3)) %>% setNames("resistance")
      c(conductance, resistance)
})

### Loss ---------------------
leakage <- reactive({
      loss <- switch(model$cs_loss_layer,
                     "none" = setValues(electrodes()$sources, 1e-6),
                     "uniform" = setValues(electrodes()$sources, 1),
                     "mean suitability" = (1 - sqrt(electrodes()$sources * electrodes()$destinations)),
                     "max suitability" = (1 - max(electrodes())),
                     "min suitability" = (1 - min(electrodes()))) %>%
            setNames("loss")
      mean_loss <- 10 ^ model$cs_loss_scale
      if(model$cs_loss_layer != "none") loss <- loss / global(loss, "mean", na.rm = T)$mean * mean_loss
      loss
})

### Circuitscape --------------
# From the result cache when these settings have been solved before
# (R/cache.R). The focal site is part of the key only when an electrode uses
# it, so moving the site doesn't discard an otherwise identical solve.
connectivity <- reactive({
      # Belt and braces: the stat option is withdrawn in community mode,
      # but a stale selection during a flush must not reach a solve.
      req(!multi_mode())

      settings <- lapply(setNames(nm = c(surfaces_setting_ids, connectivity_setting_ids)),
                         function(id) model[[id]])
      site_xy <- if ("selected site" %in% c(model$cs_sources, model$cs_targets) && !is.null(site$y))
            c(site$x, site$y)
      logged_run("connectivity", ts_connectivity_key(settings, core_data, site_xy), solve_connectivity,
                 describe = function(r) list(cells = terra::ncell(r)))
})

solve_connectivity <- function() {
      # First connectivity computation in this R process pays the Julia /
      # Circuitscape warmup; subsequent calls are a no-op. See the
      # ensure_cs_setup() definition in global setup.
      ensure_cs_setup()

      # Ground nodes
      grounds <- (electrodes()$destinations + leakage()) %>% setNames("grounds")

      # Run circuitscape
      fit_cs <- function(correction = 1e-12){
            cs_advanced(resistance()$resistance, electrodes()$sources,
                        grounds + correction, # add numerical stability for cholmod solver
                        resistance_is = "resistances", ground_is = "conductances",
                        solver = "cholmod") # cholmod is way faster than default
      }
      for(corr in seq(-12, -6)){ # larger factors are safer but introduce more spurious grounds, so start optimistic
            cs <- try(fit_cs(10 ^ corr))
            if(class(cs) != "try-error") break
      }

      # Derived variables
      # cs_advanced() returns a SpatRaster with layers 'current' and
      # 'voltage'. (Advanced mode solves a single circuit, so there is
      # no 'cumulative_current' layer -- that name applies only to the
      # multi-focal-node modes.)
      current <- cs$current %>% setNames("current_flow")
      voltage <- cs$voltage
      curr_ground <- grounds %>% "*"(voltage) %>% setNames("current_ground")
      curr_dest <- electrodes()$destinations %>% "*"(voltage) %>% setNames("current_destinations")
      curr_loss <- leakage() %>% "*"(voltage) %>% setNames("current_loss")

      # Power dissipation (pinch points)
      curr_power <- (current^2 * resistance()$resistance) %>% setNames("current_power")
      curr_dir <- gradient_bearing(-voltage) %>% setNames("current_direction")

      # Return results
      c(current, voltage, curr_dest, curr_loss, curr_power,
        resistance()$conductance, leakage(), curr_dir) %>%
            mask(aim()$bsl$response)
}
