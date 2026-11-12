## MODEL EVALUATION ----------------------------

## What the evaluation cross-validates.
##
## In single-entity mode this is just the focal entity's response. In
## community mode it is the plot-level counterpart of the selected
## community measure -- observed richness at each plot, for instance --
## so the module scores the surface the map is actually showing rather
## than a species the user cannot even see.
##
## Only linear measures can be handled this way, and richness is the one
## that qualifies. Because AIM is a weighted mean with entity-independent
## weights, reducing then interpolating equals interpolating then
## reducing, so the held-out prediction of plot richness IS the held-out
## Sum(p_s). A non-linear measure would need each entity predicted
## separately at each held-out plot, and analog_cv() scores one response
## at a time -- 143 CV runs, not the cheap multi-outcome trick.
eval_response <- reactive({
      if (!multi_mode()) return(response_values())
      req(identical(input$comm_measure, "richness"))
      req(isTruthy(model$fia_var), model$fia_var %in% multi_vars)
      comm_reduce_rows(multi_values(), "richness")
})

## Measures the evaluation cannot score, and why. Returned as a message so
## the panel explains itself rather than rendering empty.
eval_unsupported <- reactive({
      if (!multi_mode()) return(NULL)
      if (identical(input$comm_measure, "richness")) return(NULL)
      "Accuracy can only be measured for richness. Switch the community measure to richness, or switch back to a single species."
})

model_eval <- reactiveValues(comparison = NULL,
                             scatter = NULL,
                             stamp = NULL,
                             swept = NULL)

# Message shown when a scope is selected whose defining input isn't set.
# The scope list is static rather than rebuilt as inputs appear: a
# reactive choice list makes setting a focal site rebuild the selectize
# control, whose round-trip re-emits input$eval_scope and re-renders this
# section even under a scope that ignores the site. Stating the
# precondition is also more useful than hiding the option.
eval_scope_missing_msg <- function(scope){
      switch(scope,
             "polygon"  = "Choose a focal polygon under 'Focal polygon and site' in the left sidebar to evaluate over it.",
             "buffered" = "Choose a focal polygon under 'Focal polygon and site' in the left sidebar to evaluate over its buffered neighborhood.",
             "site"     = "Click the map to set a focal site, then this scope evaluates the plots within its geographic cutoff.",
             "This evaluation scope is not available.")
}

# Stamp components. Plain functions rather than reactives, so a caller
# takes a dependency on the polygon or the site only when it actually
# calls for that piece. Bundling them into one reactive made every map
# click re-render the heatmap and its warning banner, even under a scope
# that ignores the site entirely.
eval_polygon_id <- function(){
      id <- polygon_active()
      if (isTruthy(id)) id else NA_character_
}
eval_site_id <- function(){
      if (is.null(site$x) || is.null(site$y)) NA_character_
      else sprintf("%.6f,%.6f", site$x, site$y)
}

# Settings that define the response and the analog pool, for every scope.
eval_model_stamp <- function(){
      list(track      = model$track,
           focal      = focal_id(),
           # focal_id() already distinguishes community mode from a single
           # entity, but not which measure is being scored.
           measure    = if (multi_mode()) input$comm_measure else NA_character_,
           fia_var    = model$fia_var,
           bbox       = c(model$bbox_x, model$bbox_y),
           clim_vars  = paste(sort(model$clim_vars), collapse = "_"),
           theta_clim = model$theta_clim,
           theta_geog = model$theta_geog)
}

# Full stamp, recorded when a run completes.
eval_run_stamp <- function(){
      c(eval_model_stamp(),
        list(polygon = eval_polygon_id(), site = eval_site_id()))
}

# Stamped fields that have changed since the stored run, with two
# exemptions applied.
#
# SWEPT DIMENSIONS. A bandwidth or variable-set sweep varies those inputs
# by design, and "Apply optimal parameters" then moves them deliberately.
# Counting them as staleness would make the heatmap mark itself stale the
# instant the user acts on its own recommendation.
#
# SCOPE RELEVANCE. Handled by only assembling the fields the current scope
# depends on: stamp_diff() intersects names, so a field left out of `now`
# is simply not compared. theta_geog is the one field with two roles -- it
# is a swept kernel parameter, but it also sets the buffer and neighborhood
# radii, so under those scopes it is never exempt even when swept.
eval_stale_fields <- reactive({
      req(model_eval$comparison)
      scope <- if (isTruthy(input$eval_scope)) input$eval_scope else "domain"

      exempt <- if (is.null(model_eval$swept)) character() else model_eval$swept
      if (scope %in% c("buffered", "site")) exempt <- setdiff(exempt, "theta_geog")

      now <- eval_model_stamp()
      if (scope %in% c("polygon", "buffered")) now$polygon <- eval_polygon_id()
      if (scope == "site") now$site <- eval_site_id()

      stamp_diff(now, model_eval$stamp, exempt)
})

# Logical membership vector per scope, for a set of pool plot coordinates.
# `xy` must be the POST-filter coordinates so the vectors align with the
# rows analog_cv() returns.
#
# `which` limits which scopes are computed, and that is load-bearing rather
# than an optimization: reading site$x or the polygon here takes a reactive
# dependency on it, so computing every scope unconditionally makes a map
# click invalidate the scatterplot even when the domain scope is selected.
# The `&&` short-circuits before touching those reactives.
eval_scope_membership <- function(xy, which = unname(eval_scope_names)){
      keep <- list()
      if ("domain" %in% which) keep$domain <- rep(TRUE, nrow(xy))

      # Check polygon_active() directly rather than calling user_polygon(),
      # whose internal req() would silently cancel the whole evaluation
      # when no polygon is loaded. Same guard as polygon_density().
      if (any(c("polygon", "buffered") %in% which) &&
          isTruthy(polygon_active())) {
            pts <- st_as_sf(as.data.frame(xy), coords = c("x", "y"), crs = 4326)
            if ("polygon" %in% which)
                  keep$polygon <- lengths(st_intersects(pts, st_union(user_polygon()))) > 0
            if ("buffered" %in% which)
                  keep$buffered <- lengths(st_intersects(pts, user_polygon_buffered())) > 0
      }
      if ("site" %in% which && !is.null(site$x) && !is.null(site$y)) {
            dkm <- geosphere::distHaversine(cbind(site$x, site$y),
                                            as.matrix(xy[, c("x", "y")])) / 1000
            keep$site <- dkm <= model$theta_geog * 3
      }
      keep
}


eval_var_set <- function(clim_vars, vary_bandwidths = TRUE){
      # Extract historical climate from the raw raster at each pool plot,
      # matching pool().
      clim_r <- rast("data/clim_historical.tif")[[clim_vars]]

      # Build pool (x/y + climate) and response from a SINGLE plot_occ()
      # snapshot so the two are row-aligned. Capture plot_occ() once and
      # derive the response from that same object rather than re-deriving
      # it, which could in principle see a different value.
      po <- plot_occ()
      xy <- as.data.frame(po[, c("x", "y")])
      clim_hst <- terra::extract(clim_r, xy, ID = FALSE, method = "bilinear")
      d <- cbind(xy, clim_hst)
      values <- eval_response()
      req(!is.null(values), length(values) == nrow(po))

      # analog_cv() coerces `pool` internally and requires `y` to have
      # exactly nrow(pool) rows. bilinear extraction yields NA climate
      # for plots off the raster's edge; drop those rows from BOTH pool
      # and response together so the row counts stay locked. (The old
      # two-fold path tolerated mismatches because build_analog_index()
      # never row-checked y; analog_cv() checks immediately.)
      complete <- stats::complete.cases(d) & !is.na(values)
      d <- d[complete, , drop = FALSE]
      values <- values[complete]
      stopifnot(nrow(d) == length(values))

      # The CV pool in AIM coordinates, so held-out predictions use the
      # same distances as the maps. `d` stays lon/lat for scopes and plots.
      d_aim <- data.frame(aim_xy(d), d[, clim_vars, drop = FALSE])

      # Scope membership, aligned with the filtered rows -- and therefore
      # with the rows analog_cv() returns. Computed once per variable set
      # rather than per bandwidth cell, since it depends on neither.
      #
      # Only on the summarizing path. eval() is called with summarize=TRUE
      # exactly when vary_bandwidths is TRUE, and the scatter path resolves
      # its own scope in the renderer. Computing this unconditionally reads
      # site$x and the polygon, which would make eval_scatter_data() depend
      # on them -- and a single map click would then re-run the whole LOO
      # cross-validation to produce an identical plot.
      scopes <- if (vary_bandwidths) eval_scope_membership(d[, c("x", "y")]) else NULL

      # Single-configuration evaluation: leave-one-out CV via the
      # analogs package's native analog_cv(), with performance metrics
      # from cv_performance(). analog_cv() predicts each plot from a
      # neighborhood excluding its own pool row (exclude_self, handled
      # internally); cv_performance() then computes held-out metrics
      # from the obs/residual columns.
      #
      # Note on cost: analog_cv() rebuilds the analog index internally
      # on every call (it rejects pre-built indices by design), so the
      # bandwidth grid below incurs one index build per cell. A direct
      # analog_impact() call with a hoisted index cannot avoid this --
      # exclude_self requires raw pool data, which forces an internal
      # rebuild regardless. analog_cv() is retained here for clarity and
      # correctness; the per-cell rebuild is an accepted cost.
      #
      # cv_performance() auto-classifies the response as binary
      # (obs all in {0,1}, both classes present) or continuous and
      # returns the corresponding metric set in long format
      # (columns: variable, type, metric, value):
      #   continuous -> n, rmse, mae, bias, r2
      #   binary     -> n, auc, tss, tss_threshold, brier
      # `r2` is out-of-sample R-squared (1 - SS_res/SS_tot), which can
      # be negative -- unlike the old cor(pred,obs)^2. There is no AIC
      # equivalent: AIC is an in-sample criterion and cv_performance()
      # is a held-out CV tool, so `delta AIC` was dropped from the UI.
      eval <- function(theta_clim, theta_geog, summarize = TRUE){
            k <- aim_kernels(theta_clim, theta_geog)
            cv <- suppressWarnings(analogs::analog_cv(
                  fun        = analogs::analog_impact,
                  pool       = d_aim,
                  y          = values,
                  cv_method  = "loo",
                  env        = k$env,
                  geog       = k$geog,
                  stat       = "weighted_mean",
                  coord_type = "projected"))

            if(!summarize){
                  # Scatterplot path: return the cv object itself plus the
                  # plot coordinates. Scoping is deferred to the renderer
                  # for the same reason it is here -- membership depends on
                  # the polygon and the focal site, and resolving it inside
                  # this function would make every map click re-run a full
                  # LOO cross-validation. Returning cv (rather than just
                  # obs/pred) lets the renderer call cv_performance() on a
                  # row subset, exactly as the heatmap path does.
                  return(list(cv = cv, xy = d[, c("x", "y"), drop = FALSE]))
            }

            # Heatmap path: native CV performance metrics, computed once
            # per scope from the SAME cv object. Scoping is a row subset,
            # so every scope's metrics come out of one run for the cost of
            # an extra aggregation -- which is why the scope selector can
            # be a display control rather than something that recomputes.
            #
            # Rename `metric` -> `stat` so downstream consumers
            # (eval_heatmap, optimal_params) keep their existing column
            # vocabulary, and tag with the bandwidth coordinates.
            purrr::imap(scopes, function(keep, scope_name){
                  if (sum(keep) < eval_min_n) return(NULL)
                  # A small subset can hold a single class, leaving AUC and
                  # TSS undefined. Drop that scope's cell rather than
                  # failing the whole run.
                  p <- tryCatch(analogs::cv_performance(cv[keep, , drop = FALSE]),
                                error = function(e) NULL)
                  if (is.null(p)) return(NULL)
                  p %>%
                        dplyr::rename(stat = metric) %>%
                        dplyr::select(stat, value) %>%
                        dplyr::mutate(scope = scope_name,
                                      theta_clim = theta_clim,
                                      theta_geog = theta_geog)
            }) %>%
                  purrr::list_rbind()
      }

      if(vary_bandwidths){
            mult <- c(.25, .5, 1, 2, 4)
            clim_mult <- if("clim" %in% input$eval_params) mult else 1
            geog_mult <- if("geog" %in% input$eval_params) mult else 1
            params <- tidyr::expand_grid(theta_clim = model$theta_clim * clim_mult,
                                         theta_geog = model$theta_geog * geog_mult)

            e <- purrr::pmap(params, eval) %>%
                  purrr::list_rbind() %>%
                  mutate(clim_vars = paste(clim_vars, collapse = "_"))
      }else{
            e <- eval(model$theta_clim, model$theta_geog,
                      summarize = FALSE)
      }

      return(e)
}

observeEvent(input$run_eval, {
      # Refuse rather than sweeping bandwidths against a response the
      # module cannot score.
      if (!is.null(eval_unsupported())) {
            showNotification(eval_unsupported(), type = "warning", duration = 10)
            return()
      }

      if("vars" %in% input$eval_params){
            cv <- model$clim_vars
            n <- min(input$eval_n_vars[1], length(cv)):min(input$eval_n_vars[2], length(cv))
            v <- lapply(n, function(x) apply(combn(cv, x), 2, list))
            v <- unlist(unlist(v, recursive = F), recursive = F)
      }else{
            v <- list(model$clim_vars)
      }

      e <- v %>% lapply(eval_var_set) %>% bind_rows()

      model_eval$comparison <- e
      # Snapshot the settings these results describe, plus the dimensions
      # this run deliberately varied. Both are needed: the stamp detects
      # change, the swept list keeps the run from indicting itself for the
      # very dimensions it was exploring.
      model_eval$stamp <- eval_run_stamp()
      model_eval$swept <- c(
            if ("clim" %in% input$eval_params) "theta_clim",
            if ("geog" %in% input$eval_params) "theta_geog",
            if ("vars" %in% input$eval_params) "clim_vars"
      )
})

## Staleness marker. Per-scope metrics are aggregated at run time, so they
## cannot be recomputed for new settings without re-running the CV.
## Storing per-plot residuals instead would make rescoping free, but the
## eval grid can reach hundreds of configurations across ~88k pool plots --
## tens of millions of rows. Pre-aggregation is the right call; explicit
## invalidation is its cost.
output$eval_stale <- renderUI({
      stale_banner(eval_stale_fields(), eval_stamp_labels,
                   action = "Click 'Run evaluations' to update them.")
})

observeEvent(input$optimal_params, {
      req(model_eval$comparison)

      # Refuse rather than silently applying bandwidths derived from
      # settings the user has since moved away from. The warning banner is
      # visible, but a button that ignores its own panel's warning is worse
      # than one that explains itself.
      stale <- eval_stale_fields()
      if (length(stale)) {
            showNotification(
                  sprintf("These results predate a change to the %s. Re-run evaluations before applying optimal parameters.",
                          stale_phrase(stale, eval_stamp_labels)),
                  type = "warning", duration = 8)
            return()
      }

      metric <- switch(input$eval_stat,
                       "RMSE" = "rmse",
                       "R-squared" = "r2",
                       "AUC" = "auc",
                       "TSS" = "tss",
                       "Brier" = "brier")

      # Metrics where larger is better; the rest are minimized.
      higher_better <- c("r2", "auc", "tss")

      # Optimize against the SELECTED scope: the best bandwidth for a
      # particular unit is not generally the best bandwidth domain-wide,
      # and the whole point of the scope selector is to let the user ask
      # for the first one.
      opt <- model_eval$comparison %>% ungroup() %>%
            filter(scope == input$eval_scope, stat == metric) %>%
            mutate(value = ifelse(stat %in% higher_better, value, -value)) %>%
            filter(value == max(value))

      shiny::validate(need(nrow(opt) > 0,
                           "No results for this area. Re-run the evaluations, or choose a larger area."))

      # Applied as well as set: choosing optimal parameters asks for them on
      # the map, so the Run model button isn't left to press
      best <- list(theta_geog = opt$theta_geog, theta_clim = opt$theta_clim,
                   clim_vars = unlist(strsplit(opt$clim_vars, "_")))
      updateSliderInput(session, "theta_geog", value = best$theta_geog)
      updateSliderInput(session, "theta_clim", value = best$theta_clim)
      updateSelectizeInput(session, "clim_vars", selected = best$clim_vars)
      apply_model_settings(best)
})

output$eval_heatmap <- renderPlot({
      shiny::validate(need(is.null(eval_unsupported()), eval_unsupported()))
      req(model_eval$comparison)

      metric <- switch(input$eval_stat,
                       "RMSE" = "rmse",
                       "R-squared" = "r2",
                       "AUC" = "auc",
                       "TSS" = "tss",
                       "Brier" = "brier")

      # Metrics where larger is better; the rest are minimized.
      opt <- if(metric %in% c("r2", "auc", "tss")) 1 else -1

      scope_label <- names(eval_scope_names)[eval_scope_names == input$eval_scope]

      e <- model_eval$comparison %>% ungroup() %>%
            filter(scope == input$eval_scope)

      # No rows means one of three things: the scope's input isn't set,
      # fewer than eval_min_n plots fall inside it, or the last run
      # predates the selection. Name the first case specifically, since it
      # is the only one the user can act on immediately.
      scope_set <- switch(input$eval_scope,
                          "polygon"  = ,
                          "buffered" = isTruthy(polygon_active()),
                          "site"     = !is.null(site$x) && !is.null(site$y),
                          TRUE)
      shiny::validate(need(
            nrow(e) > 0,
            if (!scope_set) eval_scope_missing_msg(input$eval_scope)
            else if (input$eval_scope %in% c("polygon", "buffered") &&
                     !polygon_in_domain(input$eval_scope)) polygon_outside_msg
            else sprintf("No results for '%s'. Either the last run was done before you chose this area, or fewer than %d plots fall inside it.",
                         scope_label, eval_min_n)))

      # cv_performance() reports n alongside the requested metric. It can
      # vary across bandwidth cells, since a plot with no analog in range
      # yields no held-out prediction, so report the range rather than
      # implying a single number.
      n_vals <- e %>% filter(stat == "n") %>% pull(value)
      n_text <- if (!length(n_vals)) "n unknown"
      else if (min(n_vals) == max(n_vals)) sprintf("n = %d", as.integer(min(n_vals)))
      else sprintf("n = %d\u2013%d", as.integer(min(n_vals)), as.integer(max(n_vals)))

      e <- e %>%
            filter(stat == metric) %>%
            mutate(label = signif(value, 2),
                   optimal = value * opt == max(value * opt),
                   text_color = case_when(optimal ~ "red",
                                          scales::rescale(value) < .5 ~ "white",
                                          TRUE ~ "black"),
                   clim_vars = gsub("_", " + ", clim_vars),
                   n_vars = stringr::str_count(clim_vars, "\\+") + 1)

      shiny::validate(need(nrow(e) > 0,
                           sprintf("'%s' cannot be calculated for this variable within '%s'.",
                                   input$eval_stat, scope_label)))

      subtitle <- paste0(scope_label, " \u2014 ", n_text,
                         if (length(n_vals) && min(n_vals) < 50)
                               "  (small samples; treat these numbers as indicative)" else "")

      p <- ggplot(e, aes(theta_geog, theta_clim, fill = value, label = label)) +
            facet_wrap(~clim_vars) +
            geom_tile() +
            geom_text(aes(color = I(text_color), fontface = "bold")) +
            scale_x_log10(breaks = unique(e$theta_geog)) +
            scale_y_log10(breaks = unique(e$theta_clim)) +
            scale_fill_viridis_c() +
            guides(fill = guide_colorbar(barwidth = 12)) +
            theme(legend.position = "bottom") +
            strip_style +
            labs(y = "climate bandwidth (z)",
                 x = "geographic bandwidth (km)",
                 subtitle = subtitle,
                 fill = paste0(input$eval_stat, " (", eval_response_label(), ")"))

      # `%>%` binds tighter than `+`, so this cannot be chained onto the
      # ggplot expression above -- it would pipe only the labs() object.
      gg_stale(p, length(eval_stale_fields()) > 0)
})


eval_scatter_data <- reactive({
      eval_var_set(model$clim_vars, vary_bandwidths = FALSE)
})

# Scope resolution for the single-configuration section. Shared by the
# scatter and its stat readout so the two can never disagree about which
# plots they are describing.
eval_scatter_scope <- reactive({
      e <- eval_scatter_data()
      req(e)
      scope <- if (isTruthy(input$eval_scope)) input$eval_scope else "domain"
      keep <- eval_scope_membership(e$xy, which = scope)[[scope]]
      # NULL means the scope's defining input isn't set. Falling back to
      # every plot would label domain-wide results with the selected
      # scope's name, which is worse than showing nothing.
      shiny::validate(need(!is.null(keep), eval_scope_missing_msg(scope)))
      if (scope %in% c("polygon", "buffered")) {
            shiny::validate(need(polygon_in_domain(scope), polygon_outside_msg))
      }
      list(scope = scope,
           label = names(eval_scope_names)[eval_scope_names == scope],
           keep  = keep)
})

output$eval_scatter <- renderPlot({
      shiny::validate(need(is.null(eval_unsupported()), eval_unsupported()))
      e <- eval_scatter_data()
      s <- eval_scatter_scope()
      req(e, s)

      d <- tibble::tibble(obs      = e$cv$obs,
                          pred     = e$cv$obs - e$cv$residual,
                          in_scope = s$keep) %>%
            filter(!is.na(obs), !is.na(pred))

      shiny::validate(need(
            sum(d$in_scope) >= eval_min_n,
            sprintf("Only %d plot(s) fall within '%s' \u2014 too few to judge accuracy from.",
                    sum(d$in_scope), s$label)))

      # One calibration plot for every response type: bin the predictions,
      # and for each bin show what was actually observed.
      #
      # Predicted on the x axis, not observed, because that is the direction
      # the tool is used in. "Given a plot with 7 species, what does the
      # model say" is the inverse of anyone's question; "the map reads 3.2
      # here, what is really there" is the question, and this answers it
      # directly. It also reframes shrinkage as a usable correction rather
      # than an alarm: predictions that plateau at 3.2 against observations
      # running to 10 read as "when it says 3.2, expect about 5".
      #
      # Equal-COUNT bins, so every point carries the same standard error and
      # none is more trustworthy than another. That is also why the points
      # are not sized: sizing by bin count stretched a range of 1,142 to
      # 1,143 plots across the full radius scale, turning rounding noise
      # into apparently dramatic differences. The rug shows where the
      # predictions actually fall, which is the real density question.
      discrete_obs <- all(abs(d$obs - round(d$obs)) < 1e-8, na.rm = TRUE)
      binary_obs   <- all(d$obs %in% c(0, 1))

      # Equal-count bins via explicit quantile breaks, de-duplicated.
      # cut_number() fails outright when quantiles tie, and many responses
      # have a large spike of identical predictions -- proportion basal area
      # is zero wherever the species is absent. Falling back to equal-width
      # there produced one bin of 77,000 plots and others of 8, which is
      # what made the tail of the curve fly around.
      #
      # Merging tied breaks instead gives one wide bin covering the spike
      # and equal-count bins across the rest: counts still vary, but by a
      # factor rather than by four orders of magnitude.
      calib <- function(x) {
            nb <- max(3, min(20, floor(nrow(x) / 200)))
            qs <- unique(stats::quantile(x$pred, seq(0, 1, length.out = nb + 1),
                                         na.rm = TRUE))
            x$bin <- if (length(qs) >= 4) cut(x$pred, qs, include.lowest = TRUE)
            else cut(x$pred, min(nb, 5))
            x %>%
                  group_by(bin) %>%
                  summarize(p  = mean(pred),
                            o  = mean(obs),
                            lo = stats::quantile(obs, 0.25),
                            hi = stats::quantile(obs, 0.75),
                            n  = dplyr::n(),
                            .groups = "drop") %>%
                  # A mean of a handful of plots is noise dressed as a
                  # calibration point.
                  filter(n >= 30)
      }

      f_in  <- d[d$in_scope, , drop = FALSE]
      f_out <- d[!d$in_scope, , drop = FALSE]
      b <- calib(f_in)

      shiny::validate(need(nrow(b) >= 3,
                           "Not enough spread in the predictions to assess calibration."))

      # Equal-count binning can fail when many predictions tie -- a
      # presence/absence response has thousands at effectively zero -- in
      # which case the fallback is equal-width and the counts vary widely.
      # Report which happened rather than claiming a per-point sample size
      # that is not true.
      n_per_bin <- if (diff(range(b$n)) <= 0.02 * stats::median(b$n))
            sprintf("~%s plots", format(round(stats::median(b$n)), big.mark = ","))
      else sprintf("%s-%s plots",
                   format(min(b$n), big.mark = ","),
                   format(max(b$n), big.mark = ","))

      p <- ggplot(b, aes(p, o)) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                        color = "dodgerblue")

      # Raw points for continuous responses only. Their observed values are
      # spread out, so the cloud shows the joint distribution the bins
      # summarize; a discrete response would just stack them into stripes.
      if (!discrete_obs) p <- p +
            geom_point(data = f_in, aes(pred, obs), inherit.aes = FALSE,
                       color = "gray55",
                       # Floor at 0.12, not 0.03: the low floor was carried
                       # over from a version where points stacked on top of
                       # each other and reinforced. Here they are spread
                       # across a cloud, so each one has to be visible on
                       # its own.
                       alpha = max(0.12, min(0.5, 8000 / max(nrow(f_in), 1))),
                       size = 0.8)

      p <- p +
            geom_rug(data = f_in, aes(x = pred), inherit.aes = FALSE,
                     alpha = max(0.01, min(0.3, 2000 / nrow(f_in))),
                     sides = "b", color = "gray30")

      # When a scope narrows the evaluation, the same curve for everything
      # OUTSIDE it goes underneath in gray. That answers the question the
      # scope selector exists to raise -- is the model as well calibrated
      # here as it is elsewhere -- which a single curve cannot.
      if (nrow(f_out) >= 200) {
            b_out <- calib(f_out)
            if (nrow(b_out) >= 3) p <- p +
                        geom_line(data = b_out, color = "gray65", linewidth = 0.7)
      }

      # Interquartile spread of what was observed in each bin, as a band
      # rather than bars on the points. The point is the MEAN, which for a
      # skewed response legitimately falls outside the quartiles -- drawn as
      # bars attached to each point that reads as an error bar that has come
      # adrift. A band reads as what it is: the spread of outcomes behind
      # the curve, a separate quantity from the calibration estimate.
      #
      # Skipped for presence/absence, where the quartiles of a 0/1 variable
      # carry nothing the mean does not.
      if (!binary_obs) p <- p +
            geom_ribbon(aes(ymin = lo, ymax = hi), fill = "red", alpha = 0.12)

      p <- p +
            geom_line(color = "red", linewidth = 0.9) +
            geom_point(color = "red", size = 2.2)

      # sqrt axes tame the skew of abundance variables; on counts of 1-12 or
      # a probability they buy nothing.
      if (!discrete_obs) p <- p + scale_x_sqrt() + scale_y_sqrt()

      p +
            coord_fixed() +
            theme_bw() +
            labs(x = paste0("predicted ", eval_response_label()),
                 y = if (binary_obs) "observed frequency"
                 else paste0("observed ", eval_response_label()),
                 title = paste0(
                       "Model calibration curve\n(", s$label, "; ",
                       n_per_bin, " per point)"))
})

# Held-out accuracy for the current settings, over the selected scope.
# Reports the whole metric set cv_performance() returns rather than only
# input$eval_stat: for a single configuration there is no cost to showing
# all of them, and they answer different questions (bias and R-squared
# disagree in useful ways). Deliberately NOT faded by gg_stale -- this
# section is computed from current settings, so unlike the heatmap it is
# never stale, and that contrast is the point.
output$eval_scatter_stats <- renderUI({
      req(is.null(eval_unsupported()))
      e <- eval_scatter_data()
      s <- eval_scatter_scope()
      req(e, s)
      if (sum(s$keep) < eval_min_n) return(NULL)   # scatter says why

      p <- tryCatch(analogs::cv_performance(e$cv[s$keep, , drop = FALSE]),
                    error = function(err) NULL)
      if (is.null(p) || !nrow(p)) return(NULL)

      fmt <- function(metric, value) {
            if (identical(metric, "n")) format(as.integer(value), big.mark = ",")
            else signif(value, 3)
      }
      items <- p %>%
            filter(metric %in% names(eval_metric_labels)) %>%
            arrange(match(metric, names(eval_metric_labels)))

      div(class = "mb-2", style = "font-size: 0.9em;",
          tags$span(style = "color: #666;", paste0(s$label, ": ")),
          lapply(seq_len(nrow(items)), function(i) {
                tags$span(style = "margin-right: 1.1em;",
                          tags$span(style = "color: #666;",
                                    eval_metric_labels[[items$metric[i]]]),
                          " ",
                          tags$strong(fmt(items$metric[i], items$value[i])))
          }))
})
