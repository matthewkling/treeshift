## MULTI-ENTITY RUN AND REDUCTION -----------------
##
## One analog_impact() call per timeframe with every entity as an outcome
## column. Benchmarks: ~0.19 s for one entity and ~2 s for all 143 at the
## default domain and 10 km grain, because kernel weights are computed once
## and reused across columns. Cost is linear in cells, so this is cheap at
## coarse grain and minutes-scale at 1 km -- no button yet, deliberately,
## until a cost estimate says one is needed.
## Size and rough duration of a community run at the current settings.
## Both are read before the run rather than after, which is the point:
## the expensive configuration has to be refusable, not merely slow.
##
## Timing coefficients are measured, not guessed -- about 0.10 s fixed plus
## 0.012 s per entity at 31,624 cells, per timeframe, scaling linearly in
## cells. Marginal cost per entity climbs somewhat above ~25 entities, so
## this runs a little optimistic at large entity counts.
## The smallest grain that fits lets a refusal say what to do rather than
## leaving the user to bisect the slider (R/aim.R).
comm_run_cost <- reactive({
      req(multi_mode())
      ts_community_cost(terra::ncell(clim()$hst), ncol(multi_values()), model$clim_agg,
                        comm_max_gb)
})

multi_aim <- reactive({
      req(multi_mode())
      y <- multi_values()
      req(ncol(y) > 0)

      # req(), not validate(): this reactive is reached from the map's
      # leaflet observer, and a validation condition raised outside a
      # render context propagates as an app error rather than a message.
      # The refusal is explained by comm_cost_note in the sidebar, which
      # is a render context and can say it properly.
      cost <- comm_run_cost()
      req(cost$gb <= comm_max_gb)

      ts_multi_surfaces(clim(), pool(), y, model$theta_clim, model$theta_geog)
})

## Reduce the per-entity stack to one community surface. Reduction reads
## from the stack rather than accumulating during the run, which is what
## keeps every additional measure cheap -- re-reducing is a re-read, not a
## re-run.
## One entity -> color lookup for the whole app.
##
## Colors have to be keyed to entity IDENTITY, not to rank within a plot.
## The map ranks entities by cell count and the composition panel ranks
## them by value over the selected area, so assigning palette color 1 to
## each view's own first entity gives two entirely different mappings for
## the same species.
##
## The canonical order is domain-wide historical abundance, so the entities
## most likely to appear in any view get the first and most distinct colors.
## Past the palette length colors recycle, which can put two same-colored
## entities in one plot -- but only for entities rare enough to be far down
## this ordering, and a visible repeat is better than a mapping that
## silently disagrees between views.
comm_entity_colors <- reactive({
      m <- multi_aim()
      ent <- setdiff(names(m$bsl), "sum_weights")
      rank_val <- as.numeric(terra::global(m$bsl[[ent]], "mean", na.rm = TRUE)[, 1])
      ord <- ent[order(rank_val, decreasing = TRUE)]
      stats::setNames(cat_palette_colors(length(ord), last_cat_palette()), ord)
})

## Categorical measures: the argmax raster AND the code/label/color
## mapping, produced together.
##
## The mapping is carried here rather than as terra factor levels because
## assigning levels<- renumbers the categories to a sequential index --
## terra::cats() then returns 0..n-1 while the per-plot values are still
## column indices into the entity set, so almost nothing matches and the
## markers all collapse to "other".
comm_category_result <- reactive({
      req(categorical_layer())
      res <- ts_community_categories(multi_aim(), input$comm_measure, input$time,
                                     input$support_floor,
                                     if (isTruthy(input$comm_top_n)) input$comm_top_n else 10)
      # Looked up by entity name, so these match the composition panel
      # regardless of how the two views happen to rank entities.
      other <- res$labels == "other"
      res$colors <- c(if (any(other)) "gray80",
                      unname(comm_entity_colors()[res$labels[!other]]))
      res
})

## COMMUNITY COMPOSITION MODULE -------------------
##
## The map-level measures reduce across entities; this reduces across
## SPACE instead, so entity identity survives. Same stack, other axis --
## which is why it costs nothing beyond a global() over layers already
## computed.
comm_profile <- reactive({
      req(multi_mode())
      m <- multi_aim()
      ent <- setdiff(names(m$bsl), "sum_weights")
      scope <- if (isTruthy(input$comm_scope)) input$comm_scope else "domain"

      reduce_space <- function(r) {
            v <- r[[ent]]

            # One cell, not an average: the composition AT a place.
            if (scope == "cell") {
                  req(site$x, site$y)
                  return(as.numeric(terra::extract(v, cbind(site$x, site$y))))
            }

            mean_over <- function(g)
                  as.numeric(terra::global(terra::mask(terra::crop(v, g), g),
                                           "mean", na.rm = TRUE)[, 1])

            if (scope %in% c("polygon", "buffered")) {
                  # comm_profile_ready() explains both conditions to the user
                  req(isTruthy(polygon_active()), polygon_in_domain(scope))
                  return(mean_over(terra::vect(
                        if (scope == "polygon") user_polygon() else user_polygon_buffered())))
            }

            if (scope == "site") {
                  req(site$x, site$y, isTruthy(model$theta_geog))
                  return(mean_over(terra::buffer(
                        terra::vect(cbind(site$x, site$y), crs = "EPSG:4326"),
                        width = model$theta_geog * 3 * 1000)))
            }

            as.numeric(terra::global(v, "mean", na.rm = TRUE)[, 1])
      }

      hst <- reduce_space(m$bsl)
      fut <- reduce_space(m$fut)
      req(length(hst) == length(ent), length(fut) == length(ent))

      tibble::tibble(entity = ent, hst = hst, fut = fut) %>%
            # Entities absent in both timeframes are not part of this
            # community: they contribute zero to every reduction and only
            # crowd the label spread and the ranking.
            filter(hst > 0 | fut > 0) %>%
            mutate(change = fut - hst,
                   # Proportional change is undefined where nothing is
                   # there now, and wildly unstable just above that -- the
                   # same denominator problem the support floor addresses
                   # on the maps. Entities below the floor are dropped from
                   # the proportional ranking rather than topping it.
                   prop = ifelse(hst > 0, change / hst, NA_real_)) %>%
            filter(!is.na(hst), !is.na(fut))
})

## The entity set the panel names and colors, shared by the alluvial and
## the scatterplot so an entity has one color across both and the two show
## the same named set. Everything outside it is the "other" pool -- which
## is also exactly what the scatterplot's gray points are.
comm_bands <- reactive({
      d <- comm_profile()
      req(nrow(d) > 0)
      n <- if (isTruthy(input$comm_rank_n)) input$comm_rank_n else 15

      # Rank by the larger of the two values, so an entity arriving from
      # near zero is not excluded for having a small historical one.
      #
      # with_ties = FALSE is essential, not tidying: slice_max() otherwise
      # returns EVERY row tied at the cutoff value, and at a focal site
      # most entities sit at exactly 0. The n-th value is then 0, the tie
      # set is the whole pool, `rest` comes back empty, and nothing is ever
      # pooled into "other".
      d <- mutate(d, top = pmax(hst, fut, na.rm = TRUE))
      keep <- dplyr::slice_max(d, top, n = n, with_ties = FALSE)
      rest <- filter(d, !entity %in% keep$entity)

      # Largest first, so it sits at the bottom of the alluvial column and
      # the eye reads dominant entities against a fixed baseline.
      bands <- keep %>% arrange(desc(top)) %>% select(entity, hst, fut)

      # The tail is pooled rather than dropped: column height in the
      # alluvial is the community sum, and discarding ~128 entities would
      # make that aggregate silently wrong.
      if (nrow(rest))
            bands <- bind_rows(bands,
                               tibble::tibble(entity = "other",
                                              hst = sum(rest$hst, na.rm = TRUE),
                                              fut = sum(rest$fut, na.rm = TRUE)))

      # Same palette family as the dominant-entity map, so an entity keeps
      # a consistent identity across the map and both panel plots.
      cols <- comm_entity_colors()[bands$entity]
      cols[bands$entity == "other"] <- "gray80"
      names(cols) <- bands$entity

      list(bands = bands, cols = cols,
           named = bands$entity[bands$entity != "other"])
})

## Scope preconditions, stated rather than left as an empty panel.
comm_profile_ready <- reactive({
      # Most fundamental precondition first: someone in single-entity mode
      # with the scope on "focal polygon" should be told about the entity
      # setting, not sent off to choose a polygon that would not help.
      if (!multi_mode())
            return("This panel doesn't work in single-entity mode. Under 'Entity' in the left sidebar, switch from a single species or type to 'community composition'.")
      scope <- if (isTruthy(input$comm_scope)) input$comm_scope else "domain"
      if (scope %in% c("polygon", "buffered") && !isTruthy(polygon_active()))
            return("No focal polygon yet. Choose one under 'Focal polygon and site' in the left sidebar, or summarize over the whole bounding box instead.")
      if (scope %in% c("polygon", "buffered") && !polygon_in_domain(scope))
            return(polygon_outside_msg)
      if (scope %in% c("site", "cell") && (is.null(site$x) || is.null(site$y)))
            return("No focal site yet. Click anywhere on the map to set one.")
      NULL
})

output$comm_scatter <- renderPlot({
      shiny::validate(need(is.null(comm_profile_ready()), comm_profile_ready()))
      d <- comm_profile()
      req(nrow(d) > 0)
      cb <- comm_bands()

      # Named entities are the same set the alluvial draws, in the same
      # colors, so the two plots can be read against each other. The gray
      # points are literally the alluvial's "other" band, unpooled.
      d <- mutate(d, band = ifelse(entity %in% cb$named, entity, "other"))
      lab <- filter(d, band != "other") %>% arrange(desc(abs(change)))

      ggplot(d, aes(hst, change)) +
            geom_hline(yintercept = 0, color = "gray70", linewidth = .3) +
            geom_point(data = filter(d, band == "other"),
                       color = "gray80", size = 1.5) +
            geom_point(data = lab, aes(color = entity), size = 2.4) +
            # check_overlap rather than a repel dependency: rows are
            # pre-sorted by |change| so the most extreme entities win the
            # collisions.
            geom_text(data = lab, aes(label = entity, color = entity),
                      size = 3.1, hjust = -0.12, vjust = 0.5,
                      check_overlap = TRUE) +
            scale_color_manual(values = cb$cols, guide = "none") +
            scale_x_sqrt() +
            ggplot2::expand_limits(x = max(d$hst, na.rm = TRUE) * 1.35) +
            theme_bw() +
            labs(x = paste("historical", tolower(model$fia_var)),
                 y = "change",
                 subtitle = names(area_names)[match(
                       if (isTruthy(input$comm_scope)) input$comm_scope else "domain",
                       area_names)])
})

output$comm_ranking <- renderPlot({
      shiny::validate(need(is.null(comm_profile_ready()), comm_profile_ready()))
      d <- comm_profile()
      req(nrow(d) > 0)

      prop <- identical(input$comm_rank_by, "proportional")
      d$value <- if (prop) d$prop else d$change
      d <- filter(d, !is.na(value), value != 0)
      shiny::validate(need(nrow(d) > 0, "Nothing changes measurably over this area."))

      n <- if (isTruthy(input$comm_rank_n)) input$comm_rank_n else 15
      # Both tails, not the top of one: the losers are as much the answer
      # as the winners, and showing only gains hides half the turnover.
      # with_ties = FALSE for the same reason as comm_bands(): a tie at
      # the cutoff would return the whole pool.
      d <- bind_rows(dplyr::slice_max(d, value, n = ceiling(n / 2), with_ties = FALSE),
                     dplyr::slice_min(d, value, n = floor(n / 2), with_ties = FALSE)) %>%
            distinct(entity, .keep_all = TRUE) %>%
            mutate(entity = factor(entity, levels = entity[order(value)]))

      ggplot(d, aes(value, entity, fill = value > 0)) +
            geom_col() +
            geom_vline(xintercept = 0, color = "gray30", linewidth = .3) +
            scale_fill_manual(values = c(`TRUE` = "dodgerblue4", `FALSE` = "darkred"),
                              guide = "none") +
            theme_bw() +
            labs(y = NULL,
                 x = if (prop) "proportional change" else "change")
})

## The quantitative community measures (R/layers.R)
composition_layer <- reactive({
      req(multi_mode(), isTruthy(input$comm_measure))
      if (input$comm_measure %in% categorical_measures)
            return(comm_category_result()$raster)
      req(input$comm_measure %in% comm_measure_choices(model$track, model$fia_var, input$time))
      ts_community_layer(multi_aim(), input$comm_measure, input$time, input$support_floor)
})
## The support floor's units are the denominator's units, which change with
## the variable and with entity scope.
##
## The value is reset only when the UNITS change, not on every variable
## change: moving between two occurrence variables keeps whatever the user
## set, while moving from expected species to ft2/acre resets, because
## carrying a threshold of 1 across that boundary would mean something
## entirely different. Passing value = NULL leaves the field alone.
## Starts at the units of the settings the page opened with, so the first run
## keeps the starting value -- a linked support floor in particular.
support_floor_units <- reactiveVal(
      support_floor_spec(url_multi(link$state), link$state$fia_var, link$state$track)$label)
observeEvent(list(multi_mode(), model$fia_var, model$track), {
      req(isTruthy(model$fia_var))
      sp <- support_floor_spec(multi_mode(), model$fia_var, model$track)
      units_changed <- !identical(support_floor_units(), sp$label)
      support_floor_units(sp$label)
      updateNumericInput(session, "support_floor",
                         label = sp$label, step = sp$step,
                         value = if (units_changed) sp$default else NULL)
})

## Repopulate community measures for the current track and variable.
## Invalid combinations are removed from the menu rather than flagged
## after selection. observeEvent so the handler's read of
## input$comm_measure is isolated and cannot re-trigger its own write.
observeEvent(list(model$track, model$fia_var, input$time), {
      req(multi_mode())
      ch <- comm_measure_choices(model$track, model$fia_var, input$time)
      req(length(ch) > 0)
      cur <- if (isTruthy(input$comm_measure)) input$comm_measure else ""
      updateSelectInput(session, "comm_measure", choices = ch,
                        selected = if (cur %in% ch) cur else ch[1])
})

## A modal when the run first goes over the limit.
##
## The inline notice is not enough on its own: it renders in the FIA
## variable panel, while the grain slider that pushes a run over the limit
## is in the Raster grid panel -- and with one accordion panel open at a
## time, the notice is suspended at exactly the moment it matters. So the
## user coarsens nothing, sees the map quietly go blank, and has no idea
## why.
##
## Fires once per crossing rather than on every slider tick: the flag is
## read with isolate() so the observer does not depend on what it writes,
## and clears when the settings come back under the limit.
comm_limit_warned <- reactiveVal(FALSE)

observe({
      req(multi_mode())
      cost <- comm_run_cost()

      if (cost$gb <= comm_max_gb) {
            if (isolate(comm_limit_warned())) comm_limit_warned(FALSE)
            return()
      }
      if (isolate(comm_limit_warned())) return()
      comm_limit_warned(TRUE)

      fix <- cost$grain_that_fits
      showModal(modalDialog(
            title = "Too large to map",
            tags$p(sprintf(
                  "This query would model %s species/forest types over %s cells, which would need about %.1f GB of memory, exceeding this app's %.1f GB limit. The operation will not be run at these settings.",
                  format(cost$entities, big.mark = ","),
                  format(cost$cells, big.mark = ","),
                  cost$gb, comm_max_gb)),
            tags$p(if (is.na(fix))
                  "Even the coarsest grain size will not fit this bounding box. Shrink the box under 'Raster grid' in the left sidebar."
                  else sprintf("You could increase grain size to %d km or coarser, shrink the bounding box, or choose a single-species variable instead of `community composition`.", fix)),
            easyClose = TRUE,
            footer = tagList(
                  # Acting on the advice directly, rather than sending the
                  # user to find a slider in a panel they are not in.
                  if (!is.na(fix))
                        actionButton("comm_apply_grain",
                                     sprintf("Set grain to %d km", fix)),
                  modalButton("Close")
            )
      ))
})

observeEvent(input$comm_apply_grain, {
      cost <- comm_run_cost()
      req(!is.na(cost$grain_that_fits))
      updateSliderInput(session, "clim_agg", value = cost$grain_that_fits)
      removeModal()
})

## What the current settings will cost, shown before it is incurred rather
## than only when it is refused. Silent for runs that are quick enough not
## to warrant a warning: a notice on every community run would be noise,
## and noise is what gets ignored when it matters.
output$comm_cost_note <- renderUI({
      req(multi_mode())
      cost <- comm_run_cost()
      if (cost$gb <= comm_max_gb && cost$secs < 10) return(NULL)

      if (cost$gb > comm_max_gb) {
            div(class = "ts-notice",
                sprintf("This query would model %s species/forest types over %s cells, which would need about %.1f GB of memory, exceeding this app's %.1f GB limit. The operation will not be run at these settings.",
                        format(cost$entities, big.mark = ","),
                        format(cost$cells, big.mark = ","),
                        cost$gb, comm_max_gb,
                        if (is.na(cost$grain_that_fits))
                              "Even the coarsest grain size will not fit this bounding box; shrink the box under 'Raster grid' in the left sidebar."
                        else sprintf("Set the grain size to %d km or coarser, or shrink the bounding box, under 'Raster grid' in the left sidebar.",
                                     cost$grain_that_fits)))
      } else {
            div(class = "ts-notice",
                sprintf("Mapping %s entities over %s cells takes roughly %s and about %.1f GB. Coarsening the grain size is the quickest way to speed this up.",
                        format(cost$entities, big.mark = ","),
                        format(cost$cells, big.mark = ","),
                        if (cost$secs < 90) sprintf("%.0f seconds", cost$secs)
                        else sprintf("%.0f minutes", cost$secs / 60),
                        cost$gb))
      }
})

## Only one case remains that the menu cannot express: no measure at all
## applies, so there is nothing to select and the reason needs saying.
output$comm_note <- renderUI({
      req(multi_mode(), isTruthy(model$fia_var))
      req(length(comm_measure_choices(model$track, model$fia_var, input$time)) == 0)

      msg <- if (!model$fia_var %in% multi_vars) {
            sprintf("'%s' has no per-entity form yet \u2014 the regeneration fractions are pooled ratios, which need three coupled columns per entity. Choose another variable, or switch back to a single entity.",
                    model$fia_var)
      } else {
            sprintf("Every FIA plot is assigned a single %s, so counting how many occur together at one plot always gives one. Set Timeframe to 'change' for the turnover measures, which compare how the mix of classes shifts over time and do work here.",
                    tolower(model$track))
      }
      div(class = "ts-notice", msg)
})
