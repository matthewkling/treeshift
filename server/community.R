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
comm_run_cost <- reactive({
      req(multi_mode())
      n_ent <- ncol(multi_values())
      r <- clim()$hst
      n_cells <- terra::ncell(r)
      gb <- n_cells * n_ent * 8 * 2 / 1024^3
      secs <- (n_cells / 31624) * (0.10 + 0.012 * n_ent) * 2

      # Smallest grain that would fit, so a refusal can say what to do
      # rather than leaving the user to bisect the slider. Cells scale as
      # the inverse square of grain, so gb(g) = gb_now * (agg_now/g)^2.
      agg <- model$clim_agg
      fits <- ceiling(agg / sqrt(comm_max_gb / gb))
      list(entities = n_ent, cells = n_cells, gb = gb, secs = secs,
           # NA when even the coarsest grain the slider offers is not
           # enough -- then the bounding box is the thing to change.
           grain_that_fits = if (is.finite(fits) && fits <= 10) max(fits, 1) else NA_integer_)
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

      # sum_weights is retained even though nothing reads it yet: it is one
      # extra layer for the whole run (a weight property, not per-entity),
      # and it is the second half of the per-cell support signal. The
      # weighted count of contributing presence plots for entity s is
      # weighted_mean_s * sum_weights, which is what any reliability
      # filter or rare-entity floor will be built on.
      weight_stats <- "sum_weights"

      # Same positional rename as aim(): outcome layers are whatever isn't
      # a weight stat, returned in y-column order.
      name_outcomes <- function(r){
            is_outcome <- !names(r) %in% weight_stats
            stopifnot(sum(is_outcome) == ncol(y))
            names(r)[is_outcome] <- names(y)
            r
      }

      surf <- aim_surfaces(list(bsl = clim()$hst, fut = clim()$fut), pool(), y = y,
                           theta_clim = model$theta_clim, theta_geog = model$theta_geog,
                           stat = c("weighted_mean", weight_stats))
      lapply(surf, name_outcomes)
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
      m <- multi_aim()
      ent <- setdiff(names(m$bsl), "sum_weights")

      if (input$comm_measure %in% paired_measures) {
            # Biggest winner / loser: argmax and argmin of the per-entity
            # DIFFERENCE. These cannot be recovered from the reduced
            # surfaces -- which entity gained most is not a function of how
            # much the community gained -- so they are the categorical case
            # that genuinely needs the stack.
            va <- m$bsl[[ent]]
            vb <- m$fut[[ent]]
            dv <- vb - va
            idx <- if (identical(input$comm_measure, "biggest winner"))
                  terra::which.max(dv) else terra::which.min(dv)

            # Support over both eras together, so a cell that is empty now
            # but gains entities still qualifies. Requiring historical
            # support would mask exactly the arrivals this layer is for.
            support <- (terra::app(va, fun = "sum", na.rm = TRUE) +
                              terra::app(vb, fun = "sum", na.rm = TRUE)) / 2
      } else {
            # Dominant entity: per timeframe.
            r <- if (identical(input$time, "future")) m$fut else m$bsl
            v <- r[[ent]]
            idx <- terra::which.max(v)
            support <- terra::app(v, fun = "sum", na.rm = TRUE)
      }

      # Support floor FIRST. which.max always names a winner, including in
      # cells where essentially nothing is expected to be present -- so
      # without this the map shows a dominant tree species across
      # unforested country. Masking before the frequency count also keeps
      # those cells out of the top-N ranking, which otherwise promotes
      # whichever entity happens to win the empty ground.
      idx <- apply_support_floor(idx, support, input$support_floor)

      # Lump all but the most frequent N. Ranking by cell count rather
      # than mean value answers "what covers this map", which is the
      # question a discrete palette can address.
      n_top <- if (isTruthy(input$comm_top_n)) input$comm_top_n else 10
      f <- terra::freq(idx)
      f <- f[order(-f$count), , drop = FALSE]
      top <- utils::head(f$value, max(1, n_top))
      other <- setdiff(f$value, top)
      if (length(other)) idx <- terra::subst(idx, other, 0)
      names(idx) <- input$comm_measure

      codes  <- c(if (length(other)) 0, top)
      labels <- c(if (length(other)) "other", ent[top])
      # Looked up by entity name, so these match the composition panel
      # regardless of how the two views happen to rank entities.
      cols <- unname(comm_entity_colors()[labels[labels != "other"]])
      cols <- c(if (length(other)) "gray80", cols)

      list(raster = idx, codes = codes, labels = labels, colors = cols)
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
            return("No focal polygon yet. Choose one under 'Study area' in the left sidebar, or summarize over the whole bounding box instead.")
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

composition_layer <- reactive({
      req(multi_mode(), isTruthy(input$comm_measure))
      if (input$comm_measure %in% categorical_measures)
            return(comm_category_result()$raster)
      req(input$comm_measure %in% comm_measure_choices(model$track, model$fia_var, input$time))
      m <- multi_aim()

      reduce <- function(r){
            ent <- setdiff(names(r), "sum_weights")
            v <- r[[ent]]

            # Shared pieces. Everything below is layer-wise terra
            # arithmetic rather than an app() callback over cells, which
            # matters at 143 layers.
            s <- terra::app(v, fun = "sum", na.rm = TRUE)

            # Shannon entropy of the composition, normalized to
            # proportions q = v / sum(v):
            #   H = -sum(q log q) = log(s) - sum(v log v) / s
            # ifel() guards log(0); a zero entity contributes nothing to
            # the entropy, and 0 * -Inf would otherwise be NaN.
            shannon <- function(){
                  vlv <- terra::app(terra::ifel(v > 0, v * log(v), 0),
                                    fun = "sum", na.rm = TRUE)
                  log(s) - vlv / s
            }

            out <- switch(input$comm_measure,
                          # Expected richness. Summing occurrence
                          # probabilities is deliberate: thresholding at
                          # p > 0.5 would discard the difference between
                          # 0.51 and 0.99, and Sum(p) is the unbiased
                          # estimator if p is calibrated. Whether it IS
                          # calibrated across ~143 mostly-rare entities is
                          # the open question this layer exposes.
                          "richness" = s,

                          # Hill N1 -- the effective number of entities, in
                          # the same units as richness, so the two are
                          # directly comparable. A cell where one species
                          # dominates has low diversity however many rare
                          # species carry small probabilities, which is
                          # exactly the behavior richness lacks.
                          "diversity" = exp(shannon()),

                          # Pielou's J, against the FULL entity set rather
                          # than a count of entities present. With
                          # continuous probabilities almost every entity is
                          # nonzero somewhere, so a "number present" divisor
                          # would be degenerate. Note this makes evenness
                          # depend on the size of the entity set, which
                          # moves with theta_geog and the bbox.
                          "evenness" = shannon() / log(terra::nlyr(v)),

                          # Share of the community held by its largest
                          # member. Bounded 0-1 and needs no normalization.
                          "dominance" = terra::app(v, fun = "max", na.rm = TRUE) / s)
            req(!is.null(out))

            # app(na.rm = TRUE) returns 0 where every layer is NA, i.e.
            # outside analog range -- where a single-entity surface would
            # correctly be blank. Restore nodata there, or the map grows a
            # fringe of hard zeros that read as "no species" rather than
            # "no prediction". Every entity layer shares one NA mask, since
            # they come from one call with one set of weights and
            # multi_values() contains no NAs, so masking by any single
            # entity layer is enough.
            out <- terra::mask(out, v[[1]])

            # Support floor, applied per timeframe so a delta is blank
            # wherever EITHER timeframe is unsupported. `s` is the natural
            # denominator here: expected richness for an occurrence
            # variable, total basal area or stems for an abundance one.
            # Richness is exempt because it IS the denominator.
            if (!identical(input$comm_measure, "richness"))
                  out <- apply_support_floor(out, s, input$support_floor)

            names(out) <- input$comm_measure
            out
      }

      ## Paired reduction: both stacks at once, entity by entity. The
      ## per-entity differences cannot be recovered from the reduced
      ## surfaces -- Sum(max(f-h, 0)) is not a function of Sum(f) and
      ## Sum(h) -- which is what makes these genuinely require the stack
      ## rather than being obtainable from a cheap All-trees variable.
      ##
      ## `a` is the reference composition and `b` the comparison. For
      ## temporal pairing that is historical and future; the same formulas
      ## take a focal site's composition and every cell's for spatial
      ## pairing, which is why they are written in terms of a and b.
      reduce_paired <- function(a, b){
            ent <- setdiff(names(a), "sum_weights")
            va <- a[[ent]]
            vb <- b[[ent]]

            d      <- vb - va
            gains  <- terra::app(terra::ifel(d > 0,  d, 0), fun = "sum", na.rm = TRUE)
            losses <- terra::app(terra::ifel(d < 0, -d, 0), fun = "sum", na.rm = TRUE)
            sa     <- terra::app(va, fun = "sum", na.rm = TRUE)
            sb     <- terra::app(vb, fun = "sum", na.rm = TRUE)

            out <- switch(input$comm_measure,
                          # Bray-Curtis dissimilarity: sum|a-b| over
                          # sum(a+b), which is (gains + losses) over the
                          # combined community. 0 is identical
                          # composition, 1 is no overlap. Threshold-free,
                          # unlike Jaccard, which needs presence/absence
                          # sets and so needs a cutoff.
                          "turnover" = (gains + losses) / (sa + sb),

                          # Absolute expected arrivals and departures, in
                          # the variable's own units. Asymmetric on
                          # purpose: what must arrive is a different
                          # question from what must leave, and for
                          # sourcing or connectivity only one of them is
                          # the relevant demand.
                          "expected gains"  = gains,
                          "expected losses" = losses,

                          # Gains as a share of the future community: the
                          # fraction of what will be there that is not
                          # there now. This is the colonization demand used
                          # to weight connectivity, and it deliberately
                          # ignores losses.
                          "colonization demand" = gains / sb,

                          # Its mirror: losses as a share of the CURRENT
                          # community -- how much of what is here now is
                          # expected to depart. Normalizing by sa rather
                          # than sb is the whole difference, and it is what
                          # makes this the question a manager asks about
                          # the stand in front of them rather than about
                          # the one that replaces it.
                          #
                          # With the two fractions and turnover together,
                          # it is visible WHICH side drives turnover:
                          # Bray-Curtis is their average when the two
                          # communities are similar in size, so a high
                          # turnover with a low loss fraction is
                          # accumulation rather than replacement.
                          "loss fraction" = losses / sa)
            req(!is.null(out))

            out <- terra::mask(out, va[[1]])

            # Ratios only. Turnover is normalized by the combined
            # community, colonization demand by the future one, so each is
            # floored against its own denominator. The absolute sums stay
            # meaningful at low support and are left alone.
            if (identical(input$comm_measure, "turnover"))
                  out <- apply_support_floor(out, (sa + sb) / 2, input$support_floor)
            if (identical(input$comm_measure, "colonization demand"))
                  out <- apply_support_floor(out, sb, input$support_floor)
            if (identical(input$comm_measure, "loss fraction"))
                  out <- apply_support_floor(out, sa, input$support_floor)

            names(out) <- input$comm_measure
            out
      }

      if (input$comm_measure %in% paired_measures)
            return(reduce_paired(m$bsl, m$fut))

      hst <- reduce(m$bsl)
      fut <- reduce(m$fut)
      switch(input$time,
             "historical" = hst,
             "future"     = fut,
             "delta"      = fut - hst)
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
## is in the Study area panel -- and with one accordion panel open at a
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
                  "Even the coarsest grain size will not fit this bounding box. Shrink the box under 'Study area' in the left sidebar."
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
                              "Even the coarsest grain size will not fit this bounding box; shrink the box under 'Study area' in the left sidebar."
                        else sprintf("Set the grain size to %d km or coarser, or shrink the bounding box, under 'Study area' in the left sidebar.",
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
