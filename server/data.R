## DATA SETUP -----------------------

# Padded bbox: includes plots just outside the model bounding box so they
# can still serve as analogs for cells near the bbox edge.
# Pad tracks theta_geog * 3 (the hard cutoff used by analog_impact),
# converted from km to degrees with a conservative ratio. 1° lat ≈ 111 km;
# 1° lon at 50° N ≈ 71 km. Using /100 overestimates the pad by 10–30% at
# CONUS latitudes, which is harmless.
pool_bbox <- reactive({
      pad_deg <- input$theta_geog * 3 / 100
      list(x = c(input$bbox_x[1] - pad_deg, input$bbox_x[2] + pad_deg),
           y = c(input$bbox_y[1] - pad_deg, input$bbox_y[2] + pad_deg))
})

# all plots within padded bbox (the analog pool)
fia_plots <- reactive({
      pb <- pool_bbox()
      fia_all %>%
            select(x, y) %>%
            distinct() %>%
            filter(between(x, pb$x[1], pb$x[2]),
                   between(y, pb$y[1], pb$y[2])) %>%
            mutate(index = 1:nrow(.))
})

fia_plots_species_all <- reactive({
      fia_all %>%
            # Recode NA forest_type to the Non-forest sentinel (0) so
            # downstream filters and joins treat non-forest plots
            # uniformly with real forest types.
            mutate(forest_type = ifelse(is.na(forest_type), 0L, forest_type)) %>%
            left_join(select(fia_plots(), x, y, index),
                      by = join_by(x, y))
})

plot_occ <- reactive({
      pb <- pool_bbox()
      # Plot-level type->group lookup (one row per forest_type code)
      ft_to_fg <- distinct(forest_type_codes, forest_type, forest_group)

      fia_plots_species_all() %>%
            filter(between(x, pb$x[1], pb$x[2]),
                   between(y, pb$y[1], pb$y[2])) %>%
            mutate(focal = species == input$species & !is.na(species)) %>%
            group_by(x, y) %>%
            summarize(
                  samp_method_cd = mean(samp_method_cd),
                  pres = as.integer(any(focal)), # focal species presence/absence (any life stage)
                  pres_seed = as.integer(any(focal & n_seedlings > 0)), # focal species seedling presence
                  pres_sap = as.integer(any(focal & n_saplings > 0)), # focal species sapling presence
                  pres_tree = as.integer(any(focal & n_trees > 0)), # focal species adult (tree) presence
                  batot = sum(basal_area[focal]), # total focal species basal area
                  baprop = weighted.mean(focal, basal_area), #  proportion ba
                  baplot = sum(basal_area), # total plot ba, all species
                  nttot = sum(n_trees[focal]), # focal species tree count
                  ntplot = sum(n_trees), # total plot tree count, all species
                  satot = sum(n_saplings[focal]), # focal species sapling count
                  saplot = sum(n_saplings), # total plot sapling count, all species
                  nstot = sum(n_seedlings[focal]), # focal species seedling count
                  nsplot = sum(n_seedlings), # total plot seedling count, all species
                  # Plot species richness, by stage. A count of distinct
                  # species is an aggregation ACROSS species done at the
                  # plot, exactly like baplot and ntplot -- which is what
                  # makes it an All-trees variable rather than something
                  # only the multi-entity path can produce. Because AIM is
                  # a weighted mean (linear) with entity-independent
                  # weights, interpolating this equals summing the
                  # interpolated per-entity occurrence surfaces. The two
                  # routes should agree to floating point, which makes each
                  # a check on the other.
                  nspp      = dplyr::n_distinct(species[!is.na(species)]),
                  nspp_seed = dplyr::n_distinct(species[!is.na(species) & n_seedlings > 0]),
                  nspp_sap  = dplyr::n_distinct(species[!is.na(species) & n_saplings > 0]),
                  nspp_tree = dplyr::n_distinct(species[!is.na(species) & n_trees > 0]),
                  forest_type = dplyr::first(forest_type), # plot-level FORTYPCD
                  .groups = "drop") %>%
            left_join(ft_to_fg, by = "forest_type") %>%
            mutate(baprop = ifelse(baplot == 0, 0, baprop),
                   # Per-plot regeneration fractions (focal species), on an
                   # occurrence basis to match the interpolated surface (see
                   # ratio_pieces): 0 = adults only, 0.5 = stage and adults
                   # both present, 1 = stage only. NA where the species is
                   # absent as both that stage and adults.
                   seed_frac = ifelse((nstot + nttot) > 0,
                                      (nstot > 0) / ((nstot > 0) + (nttot > 0)), NA_real_),
                   sap_frac  = ifelse((satot + nttot) > 0,
                                      (satot > 0) / ((satot > 0) + (nttot > 0)), NA_real_))
})

# Single active response vector, aligned with plot_occ() rows.
# Computed from the (track, fia_var, focal_id) tuple.
response_values <- reactive({
      po <- plot_occ()
      req(isTruthy(input$track), isTruthy(input$fia_var))

      v <- if (input$track == "Species") {
            switch(input$fia_var,
                   "presence probability"  = po$pres,
                   "seedling presence"      = po$pres_seed,
                   "sapling presence"       = po$pres_sap,
                   "tree presence"          = po$pres_tree,
                   "proportion basal area" = po$baprop,
                   "total basal area"      = po$batot,
                   "number of trees"       = po$nttot,
                   "number of saplings"    = po$satot,
                   "number of seedlings"   = po$nstot,
                   "seedling fraction"      = po$seed_frac,
                   "sapling fraction"       = po$sap_frac)
      } else if (input$track == "Forest type") {
            req(isTruthy(input$forest_type))
            as.integer(!is.na(po$forest_type) &
                             po$forest_type == as.integer(input$forest_type))
      } else if (input$track == "Forest group") {
            req(isTruthy(input$forest_group))
            as.integer(!is.na(po$forest_group) &
                             po$forest_group == as.integer(input$forest_group))
      } else if (input$track == "All trees") {
            switch(input$fia_var,
                   "presence probability" = as.integer((po$ntplot + po$saplot) > 0 | po$nsplot > 0),
                   "seedling presence"    = as.integer(po$nsplot > 0),
                   "sapling presence"     = as.integer(po$saplot > 0),
                   "tree presence"        = as.integer(po$ntplot > 0),
                   "total basal area"     = po$baplot,
                   "number of trees"      = po$ntplot,
                   "number of saplings"   = po$saplot,
                   "number of seedlings"  = po$nsplot,
                   "species richness"          = po$nspp,
                   "seedling species richness" = po$nspp_seed,
                   "sapling species richness"  = po$nspp_sap,
                   "adult species richness"    = po$nspp_tree,
                   # Occurrence-based, matching ratio_pieces() and the
                   # Species-track fractions in plot_occ(): 0 = adults only,
                   # 0.5 = stage and adults both present, 1 = stage only.
                   "seedling fraction"    = ifelse((po$nsplot + po$ntplot) > 0,
                                                   (po$nsplot > 0) / ((po$nsplot > 0) + (po$ntplot > 0)), NA_real_),
                   "sapling fraction"     = ifelse((po$saplot + po$ntplot) > 0,
                                                   (po$saplot > 0) / ((po$saplot > 0) + (po$ntplot > 0)), NA_real_))
      }

      # Guard against transient NULL during track/fia_var transitions
      req(!is.null(v))
      v
})

# Whether the focal entity is present, per track -- the same quantity the
# "presence probability" variable measures, and independent of which
# variable is being displayed. Used by the "occupied FIA plots" filter:
# keying that off `display > 0` means "has regeneration" for the fraction
# variables, which hides plots that are occupied but simply aren't
# regenerating. (It coincides with display > 0 for the other variables.)
## Per-plot value to paint the FIA markers with: the observed counterpart
## to whatever the raster is showing. In Single mode that is the focal
## entity's response; in All mode it is the same reduction applied across
## entities at each plot, so observed richness sits next to modelled
## richness on the same color scale. That comparison is the cheapest
## available check on whether Sum(p_s) is calibrated.
display_values <- reactive({
      if (!multi_mode()) return(response_values())

      # Guard before touching multi_values(). The regeneration fractions
      # have no per-entity form, so multi_values() req()s out for them --
      # and because it would otherwise be passed to comm_reduce_rows() as
      # an unforced promise, the reactiveStop would surface from inside
      # as.matrix()'s S4 dispatch as an error rather than a silent cancel.
      # Forcing it here, in this reactive's own frame, keeps a cancellation
      # a cancellation.
      if (!isTruthy(input$fia_var) || !input$fia_var %in% multi_vars) return(NULL)
      y <- multi_values()
      req(!is.null(y), ncol(y) > 0)

      v <- comm_reduce_rows(y, input$comm_measure)
      # NULL means the measure has no per-plot counterpart (the paired
      # categoricals: a plot is surveyed at one point in time, so there is
      # no observed "biggest winner"). Return early -- subsetting NULL in
      # the lumping step below yields a ZERO-LENGTH vector rather than
      # NULL, which slips past a null check and only fails later at
      # `pd$display <- dv`.
      if (is.null(v)) return(NULL)

      # Lump the observed categories to exactly the set the raster drew, so
      # a plot and the cell under it get the same color for the same
      # entity, and anything outside the top N reads as "other" on both.
      if (input$comm_measure %in% categorical_measures) {
            cc <- comm_categories()
            v[!is.na(v) & !v %in% cc$codes] <- 0L
      }
      v
})

occupied_values <- reactive({
      po <- plot_occ()
      req(isTruthy(input$track))

      # In All mode "occupied" means any entity is present -- a forested
      # plot -- since there is no focal entity to be occupied by.
      #
      # Assigned before use for the same reason as display_values(): forcing
      # multi_values() inside an operator means a failed req() escapes from
      # dispatch as an error rather than cancelling.
      if (multi_mode()) {
            if (!isTruthy(input$fia_var) || !input$fia_var %in% multi_vars)
                  return(rep(FALSE, nrow(po)))
            y <- multi_values()
            return(rowSums(y > 0, na.rm = TRUE) > 0)
      }

      v <- if (input$track == "Species") {
            po$pres
      } else if (input$track == "Forest type") {
            req(isTruthy(input$forest_type))
            as.integer(!is.na(po$forest_type) &
                             po$forest_type == as.integer(input$forest_type))
      } else if (input$track == "Forest group") {
            req(isTruthy(input$forest_group))
            as.integer(!is.na(po$forest_group) &
                             po$forest_group == as.integer(input$forest_group))
      } else if (input$track == "All trees") {
            as.integer((po$ntplot + po$saplot) > 0 | po$nsplot > 0)
      }

      req(!is.null(v))
      as.logical(v)
})

# Plot-level data behind the FIA marker layer, with the active response
# attached as `display` and focal-entity occupancy as `occupied`. Shared by
# the map observer and the CSV export so the downloaded table is exactly
# the set of points drawn on the map.
#
# When the vector selector is on a non-point layer ("none", current
# animation/field) there is nothing on screen to mirror, so this falls back
# to the full analog pool -- the export stays available and meaningful
# regardless of what the map happens to be showing.
vector_points <- reactive({
      req(plot_occ())
      # display_values() is NULL if a community measure ever lacks a
      # row-wise counterpart. Catching it here degrades to an empty layer;
      # letting it through assigns NULL to pd$display, which silently drops
      # the column and crashes inside addCircles() instead.
      dv <- display_values()
      # Length, not just non-NULL: a zero-length vector would fail later
      # inside the data-frame assignment with a recycling error rather
      # than cancelling here.
      req(!is.null(dv), length(dv) == nrow(plot_occ()))

      if(input$vector == "focal site analogs"){
            pd <- switch(input$analog_direction,
                         "Reverse analogs" = site_analogs()$rev_analogs,
                         "Forward analogs" = site_analogs()$fwd_analogs,
                         "Contemporary analogs" = site_analogs()$bsl_analogs) %>%
                  mutate(x = analog_x, y = analog_y)
            # Index the per-plot value by analog_index rather than
            # reading a `response` y-column: pool()$values carries
            # occ/num/denom (not `response`) for the pooled-ratio
            # variables, so display_values() is the uniform source.
            pd$display <- dv[pd$analog_index]
            pd$occupied <- occupied_values()[pd$analog_index]
      } else {
            # Show the full analog pool, including plots in the padded
            # ring just outside the modeling bbox — these are still
            # contributing as analogs to edge cells in the prediction.
            pd <- plot_occ()
            pd$display <- dv
            pd$occupied <- occupied_values()
      }

      # Occupancy is the focal entity's presence, not `display > 0`:
      # for the regeneration fractions the latter would mean "has
      # regeneration" and hide adults-only plots, which are occupied.
      if(input$vector == "occupied FIA plots") pd <- filter(pd, occupied)

      pd
})

# Header-only version of the export schema, used when export_points()
# can't be built (see the download handler). Keeps the file a valid CSV
# with the documented columns rather than a Shiny error page.
export_points_empty <- function(){
      tibble::tibble(x = numeric(), y = numeric(),
                     entity = character(), variable = character(),
                     value = numeric(), occupied = logical(),
                     sampled = logical())
}

# Curated export schema for the FIA point download.
#
# Deliberately narrow, and deliberately FIXED: the header is identical no
# matter which variable is displayed, because the variable's name lives in
# a `variable` column rather than in the column name of its values. That
# means several exports can be bind_rows()'d, and a script can read `value`
# without knowing what was on screen when the file was made.
#
# Everything else in plot_occ() (batot, nsplot, samp_method_cd, ...) is
# internal shorthand. Emitting it would turn those names into a public
# interface that has to be documented here and kept stable across refactors,
# which is a bigger commitment than the columns are worth.
export_points <- reactive({
      pd <- vector_points()

      out <- tibble::tibble(
            x = pd$x,
            y = pd$y,
            entity = focal_label(),
            # The forest-type/group tracks ignore fia_var: their response is
            # membership in the type named by `entity`.
            variable = if(input$track %in% c("Forest type", "Forest group")) "presence"
            else input$fia_var,
            value = pd$display,
            # Focal entity presence, independent of `variable`. Not the same
            # question as value > 0: for the regeneration fractions that
            # would mean "has regeneration", which is false at occupied
            # adults-only plots.
            occupied = pd$occupied,
            # FALSE at plots that exist in the frame but were never surveyed
            # (samp_method_cd 2). Without this their zeros read downstream as
            # true absences.
            sampled = pd$samp_method_cd == 1
      )

      # The analog layer additionally carries the kernel geometry that
      # defines the analog set. Weight matches the plot-level decomposition
      # code: the product of the two gaussian kernel densities, unnormalized.
      if(all(c("env_dist", "geog_dist") %in% names(pd))){
            out <- out %>%
                  mutate(climate_distance = pd$env_dist,
                         geographic_distance = pd$geog_dist,
                         weight = dnorm(pd$env_dist, sd = input$theta_clim) *
                               dnorm(pd$geog_dist, sd = input$theta_geog))
      }

      out
})

# Pooled-ratio variables (proportion basal area, seedling/sapling fraction)
# are interpolated as a ratio-of-means rather than a mean-of-ratios: we
# interpolate the numerator and denominator separately and divide the
# surfaces (then gate by occurrence in aim()). This returns the per-plot
# numerator/denominator for those variables, or NULL for everything else
# (which uses the single-response path). Occurrence is denom > 0.
## MULTI-ENTITY RESPONSE MATRIX -------------------
##
## plot_occ() collapses to one row per plot around a single focal entity
## (`focal = species == input$species`), so there is nowhere for a
## multi-column `y` to come from. This builds the wide counterpart: one
## row per plot, one column per entity, holding the same quantity
## response_values() would return if that entity were the focal one.
##
## Row alignment with plot_occ() is the whole ballgame. pool()$index_hst is
## built from plot_occ()[, c("x","y")], so a mismatch here silently attaches
## every entity's response to the wrong plot and produces maps that look
## entirely plausible. Alignment is guaranteed by construction (the join
## runs FROM plot_occ()'s coordinates) and then asserted on the coordinates
## themselves rather than on row counts.
multi_values <- reactive({
      req(isTruthy(input$track))
      po <- plot_occ()
      xy <- po[, c("x", "y")]

      if (input$track %in% c("Forest type", "Forest group")) {
            # Membership is plot-level -- one type per plot -- so the wide
            # form is an indicator matrix and fia_var does not apply. Same
            # reason response_values() ignores fia_var for these tracks.
            key <- if (input$track == "Forest type") po$forest_type else po$forest_group
            lvls <- sort(unique(key[!is.na(key)]))
            req(length(lvls) > 0)
            w <- as.data.frame(
                  vapply(lvls, function(k) as.integer(!is.na(key) & key == k),
                         integer(nrow(po))),
                  check.names = FALSE)

            # Name the columns for the classes rather than their FORTYPCD
            # integers. The column names ARE the entity labels everywhere
            # downstream -- legend, category table, and any future ranking
            # -- so resolving them here fixes all of those at once rather
            # than translating codes at each display site.
            nm <- if (input$track == "Forest type") ft_name(as.integer(lvls))
            else fg_name(as.integer(lvls))
            nm[is.na(nm)] <- paste(input$track, lvls[is.na(nm)])
            names(w) <- make.unique(as.character(nm))

      } else {
            req(isTruthy(input$fia_var), input$fia_var %in% multi_vars)
            pb <- pool_bbox()

            # Same summaries plot_occ() computes, but grouped by entity as
            # well as plot. All of them at once: it is a single pass over
            # the same rows either way, and computing them together keeps
            # the definitions visibly identical to the single-entity path.
            long <- fia_plots_species_all() %>%
                  filter(between(x, pb$x[1], pb$x[2]),
                         between(y, pb$y[1], pb$y[2]),
                         !is.na(species)) %>%
                  group_by(x, y, species) %>%
                  summarize(
                        # Presence is "the entity was recorded at this
                        # plot", matching plot_occ()'s any(focal) -- it is
                        # not conditioned on any stage count.
                        pres      = 1L,
                        pres_seed = as.integer(any(n_seedlings > 0)),
                        pres_sap  = as.integer(any(n_saplings > 0)),
                        pres_tree = as.integer(any(n_trees > 0)),
                        batot     = sum(basal_area),
                        nttot     = sum(n_trees),
                        satot     = sum(n_saplings),
                        nstot     = sum(n_seedlings),
                        .groups = "drop")

            val_col <- switch(input$fia_var,
                              "presence probability"  = "pres",
                              "seedling presence"     = "pres_seed",
                              "sapling presence"      = "pres_sap",
                              "tree presence"         = "pres_tree",
                              "total basal area"      = "batot",
                              "proportion basal area" = "batot",
                              "number of trees"       = "nttot",
                              "number of saplings"    = "satot",
                              "number of seedlings"   = "nstot")

            wide <- long %>%
                  select(x, y, species, value = dplyr::all_of(val_col)) %>%
                  tidyr::pivot_wider(names_from = species, values_from = value,
                                     values_fill = 0)

            # Join FROM plot_occ()'s coordinates, so row order is inherited
            # rather than assumed. Plots where no entity was recorded are
            # dropped by the !is.na(species) filter above and must come
            # back as all-zero rows: absence is data, and dropping those
            # plots would bias every entity's surface upward.
            j <- xy %>% left_join(wide, by = dplyr::join_by(x, y))
            stopifnot(nrow(j) == nrow(xy),
                      identical(j$x, xy$x), identical(j$y, xy$y))

            w <- as.data.frame(j[, setdiff(names(j), c("x", "y")), drop = FALSE],
                               check.names = FALSE)
            w[is.na(w)] <- 0

            if (identical(input$fia_var, "proportion basal area")) {
                  # Denominator is the plot total across ALL entities, so
                  # it is shared. Plots with no basal area get 0 rather
                  # than NaN, matching plot_occ()'s baprop guard.
                  denom <- po$baplot
                  w[] <- lapply(w, function(v) ifelse(denom == 0, 0, v / denom))
            }
      }

      # Occupied-plot count per entity: the pool-level support signal, and
      # the basis for deciding a rare-entity floor. Attached rather than
      # filtered on, so that decision can be made against real numbers
      # instead of a guessed threshold.
      attr(w, "n_plots") <- vapply(w, function(v) sum(v > 0, na.rm = TRUE), integer(1))
      w
})

ratio_pieces <- reactive({
      req(isTruthy(input$track), isTruthy(input$fia_var))
      po <- plot_occ()
      # Regeneration fractions are a fraction of OCCURRENCE weight, not of
      # stem density. Seedlings and saplings are tallied on the microplot
      # (each stem = 74.965/acre) while adults come off the subplot (each =
      # 6.018/acre), so a density ratio applies a ~12.5x amplification to
      # every regen stem before dividing: one seedling + one adult reads
      # 0.93 rather than 0.5, squashing the balanced core up against the
      # pure-colonization end. Using presence indicators puts both stages on
      # the same footing and gives a readable axis: 0 = adults only
      # (relict/trailing), ~0.5 = both present (core), 1 = regen only
      # (colonizing/leading).
      #
      # Proportion basal area is deliberately left on a density basis: it is
      # focal BA over stand BA, with each stem's own TPA applied to both
      # sides, so it never crosses the microplot/subplot boundary and isn't
      # distorted the way the regen fractions are.
      if (input$track == "Species") {
            switch(input$fia_var,
                   "proportion basal area" = list(num = po$batot, denom = po$baplot),
                   "seedling fraction"     = list(num = as.integer(po$nstot > 0),
                                                  denom = as.integer(po$nstot > 0) +
                                                        as.integer(po$nttot > 0)),
                   "sapling fraction"      = list(num = as.integer(po$satot > 0),
                                                  denom = as.integer(po$satot > 0) +
                                                        as.integer(po$nttot > 0)),
                   NULL)
      } else if (input$track == "All trees") {
            switch(input$fia_var,
                   "seedling fraction" = list(num = as.integer(po$nsplot > 0),
                                              denom = as.integer(po$nsplot > 0) +
                                                    as.integer(po$ntplot > 0)),
                   "sapling fraction"  = list(num = as.integer(po$saplot > 0),
                                              denom = as.integer(po$saplot > 0) +
                                                    as.integer(po$ntplot > 0)),
                   NULL)
      } else NULL
})

clim_scales <- readRDS("data/clim_scales.rds")

# The focal polygon is a unit from the boundary library (helpers/boundaries.R),
# chosen by its stable id. Users cannot upload their own: reading user-supplied
# files on the server was ruled an IT security risk.
#
# polygon_active() is the selected unit id, or NULL when none is selected; it is
# what the rest of the app tests for "is there a polygon".
polygon_active <- reactive({
      if (isTruthy(input$polygon_unit)) input$polygon_unit else NULL
})

observeEvent(input$polygon_clear, {
      updateSelectizeInput(session, "polygon_unit", selected = "")
})

user_polygon <- reactive({
      req(polygon_active())
      boundary_units[boundary_units$id == polygon_active(), ]
})

# Summaries over the polygon (or its buffer) only see the part inside the model
# bounding box: the modeled grid and the analog pool are both cut to it. With
# the library one click away, a polygon entirely outside the box is an ordinary
# case, so every module says so rather than failing or coming up empty.
# Geometry, not extent: a polygon's bounding rectangle can touch the box when
# the polygon itself does not.
polygon_in_domain <- function(scope = "polygon") {
      g <- if (identical(scope, "buffered")) user_polygon_buffered() else user_polygon()
      terra::is.related(terra::vect(g), terra::ext(c(input$bbox_x, input$bbox_y)), "intersects")
}

polygon_outside_msg <- "The focal polygon is outside the model bounding box. Use 'Snap box to polygon' under 'Study area' to move the box to it."


## focal site
site <- reactiveValues(x = NULL, y = NULL)
