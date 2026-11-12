## DATA SETUP -----------------------

# The analog pool's bounding box: the model box padded by the kernel's reach
pool_bbox <- reactive(ts_pool_bbox(model$bbox_x, model$bbox_y, model$theta_geog))

# One row per FIA plot in the pool (R/plots.R). Its rows are the plot index
# that the pool and every response are aligned to.
plot_occ <- reactive(ts_plots(core_data, pool_bbox(), model$species))

# Single active response vector, aligned with plot_occ() rows
response_values <- reactive({
      req(isTruthy(model$track), isTruthy(model$fia_var))
      if (model$track == "Forest type") req(isTruthy(model$forest_type))
      if (model$track == "Forest group") req(isTruthy(model$forest_group))
      v <- ts_response(plot_occ(), model$track, model$fia_var,
                       model$forest_type, model$forest_group)
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
      if (!isTruthy(model$fia_var) || !model$fia_var %in% multi_vars) return(NULL)
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
      req(isTruthy(model$track))

      # In All mode "occupied" means any entity is present -- a forested
      # plot -- since there is no focal entity to be occupied by.
      #
      # Assigned before use for the same reason as display_values(): forcing
      # multi_values() inside an operator means a failed req() escapes from
      # dispatch as an error rather than cancelling.
      if (multi_mode()) {
            if (!isTruthy(model$fia_var) || !model$fia_var %in% multi_vars)
                  return(rep(FALSE, nrow(po)))
            y <- multi_values()
            return(rowSums(y > 0, na.rm = TRUE) > 0)
      }

      v <- if (model$track == "Species") {
            po$pres
      } else if (model$track == "Forest type") {
            req(isTruthy(model$forest_type))
            as.integer(!is.na(po$forest_type) &
                             po$forest_type == as.integer(model$forest_type))
      } else if (model$track == "Forest group") {
            req(isTruthy(model$forest_group))
            as.integer(!is.na(po$forest_group) &
                             po$forest_group == as.integer(model$forest_group))
      } else if (model$track == "All trees") {
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
            variable = if(model$track %in% c("Forest type", "Forest group")) "presence"
            else model$fia_var,
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
                         weight = dnorm(pd$env_dist, sd = model$theta_clim) *
                               dnorm(pd$geog_dist, sd = model$theta_geog))
      }

      out
})

# Per-entity responses for a community run, row-aligned with plot_occ()
# (R/plots.R)
multi_values <- reactive({
      req(isTruthy(model$track))
      if (!model$track %in% c("Forest type", "Forest group"))
            req(isTruthy(model$fia_var), model$fia_var %in% multi_vars)
      ts_multi_response(core_data, plot_occ(), pool_bbox(), model$track, model$fia_var)
})

# Per-plot numerator and denominator of a pooled-ratio variable, or NULL for
# the others (R/plots.R)
ratio_pieces <- reactive({
      req(isTruthy(model$track), isTruthy(model$fia_var))
      ts_ratio_pieces(plot_occ(), model$track, model$fia_var)
})

clim_scales <- core_data$clim_scales

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

# The unit list shows one collection at a time, filled from here, including at
# startup. It is searched server-side: collections run to thousands of units
# (counties). Switching collections clears the polygon, so the selector never
# shows a different unit than the map. The first fill selects the link's
# polygon, if it has one; the UI already opened on that polygon's collection.
polygon_from_link <- link$state$polygon_unit
observeEvent(input$polygon_collection, {
      updateSelectizeInput(session, "polygon_unit",
                           choices = boundary_choices(boundary_units, input$polygon_collection),
                           selected = polygon_from_link, server = TRUE)
      polygon_from_link <<- ""
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
      terra::is.related(terra::vect(g), terra::ext(c(model$bbox_x, model$bbox_y)), "intersects")
}

polygon_outside_msg <- "The focal polygon is outside the model bounding box. Use 'Snap box to polygon' under 'Focal polygon and site' to move the box to it."


## focal site
# Starts at the link's site, if it has one (server/url_state.R)
site <- reactiveValues(x = link$state$site[1], y = link$state$site[2])
