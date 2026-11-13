# Left sidebar: model and map settings, one accordion panel per group.
#
# The panels come in two groups, each under a heading: model settings, which
# take effect when the Run model button under that heading is pressed, and
# display settings, which apply right away. The button keys off
# model_setting_ids (helpers/model_settings.R), not panel membership; keep
# the two in step when adding controls.
#
# `s` is the starting state (helpers/url_state.R): the app defaults, or the
# settings in the link the app was opened with. Initial values and dependent
# menus are built from it, so a linked page opens already consistent rather
# than starting at defaults and being updated into place.

ui_left_sidebar <- function(s) {
      sidebar(
            open = "always",
            width = 320,

            ui_sidebar_group("Model settings", "sliders"),

            ## Applies the model settings (server/model_settings.R). Starts
            ## up to date: the page was built from the settings the model
            ## runs with.
            div(class = "ts-update",
                actionButton("update_model", "Model is up to date", class = "btn-dark w-100",
                             disabled = TRUE),
                uiOutput("pending_note")),

            accordion(
                  id = "model_acc",
                  open = FALSE,
                  multiple = FALSE,

                  ui_panel_fia_variable(s),
                  ui_panel_raster_grid(s),
                  ui_panel_climate_data(s),
                  ui_panel_interpolation(s),
                  ui_panel_connectivity(s)
            ),

            ui_sidebar_group("Display settings", "eye"),
            accordion(
                  id = "display_acc",
                  open = FALSE,
                  multiple = FALSE,

                  ui_panel_map_layers(s),
                  ui_panel_focal_places(s),
                  ui_panel_cartography(),
                  ui_panel_export()
            )
      )
}

# Heading over a group of panels
ui_sidebar_group <- function(title, icon_name) {
      div(class = "ts-group", tags$h2(class = "ts-group-title", icon(icon_name), title))
}

ui_panel_fia_variable <- function(s) {
      accordion_panel(
            "FIA variable",
            ## Labels lowercased to match the other parameter menus;
            ## values are unchanged, since switch() statements
            ## throughout the server match on them.
            selectInput("track", "Entity",
                        choices = track_choices,
                        selected = s$track) %>%
                  tooltip(HTML("Choose the focal entity for the model.<br><br>
                                          <dl style='margin:0;'>
                                              <dt>Species</dt> <dd>Model the presence or abundance of a single focal tree species.</dd>
                                              <dt>Forest type</dt> <dd>Model the probability that a site is classified as a focal FIA forest type (a fine-grained classification).</dd>
                                              <dt>Forest group</dt> <dd>Model the probability that a site is classified within a focal FIA forest type group (a coarser classification).</dd>
                                              <dt>All trees</dt> <dd>Model community-wide forest properties (any-species presence or total basal area).</dd>
                                          </dl>")),
            ## How MANY entities, orthogonal to which quantity is
            ## measured on each. Unlabelled, so it reads as a
            ## qualifier on the entity dropdown above rather than as
            ## a separate axis. Option text is track-specific ("single
            ## species" vs "single forest type") and set by an
            ## observer; the underlying values stay single/multi.
            ##
            ## Hidden for All trees, which already aggregates the
            ## community into a single response BEFORE interpolation
            ## -- the opposite operation to interpolating each entity
            ## and reducing afterward.
            conditionalPanel(
                  condition = "input.track != 'All trees'",
                  selectInput("entity_scope", NULL,
                              choices = entity_scope_choices(s$track),
                              selected = s$entity_scope) %>%
                        tooltip(HTML("How many entities (species or forest types) to model at once.<br><br>
                                                <dl style='margin:0;'>
                                                    <dt>Single</dt> <dd>Model one focal species, forest type, or group.</dd>
                                                    <dt>Community composition</dt> <dd>Model every entity in the region at once, reducing results to a community-level summary.</dd>
                                                </dl>
                                                The life stage and variable settings apply either way."))
            ),
            ## Single-entity pickers are hidden rather than removed, so
            ## input$species keeps its value while community mode is
            ## active.
            conditionalPanel(
                  condition = "input.track == 'Species' && input.entity_scope != 'multi'",
                  selectizeInput("species", NULL, choices = spp, selected = s$species) %>%
                        tooltip("Select a focal species by scrolling or typing")
            ),
            conditionalPanel(
                  condition = "input.track == 'Forest type' && input.entity_scope != 'multi'",
                  selectizeInput("forest_type", NULL,
                                 choices = forest_type_choices,
                                 selected = s$forest_type) %>%
                        tooltip("Select a focal FIA forest type by scrolling or typing. Labels show the FORTYPCD followed by the type name.")
            ),
            conditionalPanel(
                  condition = "input.track == 'Forest group' && input.entity_scope != 'multi'",
                  selectizeInput("forest_group", NULL,
                                 choices = forest_group_choices,
                                 selected = s$forest_group) %>%
                        tooltip("Select a focal FIA forest type group by scrolling or typing.")
            ),
            ## Life stage closes the entity group -- kind, scope,
            ## which one, which individuals -- all four unlabelled so
            ## they read as a single cluster under one heading. The
            ## scope row poses a question the picker row answers
            ## ("single species" / "which species"), so stage trails
            ## rather than splitting that pair.
            ##
            ## That placement inverts the old dependency: stage used
            ## to be filtered by metric, and now metric is filtered by
            ## stage (metric_choices takes it). A control whose
            ## options are set by a control below it reads badly, and
            ## the inverted direction is the more truthful one --
            ## "seedling basal area" is unavailable because FIA does
            ## not break basal area out by size class, not because
            ## basal area lacks a stage.
            ##
            ## Forest type and group have no stage: their entity is a
            ## plot classification, and multi_values() builds an
            ## indicator matrix that ignores the variable entirely.
            conditionalPanel(
                  condition = "input.track == 'Species' || input.track == 'All trees'",
                  selectInput("stage", NULL,
                              stage_choices,
                              selected = s$stage) %>%
                        tooltip(HTML("Which tree size class to model. Variables that FIA does not break out by size class &mdash; basal area, and its proportion &mdash; are only offered under 'all life stages'.<br><br>
                                                <dl style='margin:0;'>
                                                    <dt>Seedlings</dt> <dd>&lt; 1 in. DBH/DRC.</dd>
                                                    <dt>Saplings</dt> <dd>1.0&ndash;4.9 in. DBH/DRC.</dd>
                                                    <dt>Adults</dt> <dd>&ge; 5.0 in. DBH/DRC.</dd>
                                                    <dt>All life stages</dt> <dd>Combined.</dd>
                                                </dl>"))
            ),
            ## Declared from metric_choices() rather than a literal
            ## list. The syncing observer uses ignoreInit, so this is
            ## the startup state -- and a hard-coded list once held
            ## "abundance", which is not valid at the default life
            ## stage of "any". Selecting it made resolve_var() return
            ## NULL and the app silently do nothing.
            selectInput("metric", "Variable",
                        choices = metric_choices(s$track, s$stage),
                        selected = s$metric) %>%
                  tooltip(HTML("Choose the quantity to model.<br><br>
                                          <dl style='margin:0;'>
                                              <dt>Occurrence</dt> <dd>Probability the focal entity is present within a single FIA plot footprint (~0.07 ha) randomly located within the grid cell.</dd>
                                              <dt>Abundance</dt> <dd>Population density of the selected life stage, in stems per acre.</dd>
                                              <dt>Basal area</dt> <dd>Live tree and sapling basal area, in ft&sup2; per acre.</dd>
                                              <dt>Species richness</dt> <dd>Number of tree species present within a single FIA plot footprint (~0.07 ha).</dd>
                                              <dt>Proportion basal area</dt> <dd>The focal species' share of total stand basal area, multiplied by the probability that stand basal area is nonzero.</dd>
                                              <dt>Regeneration fraction</dt> <dd>Probability that the selected regeneration stage is present relative to adults, conditional on the species occurring on the plot. Calculated as p(regen) / (p(regen) + p(adult)). Adults only = 0; regen only = 1; adults and regen together = intermediate value (not precisely 0.5 due to differing sample areas).</dd>
                                          </dl>")),
            ## Notes on community mode: run cost, and variables that
            ## have no community form
            conditionalPanel(
                  condition = "input.track != 'All trees' && input.entity_scope == 'multi'",
                  uiOutput("comm_note"),
                  uiOutput("comm_cost_note")
            ),

            ## Canonical variable resolved from (track, metric, stage);
            ## hidden because all downstream logic switches on input$fia_var.
            tags$div(style = "display: none;",
                     selectInput("fia_var", NULL,
                                 c("presence probability", "seedling presence",
                                   "sapling presence", "tree presence",
                                   "proportion basal area", "total basal area",
                                   "number of trees", "number of saplings",
                                   "number of seedlings",
                                   "species richness", "seedling species richness",
                                   "sapling species richness", "adult species richness",
                                   "seedling fraction", "sapling fraction"),
                                 selected = s$fia_var))
      )
}

      ## The grid the model is computed on: its resolution and extent.
      ## Grain comes first because it is the resolution everything is
      ## evaluated at, including the summaries over the focal polygon and
      ## the focal-site similarity surface.
ui_panel_raster_grid <- function(s) {
      accordion_panel(
            "Raster grid",
            helpText("The grid the model is computed on: its cell size, and the rectangular area it covers."), hr(),

            sliderInput("clim_agg", "Raster grain size", url_bounds$clim_agg[1], url_bounds$clim_agg[2],
                        value = s$clim_agg, step = 1, post = " km") %>%
                  tooltip("Cell size of the model grid, and of every summary drawn from it. Larger values run faster but are less precise."),

            hr(),

            tags$strong("Model bounding box"),
            helpText("The rectangular region the model is computed over. A smaller box runs faster."),
            sliderInput("bbox_x", "Longitude range",
                        min = clim_bbox["xmin"], max = clim_bbox["xmax"],
                        value = s$bbox_x),
            sliderInput("bbox_y", "Latitude range",
                        min = clim_bbox["ymin"], max = clim_bbox["ymax"],
                        value = s$bbox_y),
            layout_column_wrap(
                  width = 1/2,
                  uiOutput("snap_focal_btn"),
                  actionButton("window_domain", "Snap to current view")
            )
      )
}

      ## Places to compare against the result: an area summarized
      ## within and a point queried at. Display settings, not model
      ## ones: they apply right away and don't change the grid. The
      ## exception is 'Snap box to polygon', which sets the bounding box
      ## under Raster grid, like 'Apply optimal parameters' in the
      ## evaluation panel. (A focal site used as a connectivity
      ## electrode does feed the model; it is left live for now.)
ui_panel_focal_places <- function(s) {
      accordion_panel(
            "Focal polygon and site",
            helpText("Places to compare against the model results: an optional focal polygon to summarize within,
                      and an optional focal site to query."), hr(),

            tags$strong("Focal polygon"),
            helpText("A management unit to summarize results within. Optional."),
            selectInput("polygon_collection", NULL,
                        choices = boundary_collections(boundary_units),
                        selected = boundary_collection_of(boundary_units, s$polygon_unit)) %>%
                  tooltip("The kind of unit to choose from: national park units, national forests and grasslands, ranger districts, BLM field offices, wilderness areas, states, counties, or EPA ecoregions. Changing it clears the current polygon."),
            # Choices are filled by the server (server/data.R)
            selectizeInput("polygon_unit", NULL, choices = NULL,
                           # Room for the whole of the largest collection, plus
                           # "None", so the list can be browsed without typing
                           options = list(placeholder = "Search by name",
                                          maxOptions = max(table(boundary_units$collection)) + 1)) %>%
                  tooltip("Type to search by name. Agency boundaries are administrative and can include land the agency doesn't manage: private inholdings, and for BLM field offices, all land within the office's area."),
            layout_column_wrap(
                  width = 1/2,
                  actionButton("polygon_snap", "Snap box to polygon") %>%
                        tooltip("Set the model bounding box (under 'Raster grid') to cover the polygon and its buffer, so summaries over them use the whole area. Press 'Run model' to apply it."),
                  actionButton("polygon_clear", "Clear") %>%
                        tooltip("Remove the polygon from the map and plots.")
            ),

            hr(),

            tags$strong("Focal site"),
            helpText("A single location to query. Click the map to set one. Optional."),
            uiOutput("site_status"),
            ## Which climates the site is matched on. Not a model setting:
            ## it applies right away, to the 'Focal site analogs' panel and
            ## the two focal site map layers.
            selectInput("analog_direction", "Analog type", c("Contemporary analogs", "Reverse analogs", "Forward analogs"),
                        selected = s$analog_direction) %>%
                  tooltip(HTML("Which climates are being matched:
                                          <ul>
                                              <li>Contemporary: plots whose present climate resembles the site's present climate. This is what the historical prediction is built from.</li>
                                              <li>Reverse (inbound): plots whose present climate resembles the site's future climate. This is the basis for future model predictions in the app.</li>
                                              <li>Forward (outbound): plots whose future climate resembles the site's present climate. Where the site's current conditions are heading; interesting, but not currently used in this app's core modeling machinery.</li>
                                          </ul>")),
            ## This group is for SETTING the site; showing which plots
            ## inform it is the 'Focal site analogs' panel's subject.
            actionButton("site_clear", "Clear") %>%
                  tooltip("Remove the focal site from the map and plots.")
      )
}

ui_panel_climate_data <- function(s) {
      accordion_panel(
            "Climate data",
            helpText("Select climate variables and future scenarios."), hr(),
            selectInput("clim_vars", "Climate variables", clim_vars(), selected = s$clim_vars, multiple = TRUE) %>%
                  tooltip(HTML("Choose one or more variables for measuring multivariate climate similarity:<br><br>
                                          <dl style='margin:0;'>
                                              <dt>CWD</dt> <dd>Climatic water deficit</dd>
                                              <dt>AET</dt> <dd>Actual evapotranspiration</dd>
                                              <dt>PPT</dt> <dd>Total annual precipitation (log-transformed)</dd>
                                              <dt>TMEAN</dt> <dd>Mean annual temperature</dd>
                                              <dt>TMINCM</dt> <dd>Avg min temp of coldest month</dd>
                                              <dt>TMAXWM</dt> <dd>Avg max temp of warmest month</dd>
                                          </dl>")),
            selectInput("scenario_ssp", "Emissions pathway",
                        choices = ssp_choices,
                        selected = s$scenario_ssp) %>%
                  tooltip("Future emissions scenario (Shared Socioeconomic Pathway). Higher numbers correspond to higher emissions."),
            selectInput("scenario_period", "Future time period",
                        choices = period_choices,
                        selected = s$scenario_period) %>%
                  tooltip("Future time period."),
            selectizeInput("scenario_model", "Climate model(s)",
                           choices = model_choices,
                           selected = s$scenario_model,
                           multiple = TRUE,
                           options = list(plugins = list("remove_button"))) %>%
                  tooltip("Select one or more global climate models (GCMs). When multiple are selected, climate layers are averaged across the model ensemble.
                                      See the 'Climate models' panel in the right sidebar for help comparing GCMs and choosing among them.")
      )
}

ui_panel_interpolation <- function(s) {
      accordion_panel(
            "Interpolation",
            helpText("Predictions for each grid cell are built from nearby FIA plots with similar climate. Plots are weighted
                                 more heavily the closer they are, in both senses, and these two settings control how quickly that
                                 falloff happens. They are standrd deviations of Gaussian weighting kernels.
                                 Plots more than three times either setting away are not used at all."), hr(),
            sliderInput("theta_geog", "Geographic bandwidth (km)", url_bounds$theta_geog[1], url_bounds$theta_geog[2],
                        value = s$theta_geog) %>%
                  tooltip("Controls how plot weight decays with geographic distance from a grid cell.
                                      This distance also defines two of the areas used in some of the righthand modules: the 'focal site neighborhood' is everything within three bandwidths of the site you clicked, and the 'buffered polygon' is your polygon widened by the same amount."),
            sliderInput("theta_clim", "Climate bandwidth (z)", url_bounds$theta_clim[1], url_bounds$theta_clim[2],
                        value = s$theta_clim) %>%
                  tooltip("Measured in standard deviations of climate across all FIA plots in the contiguous US, so 0.1 means plots within a tenth of a standard deviation count almost fully.")
      )
}

# Connectivity settings are model settings: they wait for Run model like the
# rest (model_setting_ids). The switch itself is not one -- it only offers the
# connectivity layers in the Map layers menus -- so it stays out of the
# "what changed" list.
ui_panel_connectivity <- function(s) {
      accordion_panel(
            "Landscape connectivity",
            helpText("Model circuit-based landscape connectivity between historical (source) and future (target) FIA surfaces.
                      This approach can identify key movement corridors and pinch points for large-scale climate connectivity.
                      Turn it on to show its settings and to add its layers to the 'Map layers' menus under Display settings.
                      To view connectivity results, choose 'connectivity' as the raster variable there.
                      The first connectivity result may take a couple of minutes while the solver starts up; later ones are faster."), hr(),
            uiOutput("cs_note"),
            conditionalPanel(
                  condition = "input.track == 'All trees' || input.entity_scope != 'multi'",
                  input_switch("cs_enable", "Enable connectivity", s$cs_enable) %>%
                        tooltip("Adds the connectivity raster and current animation/field vector layers
                                 to the 'Map layers' menus under Display settings.")
            ),
            conditionalPanel(
                  condition = "input.cs_enable && (input.track == 'All trees' || input.entity_scope != 'multi')",
                  selectInput("cs_sources", "Origin electrodes (sources)", c("grid cells", "FIA plots", "selected site"),
                              selected = s$cs_sources) %>%
                        tooltip("Sites where current is injected. If `grid cells` or `FIA plots` is selected, current is injected
                                            in proportion to modeled or observed (respectively) values of the selected FIA variable.
                                            If `selected site` is chosen, current is injected only at the focal site."),
                  selectInput("cs_targets", "Destination electrodes (target grounds)", c("grid cells", "FIA plots", "selected site"),
                              selected = s$cs_targets) %>%
                        tooltip("Grounding sites that 'pull' current and allow it to exit the circuit.
                                            If `grid cells` or `FIA plots` is selected, ground conductance is proportional to modeled values
                                            of the selected FIA variable. If `selected site` is chosen, the focal site."),
                  selectInput("cs_conductance_layer", "Horizontal conductance", c("suitability", "wind", "suitability * wind", "uniform"),
                              selected = s$cs_conductance_layer) %>%
                        tooltip("The inverse of landscape resistance, horizontal conductance determines how easily current can flow across a site.
                                            `Suitability` models conductance proportional to the mean of historical and future FIA raster values.
                                            `Wind` models conductance proportional to the frequency and speed with which wind blows toward targets."),
                  selectInput("cs_loss_layer", "Loss (non-target grounds)", c("none", "uniform", "suitability" = "mean suitability"),
                              selected = s$cs_loss_layer) %>%
                        tooltip("Choose whether and how to model current 'leakage' representing dispersers that are removed from the pool
                                            (via mortality or permanent settlement) before reaching target grounds."),
                  sliderInput("cs_conductance_scale", "log10 conductance scale", url_bounds$cs_conductance_scale[1],
                              url_bounds$cs_conductance_scale[2], value = s$cs_conductance_scale, step = .1) %>%
                        tooltip("Relative strength of horizontal conductance, compared to source, target, and loss components."),
                  sliderInput("cs_loss_scale", "log10 loss scale", url_bounds$cs_loss_scale[1],
                              url_bounds$cs_loss_scale[2], value = s$cs_loss_scale, step = .1) %>%
                        tooltip("Relative strength of loss conductance, compared to source, target, and horizontal components.
                                                        Ignored if loss is set to 'none.'")
            )
      )
}

ui_panel_map_layers <- function(s) {
      accordion_panel(
            "Map layers",
            helpText("Choose which model results are drawn on the map."), hr(),
            ## Connectivity entries are offered only when it is on; see
            ## the cs_enable observer in server/controls.R
            selectInput("stat", "Raster variable",
                        if (url_cs_on(s)) stat_choices_cs else stat_choices_base,
                        selected = s$stat) %>%
                  tooltip(HTML("<dl style='margin:0;'>
                                              <dt>FIA variable</dt> <dd>Predicted values of the variable selected above.</dd>
                                              <dt>Analog evidence</dt> <dd>Measures of the FIA evidence behind each cell's estimate: its support, closeness, or spread. Choose one in the menu that appears below.</dd>
                                              <dt>Focal site similarity</dt> <dd>How climatically similar each cell is to the site you clicked on the map. Which climates are compared is set by the analog type, under 'Focal polygon and site'.</dd>
                                              <dt>Connectivity</dt> <dd>Circuit-based connectivity outputs between historical and future surfaces (available when connectivity is turned on under 'Landscape connectivity'); choose the specific connectivity variable below</dd>
                                          </dl>")),
            ## The analog evidence entry's measure; see evidence_choices
            ## (helpers/menus.R)
            conditionalPanel(
                  condition = "input.stat == 'evidence'",
                  selectInput("evidence_stat", "Evidence measure", unname(evidence_choices),
                              selected = s$evidence_stat) %>%
                        tooltip(HTML("<dl style='margin:0;'>
                                              <dt>Support</dt> <dd>How much FIA evidence is behind each cell's estimate, counted as the number of exactly matching plots it's worth. A plot that closely matches the cell's climate and location counts almost fully; a weaker match counts for less. Low values mean thin or distant evidence.</dd>
                                              <dt>Closeness</dt> <dd>How well the plots behind each estimate match the cell's climate and location, from near 0 (only weak matches) to 1 (exact matches). Low support with high closeness means a few good matches; low support with low closeness means only weak ones.</dd>
                                              <dt>Spread</dt> <dd>How much the values at the plots behind each estimate differ from one another, in the units of the variable. High spread means those plots disagree, so the estimate averages over very different stands. For presence variables, spread is highest where about half the plots have the species. Not shown for ratio variables.</dd>
                                              <dd>Support, closeness and spread depend on the bandwidth settings, so compare places within one map, not across settings.</dd>
                                          </dl>"))
            ),
            uiOutput("stat_note"),
            conditionalPanel(
                  condition = "input.stat == 'connectivity'",
                  selectInput("cs_stat", "Connectivity variable", cs_stat_names, selected = s$cs_stat) %>%
                        tooltip(HTML("<dl style='margin:0;'>
                                              <dt>Current flow</dt> <dd>Total current passing through each cell; high values indicate important range shift corridors</dd>
                                              <dt>Current on target</dt> <dd>Current delivered to destination grounds; shows where dispersers successfully reach future-suitable habitat</dd>
                                              <dt>Current loss</dt> <dd>Current lost to non-target grounds; represents dispersers removed via mortality or settlement before reaching destinations</dd>
                                              <dt>Current pinch points</dt> <dd>Power dissipation (I²R); highlights bottlenecks where connectivity is both concentrated and tenuous — priorities for restoration or protection</dd>
                                              <dt>Current direction</dt> <dd>Bearing of net current flow through each cell</dd>
                                              <dt>Voltage</dt> <dd>How connected each cell is to source populations; high values indicate greater accessibility for dispersers</dd>
                                              <dt>Conductance</dt> <dd>Model input: ease of movement across the landscape (inverse of resistance)</dd>
                                              <dt>Loss</dt> <dd>Model input: rate at which dispersers are removed from the migrating pool at each cell</dd>
                                          </dl>"))
            ),
            ## Timeframe comes first among these because it is the
            ## first cut on the community measures below: "state or
            ## change" narrows the list before the user reads it.
            ##
            ## "change" covers both a difference of a static measure
            ## (change in richness) and the paired measures, which are
            ## change quantities in their own right and not deltas of
            ## anything -- net richness change of zero is entirely
            ## compatible with complete species replacement.
            selectInput("time", "Timeframe",
                        c("historical", "future", "change" = "delta"),
                        selected = s$time) %>%
                  tooltip("Choose whether to display historical (2000-2025) or future (scenario selected under 'Climate data') model predictions, or the change between them.
                                      Only relevant for the `FIA variable` raster variable and the analog support, closeness and spread layers."),
            ## Community mode only: how the per-entity surfaces are
            ## combined into one map. It sits below the timeframe,
            ## which narrows the list, and is filtered by the
            ## variable (comm_measure_choices()).
            conditionalPanel(
                  condition = "input.track != 'All trees' && input.entity_scope == 'multi'",
                  ## Outside community mode the variable may have no
                  ## community measures (the ratio variables); the menu
                  ## still needs a valid value, so it falls back to the
                  ## occurrence set.
                  selectInput("comm_measure", "Community measure",
                              choices = if (length(comm_measure_choices(s$track, s$fia_var, s$time)))
                                    comm_measure_choices(s$track, s$fia_var, s$time)
                              else comm_measure_choices(s$track, "presence probability", s$time),
                              selected = s$comm_measure) %>%
                        tooltip(HTML("How the per-entity (species or forest types) surfaces are combined into one map. Only combinations that are well defined for the current entity type and variable are listed.<br><br>
                                                <dl style='margin:0;'>
                                                    <dt>Richness</dt> <dd>Expected number of entities present, summed across per-entity occurrence probabilities. No threshold: a cell with 0.5 probability for each of two species has expected richness 1.</dd>
                                                    <dt>Diversity</dt> <dd>Effective number of entities (Hill N1 = exp of Shannon entropy), in the same units as richness. Low where one entity dominates, however many rare ones are present.</dd>
                                                    <dt>Evenness</dt> <dd>Pielou&#39;s J (Shannon entropy relative to its maximum given the species pool). 0 where a single species dominates, with higher values for more balanced communities.</dd>
                                                    <dt>Dominance</dt> <dd>Share of the community held by its single largest member (0&ndash;1).</dd>
                                                    <dt>Dominant entity</dt> <dd>Which entity is has the highest value at each cell. Drawn as discrete classes; all but the most frequent few are pooled into &quot;other&quot;. Most intuitive over an abundance variable &mdash; over occurrence it shows the most widely probable entity.</dd>
                                                </dl>
                                                The measures below compare the historical and future communities directly (they are not deltas of the measures above).<br>
                                                <dl style='margin:0;'>
                                                    <dt>Turnover</dt> <dd>Bray&ndash;Curtis dissimilarity between the historical and future communities. 0 is identical composition, 1 is no overlap.</dd>
                                                    <dt>Expected gains / losses</dt> <dd>How much arrives and how much departs, in the variable&#39;s own units, counted separately. Asymmetric on purpose &mdash; what arrives is a different question from what leaves.</dd>
                                                    <dt>Colonization demand</dt> <dd>Gains as a share of the future community: the fraction of what will be suitable there that is not there now.</dd>
                                                    <dt>Loss fraction</dt> <dd>The opposite: losses as a share of the community that is there now. Read alongside turnover, the two fractions show which side is driving it &mdash; high turnover with a low loss fraction is accumulation rather than replacement.</dd>
                                                    <dt>Biggest winner / loser</dt> <dd>Which single entity gains or loses the most at each cell, drawn as discrete classes. Says which entity is moving, where the measures above say how much. Use the Community composition panel to see the full ranking behind any one cell.</dd>
                                                </dl>"))
            ),
            ## Number of categories drawn for the categorical
            ## measures. 143 entities is far past what any palette or
            ## legend can carry, so all but the most frequent N are
            ## lumped into a single "other" class. The right N differs
            ## a lot between a polygon-scoped domain and a
            ## region-wide one, so it is a control rather than a
            ## constant.
            conditionalPanel(
                  condition = paste0(
                        "input.track != 'All trees' && input.entity_scope == 'multi' && ",
                        "['", paste(categorical_measures, collapse = "','"),
                        "'].indexOf(input.comm_measure) >= 0"),
                  numericInput("comm_top_n", "Categories shown", value = s$comm_top_n,
                               min = url_bounds$comm_top_n[1], max = url_bounds$comm_top_n[2], step = 1) %>%
                        tooltip("How many entities get their own color. The rest are pooled into 'other', ranked by how much of the mapped area each one covers.")
            ),
            ## Support floor for ratio-valued layers. Shown only when
            ## the selected layer is a ratio -- the community
            ## measures other than richness, and the regeneration
            ## fractions. Richness IS the denominator, so masking it
            ## by itself would be circular.
            conditionalPanel(
                  condition = "(input.track != 'All trees' && input.entity_scope == 'multi' && input.comm_measure != 'richness' && input.comm_measure != 'expected gains' && input.comm_measure != 'expected losses') || input.metric == 'regeneration fraction'",
                  numericInput("support_floor", "Minimum expected species", value = s$support_floor, min = 0, step = 0.5) %>%
                        tooltip(HTML("Blank out cells where there is too little predicted mass present for the ratio to mean anything.<br><br>
                                                Ratios stay well behaved as their denominator approaches zero but stop being informative: evenness and the regeneration fractions take their most extreme values where the community or species is barely present.<br><br>
                                                The label shows which quantity the threshold applies to &mdash; expected richness for occurrence variables, the summed response for abundance ones, and the species' own occurrence probability for the regeneration fractions. Set 0 for no masking.<br><br>
                                                This also rescales the color palette, since the extreme values from near-empty cells otherwise compress the rest of the map into a narrow band."))
            ),
            selectInput("vector", "Vector variable",
                        if (url_cs_on(s)) vector_choices_cs else vector_choices_base,
                        selected = s$vector) %>%
                  tooltip(HTML("Select which point or line layer is displayed:<br><br>
                                          <dl style='margin:0;'>
                                              <dt>Occupied FIA plots</dt> <dd>Plots where the selected species was recorded</dd>
                                              <dt>All FIA plots</dt> <dd>All plots including non-forested plots not surveyed</dd>
                                              <dt>Focal site analogs</dt> <dd>Plots within climatic and geographic range of the selected site. Which climates are matched &mdash; contemporary, reverse or forward &mdash; is set by the analog type, under 'Focal polygon and site'.</dd>
                                              <dt>Current animation</dt> <dd>Animation with particle movement showing direction and magnitude of mean current flow (available when connectivity is turned on)</dd>
                                              <dt>Current field</dt> <dd>Vector field showing showing direction and magnitude of mean current flow (available when connectivity is turned on)</dd>
                                          </dl>"))
      )
}

ui_panel_cartography <- function() {
      accordion_panel(
            "Cartography",
            helpText("Adjust the styling of the map."), hr(),
            sliderInput("vector_scale", "Vector point/line size", -1, 1, 0, step = .1) %>%
                  tooltip("Use this slider to adjust the size of FIA plot markers or the length of current lines"),
            sliderInput("rast_opacity", "Raster opacity", 0, 1, .8, step = .01),
            selectInput("rast_transform", "Scale transformation",
                        choices = c("identity", "log10", "sqrt", "uniform")) %>%
                  tooltip("Choose how raster values are spread across the colors. The legend stays in the variable's own units; under 'uniform', colors are evenly spread across the map's values, so legend ticks fall at uneven values. Not applied to analog support, which is shown in fixed classes."),
            ## Choices swap between continuous ramps and qualitative
            ## sets depending on the layer, rather than showing both
            ## lists at once and relying on the user to pick the
            ## applicable one.
            selectInput("rast_palette", "Color palette",
                        choices = rast_palettes_continuous,
                        selected = "default") %>%
                  tooltip("Override the default palette for the current raster variable. When the raster shows the historical or future FIA variable, the FIA plot markers and focal-site scatterplots adopt the same palette. For layers showing entity identities the options become qualitative palettes."),
            input_switch("show_site", "Show focal site", TRUE),
            uiOutput("site_color_ui"),
            input_switch("show_polygon", "Show focal polygon", TRUE),

            ## Hidden until a polygon is selected -- an outline color
            ## for nothing is just clutter.
            uiOutput("polygon_color_ui")
      )
}

ui_panel_export <- function() {
      accordion_panel(
            "Export",
            helpText("Share a link to these settings, or download the raster data displayed on the map as a .tif,
                                  the underlying FIA plot data as a .csv, or a record of the settings and versions behind them as a .txt."), hr(),
            ## The address bar always holds the link to the current settings
            ## (server/url_state.R); this copies it. The clipboard call is
            ## browser-side, and reports back so the server can confirm it.
            actionButton("copy_link", "Copy link to these settings", icon = icon("link"),
                         onclick = paste0(
                               "navigator.clipboard.writeText(window.location.href).then(",
                               "() => Shiny.setInputValue('copy_link_done', 'ok', {priority: 'event'}), ",
                               "() => Shiny.setInputValue('copy_link_done', 'failed', {priority: 'event'}));")) %>%
                  tooltip("The link reopens the app with these model settings and map layers, so others can see the same results.
                           It works with this version of TreeShift and its data; results can change when either is updated.
                           Your browser's address bar holds the same link."),
            br(), br(),
            textInput("download_filename", "Filename (without extension)", "fia_aim_result"),
            downloadButton("download_raster", "Download raster"),
            br(),
            downloadButton("download_points", "Download FIA points") %>%
                  tooltip(HTML("Plot-level FIA data for the points currently drawn on the map.<br><br>
                                           <dl style='margin:0;'>
                                               <dt>x, y</dt> <dd>Plot location (publicly released, fuzzed ~1.6 km)</dd>
                                               <dt>entity</dt> <dd>Selected species, forest type, or forest group</dd>
                                               <dt>variable</dt> <dd>Selected FIA variable</dd>
                                               <dt>value</dt> <dd>Observed historical value of that variable</dd>
                                               <dt>occupied</dt> <dd>Whether the entity is present at the plot</dd>
                                               <dt>sampled</dt> <dd>FALSE at plots that were never surveyed</dd>
                                           </dl>
                                           The analog layer adds climate/geographic distance from the focal site and
                                           the analog's kernel weight.")),
            br(),
            downloadButton("download_record", "Download settings record") %>%
                  tooltip("A text file for citing these results: the link that reopens them, the settings in plain language,
                           and the versions of TreeShift and its data that produced them.")
      )
}
