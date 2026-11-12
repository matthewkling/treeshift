# Right sidebar: analysis panels that summarize the mapped result. Settings
# that change the map live in the left sidebar; the controls here only shape
# their own panel's figures, so none are carried in links.

ui_right_sidebar <- function() {
      sidebar(
            # title = "Focal site analogs",
            width = "40%",
            position = "right",
            open = TRUE,

            accordion(
                  open = FALSE,
                  multiple = FALSE,

                  ui_panel_community_composition(),
                  ui_panel_value_distributions(),
                  ui_panel_climate_relationships(),
                  ui_panel_model_evaluation(),
                  ui_panel_climate_models(),
                  ui_panel_focal_site_analogs()
            )

      )
}

      ## First panel: the map-level community measures
      ## reduce ACROSS entities, which is exactly the
      ## information this panel restores. A turnover or
      ## richness surface says how much changes; this
      ## says which entities changed, and is the only
      ## view in the app where entity identity survives
      ## the reduction.
ui_panel_community_composition <- function() {
      accordion_panel(
            "Community composition",
            helpText("When you select the 'community composition' FIA variable at left,
                                             the map summarizes across species or forest types into one summary number per cell.
                                             This module lets you unpack the individual species behind those summaries. The two settings
                                             below apply to all three plots."), hr(),

            ## Both controls apply to every plot below, so
            ## they sit together above the rule rather
            ## than beside whichever plot they happen to
            ## precede. "Entities shown" is the count
            ## NAMED rather than pooled into "other", and
            ## it feeds all three.
            ## The full set: unlike the other modules, a
            ## single cell is meaningful here -- reading
            ## the composition at one place is the point
            ## of drilling in.
            selectInput("comm_scope", "Summarize over", area_names) %>%
                  tooltip("Which part of the map these plots describe. All but the last average over many cells; 'Focal site' reads just the cell you clicked, so it shows what is at one place rather than a regional average."),
            sliderInput("comm_rank_n", "Entities shown", 5, 40, value = 15, step = 1) %>%
                  tooltip("How many species get named individually. The rest are combined into a single 'other' group rather than left out, so the totals stay right."),
            hr(),

            ## Leads the panel: it is the only view here
            ## that shows the community as a whole. Band
            ## thickness is one entity and total column
            ## height is the community sum, so a
            ## community that reshuffles without changing
            ## size looks different from one that shrinks.
            helpText("Current and predicted future forest composition. Each band is one species; the
                                             height of the whole column is the community."),
            plotOutput("comm_composition") %>%
                  tooltip("Each band runs from where a species is today, on the left, to where it is expected to be in the future, on the right. Labels sit on whichever side the band is thicker: on the left for species that decline, on the right for those that grow."),
            hr(),

            ## Second: the same numbers spread out, so
            ## rare-but-gaining stays distinct from
            ## abundant-but-declining. Absolute and
            ## proportional change rank differently and
            ## neither is wrong; the scatter shows both
            ## at once.
            helpText("A different view of the same data, showing whether gains are predicted in species that are already common or in ones
                                             that are currently rare. Rare and gaining sits at the top left; common and declining sits
                                             at the bottom right. Colored points are the species named above, gray are the rest."),
            plotOutput("comm_scatter"), hr(),

            ## Last: change without levels, for the cases
            ## where the ordering itself is the answer.
            ## Its "Rank by" control stays here because,
            ## unlike the two above, it applies to this
            ## plot alone.
            helpText("Species with the largest projected gains and losses."),
            selectInput("comm_rank_by", "Rank by",
                        c("absolute change" = "absolute",
                          "proportional change" = "proportional")) %>%
                  tooltip("Absolute change ranks by how much is gained or lost outright, so common species tend to lead. Proportional change measures that against how much is there now, so a rare species that doubles outranks a common one that shifts a little. Neither ordering is more correct than the other."),
            plotOutput("comm_ranking")
      )
}

      ## Named for its contents rather than its input:
      ## both plots are views of site_analogs(), the
      ## analog SET of the focal site. The site itself is
      ## now set under "Focal polygon and site" in the left sidebar.
      ##
      ## The composition plot that used to live here has
      ## been replaced by a stack-based one in the
      ## Community composition panel: it was the only view
      ## here where entity identity survived, and it
      ## belonged with the other cross-entity displays.
      ## Renamed from "Study area", the name the left panel
      ## that owned the polygon had at the time. This
      ## panel is named for what it renders, matching the
      ## rule that right-sidebar panels are local to one
      ## output.
ui_panel_value_distributions <- function() {
      accordion_panel(
            "Value distributions",
            helpText("This module summarizes how the predicted values on the map differ across geographic areas.
                                             It shows one curve (or bar stack, for categorical variables) for the whole bounding box,
                                             and one for each area you have set: your focal polygon, that polygon
                                             plus a surrounding buffer, and the neighborhood around your focal site."), hr(),

            plotOutput("polygon_density"),
            selectInput("density_transform", "X-axis scale",
                        choices = c("identity", "uniform", "sqrt", "log10"),
                        selected = "identity") %>%
                  tooltip("Optionally transform the x-axis for skewed variables. 'uniform' rank-transforms the values across the bounding box, spreading skewed distributions evenly.")
      )
}

ui_panel_climate_relationships <- function() {
      accordion_panel(
            "Climate relationships",
            helpText("This plot shows how the mapped variable changes along a climate gradient. The curve is a smoothed
                                             relationship fitted to either the FIA plot observations or the modeled grid cells; the
                                             density histogram above it shows how much data sits at each point along the gradient.
                                             For categorical variables, this becomes one distribution per species instead, showing where along the
                                             gradient each one comes out on top."), hr(),

            ## Two columns of two, paired by what they
            ## control: the left column defines the x
            ## axis, the right column defines which data
            ## the curve is fitted to.
            layout_column_wrap(
                  width = 1/2,
                  selectInput("cr_v1", "X-axis variable", display_vars()) %>%
                        tooltip("The variable the response is plotted against. It need not be one the model was fitted on; the plot says so when it isn't."),
                  selectInput("cr_data", "Dataset", c("modeled grid cells", "FIA plots")) %>%
                        tooltip("Whether to fit the curve to the model's predictions across the grid, or to the raw FIA plot observations."),
                  selectInput("cr_clim", "Climate timeframe", c("historical", "future")) %>%
                        tooltip("Which climate the x-axis is drawn from. Use 'future' to pair a future prediction with the future climate that produced it. FIA plot observations are historical, so pairing them with future climate relates an observed response to a climate the plot has not experienced."),
                  ## Same scope vocabulary as the evaluation and composition
                  ## modules, from the same named vector, so an area means the
                  ## same thing wherever it is offered.
                  selectInput("cr_scope", "Area", eval_scope_names) %>%
                        tooltip("Which part of the map the relationship is fitted to. Narrowing it shows how the variable behaves in your area specifically, which can differ from the regional pattern.")
            ),


            plotOutput("climate_response", height = "600px")
      )
}

ui_panel_model_evaluation <- function() {
      accordion_panel(
            "Model evaluation",
            helpText("This module uses 'leave-one-out cross validation' to quantify the accuracy of the model's present-day predictions.
                                             Each FIA plot is held out in turn, predicted from the others, and compared against what was actually measured there.
                                             The first plot below shows accuracy for the currently-displayed results, while the second section shows how accuracy
                                             varies over a range of alternative model parameters you could select.
                                             You can choose the spatial area ('evaluation scope') that accuracy is measured over."),
            br(),
            helpText("Note that these accuracy estimates and optimal parameter recommendations DO NOT apply to future predictions."),
            br(),

            selectInput("eval_scope", "Evaluation scope",
                        choices = eval_scope_names,
                        selected = "domain") %>%
                  tooltip(HTML("Which held-out plots the accuracy metrics are computed from. The model itself is unchanged.<br><br>
                                                       <dl style='margin:0;'>
                                                           <dt>Model bounding box</dt> <dd>Every plot in the analog pool.</dd>
                                                           <dt>Focal polygon</dt> <dd>Plots inside the focal polygon. Often few; check the reported n.</dd>
                                                           <dt>Buffered polygon</dt> <dd>Plots within 3 geographic bandwidths of the polygon &mdash; the neighborhood predictions there are actually drawn from. Usually the better choice for a small unit.</dd>
                                                           <dt>Focal site neighborhood</dt> <dd>Plots within the geographic cutoff of the clicked site.</dd>
                                                       </dl>
                                                       Scopes whose input is not yet set under 'Focal polygon and site' in the left sidebar will say so when selected.")),

            hr(),

            helpText("This plot compares model predictions against what was actually measured at each plot, using your current
                                             settings. Each red point is a group of plots with similar predicted values, plotted against what
                                             was really found there. Points on the dashed 1:1 line mean the predictions are unbiased, while
                                             deviations from the 1:1 line indicate bias; the shaded band is the IQR of observed values. Marks along the
                                             bottom show where the predictions actually fall. When you narrow the spatial scope, a gray line
                                             shows the same curve for plots in the model domain that are outside your scope."), br(),
            plotOutput("eval_scatter"),
            uiOutput("eval_scatter_stats"), hr(),

            helpText("This tool helps you explore how different parameter values affect predictive accuracy.
                                             For bandwidth parameters, it tests five values bracketing the value selected under 'Analog bandwidths'.
                                             For climate variables, it tests all combinations of the variables selected under 'Climate data'.
                                             The scatterplot above always reflects your current settings; the heatmap below shows the last completed
                                             run, and is faded and flagged if those settings have since changed."), br(),

            layout_column_wrap(
                  width = 1/2,
                  checkboxGroupInput("eval_params", NULL,
                                     choiceNames = c("Test geographic bandwidths", "Test climate bandwidths", "Test climate variables"),
                                     choiceValues = c("geog", "clim", "vars")),
                  ## Choices swap with the response type; declared with the
                  ## binary set because that matches the default variable.
                  selectInput("eval_stat", "Evaluation statistic", c("AUC", "TSS", "Brier")),
                  sliderInput("eval_n_vars", "Max n climate vars", 1, 6, value = c(1, 2), step = 1),
                  actionButton("optimal_params", "Apply optimal parameters") %>%
                        tooltip("Sets the bandwidths and climate variables under Model settings to the best-scoring combination for this area and statistic, and runs the model with them."),
                  actionButton("run_eval", "Run evaluations")
            ),
            ## The one control here that reaches back into the model
            ## settings, so it says so where it can be seen
            helpText("'Apply optimal parameters' changes your model settings and runs the model with them."),

            ## Directly above the heatmap, which is the
            ## only thing it describes. Sitting under the
            ## scope selector, it read as a caution about
            ## the scatterplot too -- and that plot always
            ## reflects current settings, so it is never
            ## the stale one.
            uiOutput("eval_stale"),
            plotOutput("eval_heatmap", height = "600px")
      )
}

ui_panel_climate_models <- function() {
      accordion_panel(
            "Climate models",
            helpText("This module helps you choose among, and contextualize, climate models.
                                             It can take a couple minutes to run.
                                             For every available emissions scenario and GCM, it plots the projected climate change
                                             (future minus historical) as a trajectory across time periods, so
                                             you can see which models are hot/dry, warm/wet, and so on, and how
                                             the spread fans out over eras and emissions scenarios. Models in your
                                             current ensemble selection (set in the left sidebar) are drawn at full
                                             opacity, others faded. The x and y variables can be changed
                                             without re-clicking 'Update'."), hr(),
            layout_column_wrap(
                  width = 1/2,
                  selectInput("model_delta_region", "Summarize over",
                              choices = area_names,
                              selected = "domain") %>%
                        tooltip("Region over which deltas are summarized. 'Model bounding box'
                                                        averages over the current bounding box; 'Focal site' samples
                                                        at the clicked site; 'Focal polygon' averages over the selected
                                                        polygon."),
                  div(style = "padding-top: 1.9rem;",
                      actionButton("model_delta_update", "Update chart") %>%
                            tooltip("Recompute the chart for the current axis variables and
                                                            summary region. The chart does not update automatically,
                                                            since computing deltas across the full model ensemble is
                                                            slow."))
            ),
            plotOutput("model_delta_plot"),
            layout_column_wrap(
                  width = 1/2,
                  selectInput("model_delta_x", "X variable",
                              choices = clim_vars()) %>%
                        tooltip("Climate variable for the x-axis (delta = future - historical)."),
                  selectInput("model_delta_y", "Y variable",
                              choices = clim_vars()) %>%
                        tooltip("Climate variable for the y-axis (delta = future - historical).")
            )
      )
}

      ## Last, and deliberately so. Every other panel here
      ## reports a result; this one shows the mechanism
      ## behind them -- which plots inform a prediction,
      ## and how the two bandwidths weight them. That is
      ## worth having for anyone asking where the numbers
      ## come from, but it is not what a manager opens
      ## first.
      ##
      ## Named for its subject rather than its purpose,
      ## even though the purpose is what the intro text
      ## explains. The analog type it shows is set with
      ## the focal site in the left sidebar, since it also
      ## drives the 'focal site similarity' raster and the
      ## 'focal site analogs' vector layer.
      ##
      ## The analog SET does not depend on the entity --
      ## the same plots are near and climatically similar
      ## whichever species is selected -- so only the
      ## point coloring is entity-specific, and it follows
      ## display_values() like the map markers do.
ui_panel_focal_site_analogs <- function() {
      accordion_panel(
            "Focal site analogs",
            helpText(
                  "This module helps you explore where a model prediction at a given location actually comes from:
                                          the 'analog' FIA plots that are geographically close and climatically similar.
                                          When you click the map to set a focal site, these figures show the FIA plots (drawn as points) that inform the estimate there.
                                          'Contemporary' analogs compare the focal cell's and plots' current climates, while 'reverse' and 'forward' analogs compare climates across time to make future projections.

                                          Reference lines mark the two bandwidths (dashed) and the hard cutoffs at three times each (solid).
                                          Weights, represented by shading in the righthand figure, fall off with both distances.

                                          Use the button below to see the same analog plots on the map.

                                          The analog type is set with the focal site, under 'Focal polygon and site' at left. It also
                                          governs two map layers: the 'focal site similarity' raster and the 'focal site
                                          analogs' points. Nothing here changes the model itself."), hr(),
            ## Kept, and relabeled. The old name
            ## described the mechanism rather than
            ## what turning it off would show. In a
            ## panel about how predictions are made,
            ## seeing which nearby plots the climate
            ## kernel rules out is half the point --
            ## otherwise the cutoff is a line with
            ## nothing outside it.
            input_switch("filter_clim", "Hide plots beyond the climate cutoff", TRUE) %>%
                  tooltip("Turn this off to also show the plots that are close by but too climatically different to count. Those are what the climate bandwidth excludes. Plots beyond three geographic bandwidths are never shown either way."),

            actionButton("highlight_site", "Show these plots on the map") %>%
                  tooltip("Switch the map's vector layer to this site's analog plots and zoom to its geographic cutoff, so the same plots can be seen in geographic space."),
            br(), br(),

            layout_column_wrap(
                  width = 1/2,
                  plotOutput("clim_plot"),
                  plotOutput("dist_plot")
            )
      )
}
