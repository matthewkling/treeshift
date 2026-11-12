## Uncertainty

Every map value is an estimate, and is subject to several sources of uncertainty. TreeShift doesn't roll them into one error bar, because they differ in kind: no single number can capture the uncertainty from the evidence near a site, the model's accuracy, and future climate together. Instead, it gives you a way to look at each one.

**How much evidence is behind a value.** Three map layers, under `Raster variable`, describe the FIA plots behind each cell's value:

- `Analog support`: how much evidence there is, counted as the number of exactly matching plots it's worth. Ten closely matching plots (geographically close and climatically similar) count for about 10; ten weak matches count for much less.
- `Analog closeness`: how well those plots match the cell's climate and location, from near 0 to 1. Low support with high closeness means a few good matches; low closeness means only weak ones.
- `Analog spread`: how much those plots' values differ from one another, in the units of the variable. High spread means the plots disagree, so the value averages over very different stands. For presence variables, spread is highest where about half the plots have the species. Not shown for ratio variables.

Low support or low closeness means thin evidence. These layers matter most for future projections: where a projected climate has few close matches on today's landscape (a "non-analog" climate), support and closeness drop, and the projection should be read with extra caution. The layers depend on the bandwidth settings, so compare places within one map, not across settings. They describe the evidence, not the size of the error.

**How well the model reproduces today's forests.** The `Model evaluation` panel tests the model against the FIA data using leave-one-out cross-validation. It holds out each plot in turn, estimates it from the other plots, and compares the estimate with what was measured there. Presence variables are scored by how well the model separates plots where the species was found from plots where it wasn't (AUC, TSS and Brier score). Other variables are scored by how close the estimates come to the measured values (R² and RMSE). You can score the whole study area or just the plots in and around a polygon or focal site, and compare a range of bandwidths and climate variables to see which settings fit the data best.

These scores apply to the present day only. Future forests can't be measured yet, so there is no direct test of a projection, and a model that reproduces today's forests well can still be wrong about the future.

**Climate models and scenarios.** A future projection depends on the emissions pathway, the time period, and the climate models selected under `Climate data`. When several climate models are selected, their climates are averaged before the forest model runs, so the map shows a single projection for the average climate, not the range across models. To see how much a result depends on that choice, run the climate models one at a time, or compare emissions pathways, and look for conclusions that hold across them. The `Climate models` panel shows how the models differ in their projected climate change for your area, for example which are warmer, drier or wetter, so you can choose models that span the range.

**What these tools don't cover.** None of these measures capture the processes the model leaves out (see [Reading the results](#about-reading_results)), or the assumption that forests will respond to future climate the way they respond across today's landscape.

**Putting it together.** A result deserves more confidence when the evidence behind it is strong (good support and closeness), the model scores well for your area, and the result holds across reasonable bandwidths, climate models and emissions pathways.
