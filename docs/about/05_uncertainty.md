## Uncertainty

Every map value is an estimate, and is subject to several sources of uncertainty. TreeShift doesn't roll them into one error bar, because they differ in kind: no single number can capture the uncertainty from the evidence near a site, the model's accuracy, and future climate together. Instead, it gives you a way to look at each one.

### How much evidence is behind a value

Three map layers, under `Raster variable`, describe the FIA plots behind each cell's value:

- `Analog support`: how much evidence there is, counted as the number of exactly matching plots it's worth. Ten closely matching plots (geographically close and climatically similar) count for about 10; ten weak matches count for much less.
- `Analog closeness`: how well those plots match the cell's climate and location, from near 0 to 1. Low support with high closeness means a few good matches; low closeness means only weak ones.
- `Analog spread`: how much those plots' values differ from one another, in the units of the variable. High spread means the plots disagree, so the value averages over very different stands. For presence variables, spread is highest where about half the plots have the species. Not shown for ratio variables.

Low support or low closeness means thin evidence. These layers matter most for future projections: where a projected climate has few close matches on today's landscape (a "non-analog" climate), support and closeness drop, and the projection should be read with extra caution. The layers depend on the bandwidth settings, so compare places within one map, not across settings. They describe the evidence, not the size of the error.

### How well the model reproduces today's forests

The `Model evaluation` panel tests the model against the FIA data using leave-one-out cross-validation. It holds out each plot in turn, estimates it from the other plots, and compares the estimate with what was measured there. Presence variables are scored by how well the model separates plots where the species was found from plots where it wasn't (AUC, TSS and Brier score). Other variables are scored by how close the estimates come to the measured values (R² and RMSE). You can score the whole study area or just the plots in and around a polygon or focal site, and compare a range of bandwidths and climate variables to see which settings fit the data best.

These scores apply to the present day only. Future forests can't be measured yet, so there is no direct test of a projection, and a model that reproduces today's forests well can still be wrong about the future.

### Climate models and scenarios

A future projection depends on the emissions pathway, the time period, and the climate model selected under `Climate data`. Each run uses a single climate model, and different models project different futures. You can pick a model by name, or by how much drying it projects for the bounding box: the least, median or most drying of the available models, ranked by their change in climatic water deficit (a measure of drought stress). The median is the default. To see how much a result depends on the choice of model, run the least and most drying models as well, or compare emissions pathways, and look for conclusions that hold across them. For a single site, the `Model sensitivity` panel shows this without rerunning the map. The `Climate futures` panel shows how the models differ in their projected climate change for your area, for example which are warmer, drier or wetter.

### How much a site's estimate depends on your choices

The `Model sensitivity` panel takes the estimate at a focal site and varies one choice at a time, holding the others at your current settings: the climate model, the emissions pathway, each bandwidth, the climate variables (swapping one variable at a time for another), and the sample of FIA plots near the site (by resampling them). A wide row means the estimate there depends a lot on that choice; a narrow one means it holds up. Points are colored as on the map and sized by their `Analog support`, so you can see when an alternative rests on little evidence. The rows aren't on a common scale, so the panel doesn't add them into one range, and it doesn't show how the choices interact. The plot-sample range is somewhat narrow, because nearby plots aren't independent of one another. It is a quick way to check one place without rerunning the whole map for each alternative.

### What these tools don't cover

None of these measures capture the processes the model leaves out (see [Reading the results](#about-reading_results)), or the assumption that forests will respond to future climate the way they respond across today's landscape.

### Putting it together

A result deserves more confidence when the evidence behind it is strong (good support and closeness), the model scores well for your area, and the result holds across reasonable bandwidths, climate models and emissions pathways.
