## How the model works

**Similar places, weighted.** For each grid cell, TreeShift gives every FIA plot a weight based on two things: how far away it is, and how different its climate is. Close plots with similar climate count almost fully; distant or climatically different plots count for little. The map value for a cell is a weighted summary of those plots.

**The two bandwidths.** The *geographic bandwidth* (in kilometers) sets how quickly weight falls off with distance. The *climate bandwidth* sets how quickly it falls off with climate difference, measured in standard deviations of each climate variable relative to variation across all FIA plots. Narrow bandwidths stay close to local plots but rest on fewer of them; wide bandwidths draw on more plots but blur local differences. There is no single correct setting, and it is worth checking whether a conclusion holds across a range of them. The *Model evaluation* panel helps with this.

**Future climate.** For a future scenario, each cell's projected climate is compared with the recent climate of FIA plots. The plots that best match that future climate, and are within the geographic bandwidth, shape the projection. This assumes forests respond to climate the way they do across today's landscape, and it does not model how quickly trees can actually move, establish, or die back.

**A generalization, not a copy of the plots.** Because each cell blends many plots, the map generalizes the plot data rather than exactly interpolating it. A cell can show a species as moderately likely even where the nearest plot lacks it.

**Methods.** The weighting is the analog impact model (AIM) implemented in the [`analogs`](https://matthewkling.github.io/analogs/) R package. Depending on the scale of your query, computations may be based on the exact analogs method or a nearly-exact "fast Fourier transform." Landscape connectivity uses circuit theory, solved with Circuitscape through the [`circuitscaper`](https://matthewkling.github.io/circuitscaper/) R package.
