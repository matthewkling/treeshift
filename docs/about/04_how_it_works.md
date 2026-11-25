## How the model works

### Similar places, weighted

TreeShift predictions are based on "spatial climate analogs." For each grid cell, it gives every FIA plot a weight based on two things: how far away it is, and how different its climate is. Close plots with similar climate count almost fully; distant or climatically different plots count little. A cell's mapped value is a weighted summary of those plots. Because each cell blends many plots, the map generalizes the plot data rather than exactly interpolating it; a cell can show a species as moderately likely even where the nearest plot lacks it.

The weighting is the analog impact model (AIM) implemented in the [`analogs`](https://matthewkling.github.io/analogs/) R package. Depending on the scale of your query, computations may be based on the exact analogs method or a nearly-exact "fast Fourier transform." 

### The two bandwidths

The *geographic bandwidth* (in kilometers) sets how quickly weight falls off with distance. The *climate bandwidth* sets how quickly it falls off with climate difference, measured in standard deviations of each climate variable relative to variation across all FIA plots. Narrow bandwidths stay close to local plots but rest on fewer of them; wide bandwidths draw on more plots but blur local differences. There is no single correct setting, and it is worth checking whether a conclusion holds across a range of them. The *Model evaluation* panel shows which settings fit the FIA data best, and the *Model sensitivity* panel shows how much the estimate at a site changes across them.

### Future climate

For a future scenario, each cell's projected climate is compared with the recent climate of FIA plots. The plots that best match that future climate, and are within the geographic bandwidth, shape the projection. This assumes forests respond to climate the way they do across today's landscape, and it does not model how quickly trees can actually move, establish, or die back.


### Landscape connectivity

TreeShift's landscape connectivity models run as a second step on top of the gridded FIA maps. They estimate how the landscape constrains movement routes between mapped habitat patches. You can model omnidirectional connectivity within a single period's habitat by setting `timeframe` to *historical* or *future*, or you can model directional pathways linking historical and future habitats by setting `timeframe` to  *change*.

Connectivity is modeled with circuit theory, which treats the landscape like an electrical circuit: movement is easiest through places the model rates as good habitat, and hardest through poor habitat or open water. The result is a measure of "current flow" through each cell, which represents the cumulative number of movement routes that pass through a given cell, summing over all possible routes between pairs of cells on the landscape. Cell pairs are weighted by their mapped FIA value (so that the model reflects movement between suitable habitats) and are included only if they are within `dispersal reach` of each other (letting you model various scales of dispersal potential). Ease of movement through cells, i.e. resistance, is based on climate suitability, human land use, or a combination.

The connectivity maps compare the movement through each place with what it would be on an even landscape, with the same habitat and reach but every place equally easy to cross (a.k.a. "normalized current"). Values well above 1 mark places where movement is funneled: corridors and pinch points where movement squeezes through a narrow gap or neck of habitat. Values near 1 mean movement spreads out as it would on an even landscape, and values below 1 mean it is impeded, so it goes around. Like the other maps, connectivity stops at the US borders: routes through Canada or Mexico are not counted, so near a border movement may look more confined than it actually is.

Two habitat measures answer a different question: not where movement goes, but how much habitat each place of habitat can reach. Both are areas, in square kilometers. *Proximity of suitable habitat* is the area of habitat within dispersal reach, whatever lies between. *Connected habitat* counts that habitat in full where it is as easy to reach as it would be through continuous, fully suitable habitat on unmodified land, and less where it is harder to reach: it is lower where there is less habitat within reach, where that habitat is farther away, or where the land between is hard to cross. The gap between the two shows how much the land in between matters. A longer dispersal reach takes in more habitat, so both give larger areas; compare maps made with the same reach. Under *change*, each has two directions: the future habitat around each place of historical habitat (could its population follow its shifting habitat?), and the historical habitat around each place of future habitat (could it be colonized?). Connected habitat is estimated from a limited number of random samples, so a single cell's value is uncertain by roughly a fifth; read the map's broad classes rather than individual cells. Both measures need a dispersal reach other than *no limit*.

Landscape connectivity is estimated with the random-dipole method of the [`circuitscaper`](https://matthewkling.github.io/circuitscaper/) R package. Large or fine-grained maps can take a minute or more.
