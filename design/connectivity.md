# TreeShift connectivity: design note

Status: agreed design, 2026-11-17, from a design discussion in claude.ai. Not yet implemented. This replaces the current experimental connectivity module (Circuitscape advanced mode with user-set electrodes, loss and scales). The underlying "dipole methods" are being built in circuitscaper (`dp_*()` functions; see circuitscaper's `dev/dipole/SPEC.md`).

## 1. Principles

- **Who it's for.** USFS planners. Connectivity is an explicit requirement of the 2012 Planning Rule, but planning direction for it tends to be vague. TreeShift's distinct contribution is species-specific, climate-aware connectivity for trees, now and under future climate.
- **One question per map, answerable in one sentence.** Method concepts (normalization, squared current, draws) stay hidden behind defaults.
- **One main setting: dispersal reach.** Everything else gets a default; a few expert options are tucked away.
- **Relative outputs.** Show classes or percentiles, never raw circuit values.
- **Labelled "experimental"** until the dipole-methods paper is at least a preprint.
- **Fast.** Seconds per map at the app's 2 km grain, reusing factorizations within a session where the cache key is airtight.

## 2. Inputs

| Symbol | Meaning | Source |
|---|---|---|
| O_i | Present occupancy at cell i | Baseline `response` surface. It is geographically local because AIM uses the θ_geog kernel to plots, so it reflects where the species actually is, not just where the climate suits it. |
| S_j | Future suitability at cell j | Future `response` surface: occurrence among nearby plots whose present climate matches j's future climate. **Mask where effective weight is low** (future climate has no nearby analogs) and show those cells as a "locally novel future climate" class rather than scoring them. |
| K(d) | Dispersal kernel (reach) | New setting, **independent of θ_geog**. θ_geog is the model's spatial scale for the niche, not dispersal. |
| E(a, b) | Climate-match kernel | For the genetic variants only. It compares a source's present climate a with a destination's future climate b. Bandwidth default to be decided; θ_clim is one candidate. |
| M_j | Climate departure at j | 1 − E(present climate at j, future climate at j), for genetic rescue. |
| g | Circuit conductance | Climate suitability: √(present × future suitability), as now. Optional human-modification toggle: × (1 − HM)^k, with HM from the Global Human Modification index v3 (2022, 300 m), averaged to 2 km. Consider excluding the timber-harvest stressor if stressor layers are available, since harvested forest is largely permeable to trees. |
| c | Propagule supply | Expert setting, used in the arrival probability 1 − exp(−c·P). It changes how fast scores saturate, not their ranking. |

Long-distance dispersal: trees can skip barriers narrower than their long-distance dispersal, which cell-to-cell circuits can't represent. For v1, choose the grain and conductance floor so narrow unsuitable gaps are passable and wide ones remain barriers, and say so in the help text. Better models (sparsified long-range links, stepping-stone graphs) are research items in the circuitscaper project's `dev-ideas.md`.

## 3. Maps

The module offers a selector of maps. All share the reach setting.

### A. Current connectivity (within era)

- **Question:** where is movement among the current range concentrated, and where is it fragile?
- **Method:** random-dipole map (`dp_omni`), with pair weight O_i · O_j · K(d_ij) and conductance g.
- **Default view:** normalized current in four classes: impeded, diffuse, intensified, channelized. Omniscape's vocabulary, but the thresholds must be recalibrated for our normalized distribution, not copied.
- **Pinch-point option:** power dissipation, I²·r per link summed per cell, shown as top-percentile classes. This is the sensitivity of kernel-weighted total resistance to a proportional loss of a link's permeability: where degradation would hurt connectivity most.

### B. Range-shift variants

Each variant scores a focal site as a **state term** at that site times a **flow term** summed over partner cells. The same pair weights then drive a **routes** view.

| Variant | Focal | State term | Flow term P | Score shown | Reads as |
|---|---|---|---|---|---|
| Expansion centers | future cell j | (1 − O_j) · S_j | Σ_i O_i K(d_ij) | (1 − O_j) · S_j · (1 − e^(−cP)) | Likely new population nuclei: absent now, suitable later, with propagules arriving |
| Contraction | present cell i | O_i · (1 − S_i) | Σ_j S_j K(d_ij) | Shown inverted, as **stranded populations**: O_i · (1 − S_i) · e^(−cP) | Declining populations with no future habitat in reach |
| Genetic rescue | cell j | O_j · M_j | Σ_i O_i K(d_ij) E(c_i, c_j) | O_j · M_j · (1 − e^(−cP)) | At-risk populations that natural gene flow from pre-adapted sources could rescue |
| Genetic escape | cell i | O_i | Σ_j S_j K(d_ij) E(c_i, c_j) | Shown inverted, as **stranded genotypes**: O_i · e^(−cP) | Locally adapted populations whose climate will exist nowhere in reach |

In the table, S_i is future suitability at the source cell i; c_i is the source's present climate and c_j the destination's future climate.

- **Flow terms.** These are `analog_density()` from the analogs package: present grid to future grid for the inbound terms (expansion, rescue), future grid to present grid for the outbound terms (contraction, escape). Use raw sums (`normalize = FALSE`) before the arrival transform.
- **Planner reading.** The inverted forms (stranded populations and genotypes) point at places where assisted migration or ex situ conservation may be needed. They are probably the most actionable outputs.
- **Routes toggle (per variant).** Run the random-dipole pathways with that variant's full pair weight: source term × destination term × K (× E). For example, expansion uses O_i · (1 − O_j) S_j · K. The density map and the routes map are then two views of one model: whether propagules get there, and which routes the flow takes. Pinch points are available here too.
  - **Barriers.** The flow terms ignore barriers, which is defensible for trees given long-distance dispersal; the routes see them.
  - **Dependency.** Climate-matched pair weights for the genetic variants' routes need a pair-weight hook in circuitscaper's `dp_pathways` (deferred there as "continuous climate-similarity pair weights").

## 4. Settings

- **Main:** map selector (current connectivity; expansion; contraction; rescue; escape); routes on/off; pinch points on/off; dispersal reach as presets named by dispersal syndrome (e.g. wind-dispersed, animal-cached, gravity), with values grounded in literature and a default that could follow the species' dispersal mode; human-modification on/off.
- **Expert:** reach in km; propagule supply c; climate-match bandwidth; conductance option (suitability, uniform, perhaps wind).

## 5. What to retire from the current module

- **Advanced-mode electrodes:** source and target choices (grid cells, FIA plots, selected site).
- **Loss settings:** loss layer, loss scale and conductance scale.
- **Outputs:** current on target, current loss, voltage, and the loss input layer.
- **Code:** `server/connectivity.R`'s cs_advanced solve and its correction loop. (Its "none" loss setting also leaks about 27% of current through a 10⁻⁶ ground floor, one more reason to retire it.)
- **Decide:** current direction, plus the animation and field vector layers. The design has no default directional output. Net flux (one oriented solve) could feed these as an expert option, but it has no dispersal limit; keep it only if the animation earns its explanatory cost.

## 6. Implementation notes

- **Dependencies.** Circuit maps need circuitscaper's `dp_omni` and `dp_pathways` (in development). Flow terms need only `analog_density()`, so the four density maps can ship first, before the circuit routes exist.
- **Factorization reuse.** Factor once per box, grain and conductance option; reuse across variants, reach settings and the uniform flow-potential circuit. Per CLAUDE.md, cache only with an airtight key.
- **Masking.** Exclude low-effective-weight future cells from S everywhere: in scores, as partners in flow terms and as dipole endpoints.
- **Scale.** At the default box and 2 km, expect about 5 s per circuit map, less with cached factorizations. Density maps are fast (an FFT path is available in analogs).

## 7. Open questions

1. Reach preset values, and whether the default follows a species' dispersal syndrome; this needs a short literature pass.
2. Climate-match bandwidth for E and the departure term M; reuse θ_clim or set separately?
3. Default for c, and whether to expose it at all.
4. Class thresholds for normalized current and pinch points (recalibrate on a few species).
5. Effective-weight threshold for masking S.
6. Keep wind conductance?
7. Keep any directional output (see §5)?
8. Human-modification exponent k, and the timber stressor question.

## 8. Checks before release

- **Sanity checks:** expansion scores near zero where occupancy is high; stranded populations only where present habitat has little future habitat in reach; rescue near zero where the climate barely changes.
- **Earlier analyses:** compare against the NAFEW `analog_density()` maps for the same species (similar, not identical: these definitions are revised from first principles).
- **Planner test:** show the maps with one-line captions to one or two planners. If they can say what they'd do differently because of a map, it stays.
