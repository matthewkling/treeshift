library(data.table)
library(tidyverse)
library(janitor)
library(terra)
select <- dplyr::select


# load FIA tables -------------------------------

dir <- "/Volumes/ultramarine/FIA/2025_06_27"

plot <- fread(paste0(dir, "/CSV_FIADB_ENTIRE/ENTIRE_PLOT.csv"), stringsAsFactors = F,
              colClasses = c(CN = "character")) %>%
      select(PLT_CN = CN, MANUAL, DESIGNCD, PLOT_STATUS_CD, SAMP_METHOD_CD,
             STATECD, INVYR, MEASYEAR, UNITCD, COUNTYCD, PLOT,
             LAT, LON, ELEV) %>% as_tibble()

# SUBPLOT table
# primary key: CN
# unique key: PLT_CN, SUBP
# When PLOT.MANUAL <1.0, the field crew measured slope at a condition level
# When PLOT.MANUAL >=1.0, slope is collected at a subplot level, and then the aspect
# from the subplot representing the greatest proportion of the condition is assigned as a surrogate
subplot <- fread(paste0(dir, "/CSV_FIADB_ENTIRE/ENTIRE_SUBPLOT.csv"), stringsAsFactors = F,
                 colClasses = c(CN = "character", PLT_CN = "character")) %>%
      select(PLT_CN, SUBP, SUBP_CN = CN, PREV_SBP_CN,
             SUBP_STATUS_CD, POINT_NONSAMPLE_REASN_CD,
             SLOPE, ASPECT) %>%
      as_tibble()

# TREE table
# primary key: CN
# unique key: PLT_CN, SUBP, TREE
# Tree numbers (TREE) can be used to track trees when PLOT.DESIGNCD is the same between inventories.
tree <- fread(paste0(dir, "/CSV_FIADB_ENTIRE/ENTIRE_TREE.csv"),
              stringsAsFactors = F,
              select = c("PLT_CN", "PLOT", "SUBP", "INVYR", "TREE", "SPCD",
                         "PREVCOND", "CONDID",
                         "DIA", "DIAHTCD", "STATUSCD", "DIACHECK",
                         "RECONCILECD", "TPA_UNADJ"
              ),
              colClasses = c(PLT_CN = "character", DIST = "numeric")) %>%
      as_tibble()

# SEEDLING table
# foreign key: PLT_CN (links to CN from plot table)
seedling <- fread(paste0(dir, "/CSV_FIADB_ENTIRE/ENTIRE_SEEDLING.csv"), stringsAsFactors = F,
                  colClasses = c(CN = "character", PLT_CN = "character")) %>%
      select(PLT_CN, SUBP, SPCD, seed_tpa = TPA_UNADJ) %>%
      as_tibble() %>%
      # SEEDLING can carry multiple rows per subplot-species (one per condition);
      # collapse to one row so the join below doesn't multiply seedling density.
      group_by(PLT_CN, SUBP, SPCD) %>%
      summarize(seed_tpa = sum(seed_tpa, na.rm = TRUE), .groups = "drop")

# SPECIES table (WOODLAND flag = species measured at diameter-at-root-collar, not DBH)
species <- fread(paste0(dir, "/FIADB_REFERENCE/REF_SPECIES.csv"), stringsAsFactors = F) %>%
      select(SPCD, COMMON_NAME, GENUS, SPECIES, WOODLAND) %>%
      as_tibble()

# CONDITION table (used to identify and exclude planted plots, and for forest type)
# primary key: CN
# unique key: PLT_CN, CONDID
cond <- fread(paste0(dir, "/CSV_FIADB_ENTIRE/ENTIRE_COND.csv"), stringsAsFactors = F,
              colClasses = c(CN = "character", PLT_CN = "character")) %>%
      select(PLT_CN, CONDID, PROP_BASIS, STDORGCD, PHYSCLCD, COND_STATUS_CD, CONDPROP_UNADJ,
             INVYR, FORTYPCD, FLDTYPCD, CONDPROP_UNADJ,
             SLOPE_COND = SLOPE, ASPECT_COND = ASPECT) %>%
      as_tibble()
planted <- cond %>%
      filter(STDORGCD == 1) %>%
      select(PLT_CN) %>%
      distinct() %>%
      mutate(planted = TRUE)


# forest types -------------------------------

forest_type <- cond %>%
      group_by(PLT_CN) %>%
      filter(INVYR == max(INVYR)) %>% # latest survey only
      slice_max(CONDPROP_UNADJ, n = 1, with_ties = FALSE) %>% # dominant condition only (deterministic on ties -> one row per plot)
      ungroup() %>%
      select(PLT_CN, forest_type = FORTYPCD)
type_group <- read_csv(paste0(dir, "/FIADB_REFERENCE/REF_FOREST_TYPE_GROUP.csv")) %>%
      select(forest_group = VALUE, forest_group_name = MEANING, forest_group_abbr = ABBR)
type <- read_csv(paste0(dir, "/FIADB_REFERENCE/REF_FOREST_TYPE.csv")) %>%
      select(forest_type = VALUE, forest_type_name = MEANING, forest_group = TYPGRPCD) %>%
      left_join(type_group)
write_csv(type, "data/forest_type_codes.csv")


# munge subplot topography -------------------------------

# subplot <- subplot %>%
#       rename(slope_pct = SLOPE,
#              aspect_deg = ASPECT) %>%
#       mutate(# slopes below 5% have aspect set to 360; set these to zero to avoid bias
#             slope_pct = ifelse(slope_pct < 5, 0, slope_pct),
#             slope_rad = atan(slope_pct/100),
#             slope_deg = slope_rad * 180 / pi,
#             aspect_rad = aspect_deg / 180 * pi) %>%
#       mutate(northness = cos(aspect_rad) * sin(slope_rad),
#              eastness = sin(aspect_rad) * sin(slope_rad)) %>%
#       select(PLT_CN:POINT_NONSAMPLE_REASN_CD, northness, eastness)


# join tables -----------------------------------------

# Tree and seedling records are attached to the plot scaffold *independently*,
# rather than threading seedlings through the tree table. Joining them to each
# other (on PLT_CN-SUBP-SPCD) has three failure modes: a seedling record for a
# species with no >=1" stem on that subplot has no tree row to match, so it was
# appended with NA plot attributes and silently dropped by the plot filters
# (losing exactly the regeneration-without-adults signal); a seedling riding on
# dead-tree rows was dropped by the live-stem filter; and a seedling value was
# replicated across each of that species' tree rows, inflating the sum. Keeping
# the two paths separate removes all three by construction.

# Plot scaffold: every plot-level filter applied once, and the latest visit
# selected per plot location. Keyed by plt_cn.
plot_base <- plot %>%
      left_join(planted, by = "PLT_CN") %>%
      left_join(forest_type, by = "PLT_CN") %>%
      filter(is.na(planted),
             PLOT_STATUS_CD %in% 1:2, # keep sampled plots (forest or nonforest); drop non-sampled (status 3)
             DESIGNCD %in% c(1, 111:118, 230:242, 311:323, 328, 501:506), # 4-subplot design
             MANUAL >= 2.0) %>%
      clean_names() %>%
      mutate(lon = ifelse(lon > 0, lon - 360, lon),
             yr = measyear, # measurement year (not invyr) is the correct sampling date
             plot_id = paste(statecd, unitcd, countycd, plot)) %>%
      filter(lat > 24, lat < 50,
             is.finite(lon), is.finite(lat)) %>%
      # Most recent visit per plot location. Grouped by plot_id (not subplot) so
      # each plot is represented by exactly one visit, for every life stage.
      group_by(plot_id) %>%
      filter(yr == max(yr)) %>%
      ungroup() %>%
      select(plt_cn, plot_id, lon, lat, forest_type, samp_method_cd, yr)

# Species lookup: clean binomial + woodland flag, resolved once. Records not
# identified to species (missing genus/epithet, or genus-level / unknown codes
# like "spp.") get a NA name so they don't surface as fake species -- their stems
# still count toward all-trees totals.
species_lut <- species %>%
      mutate(species_name = ifelse(is.na(GENUS) | GENUS == "" |
                                         is.na(SPECIES) | SPECIES == "" |
                                         tolower(SPECIES) %in% c("spp.", "spp", "unknown"),
                                   NA_character_, paste(GENUS, SPECIES)),
             is_woodland = !is.na(WOODLAND) & toupper(WOODLAND) == "Y") %>%
      select(SPCD, genus = GENUS, species = species_name, is_woodland)

# Live stems (>=1" DBH/DRC), attached to their plot by plt_cn.
stems <- tree %>%
      filter(STATUSCD == 1, SUBP <= 4) %>% # live trees on the 4 subplots
      select(PLT_CN, SUBP, SPCD, DIA, TPA_UNADJ) %>%
      left_join(species_lut, by = "SPCD") %>%
      clean_names() %>%
      inner_join(plot_base, by = "plt_cn")

# Seedlings (<1" DBH/DRC), attached to their plot by plt_cn -- independently of
# the tree records, so a species regenerating where it has no established stem
# survives.
seeds <- seedling %>%
      filter(SUBP <= 4) %>%
      left_join(species_lut, by = "SPCD") %>%
      clean_names() %>%
      inner_join(plot_base, by = "plt_cn")


# summarize to plot-species level -------------------

# Woodland species (junipers, pinyons, dryland oaks) are measured at
# diameter-at-root-collar, not DBH, which overstates a DBH-basis basal
# area. Convert DRC -> DBH (inches) via Chojnacky & Rogers 1999. Their
# equations are fit for pinyon, Gambel oak, and Utah juniper; we extend
# them by genus (the growth form driving the DRC/DBH offset is broadly
# genus-conserved): Pinus -> pinyon slope, Quercus -> oak slope, all other
# woodland -> juniper baseline. Juniper is the fallback because the
# remaining woodland genera (mesquite, mountain-mahogany, etc.) share its
# low-branching, multi-stemmed growth form -- the largest DRC/DBH offset of
# the three -- more than the single-stem pinyon or oak forms. Only
# woodland-flagged species are touched, so timber pines/oaks (DBH-measured)
# are unaffected. stm assumed multi-stem (0); floored at 0 where the linear
# model goes negative for small stems. Used for BA only -- the sapling/tree
# split below stays on the field-measured diameter.
stems <- stems %>%
      mutate(dia_ba = case_when(
            !is_woodland       ~ dia,
            genus == "Pinus"   ~ pmax(-1.9374 + 0.9823 * dia, 0), # pinyons
            genus == "Quercus" ~ pmax(-1.4599 + 0.9533 * dia, 0), # dryland/scrub oaks
            TRUE               ~ pmax(-2.6843 + 1.0222 * dia, 0))) # junipers + other woodland

tree_summ <- stems %>%
      group_by(plt_cn, species) %>%
      summarize(basal_area = sum(0.0054542 * dia_ba^2 * tpa_unadj, na.rm = T), # ft2/acre, all live stems (DBH basis; woodland DRC->DBH converted)
                n_saplings = sum(tpa_unadj[dia < 5], na.rm = T),  # stems/acre, saplings 1.0-4.9" DBH (microplot)
                n_trees = sum(tpa_unadj[dia >= 5], na.rm = T),    # stems/acre, trees >=5.0" DBH (subplot/macroplot)
                .groups = "drop")

# Summed over subplots directly from the seedling records -- no tree rows to
# replicate across, so no de-duplication needed.
seed_summ <- seeds %>%
      group_by(plt_cn, species) %>%
      summarize(n_seedlings = sum(seed_tpa, na.rm = T), # stems/acre, seedlings <1" DBH (microplot)
                .groups = "drop")

# Combine the two life-stage paths at plot-species level. A species present in
# only one of them (adults but no regen, or regen but no adults) gets 0 for the
# other, rather than being dropped.
plot_species <- full_join(tree_summ, seed_summ, by = c("plt_cn", "species")) %>%
      mutate(across(c(basal_area, n_saplings, n_trees, n_seedlings),
                    ~ ifelse(is.na(.x), 0, .x)))

# Attach plot attributes. Plots with no tree or seedling records survive as a
# single species = NA row of zeros -- these are the non-forest absences.
d <- plot_base %>%
      left_join(plot_species, by = "plt_cn") %>%
      mutate(across(c(basal_area, n_saplings, n_trees, n_seedlings),
                    ~ ifelse(is.na(.x), 0, .x))) %>%
      select(lon, lat, forest_type, species, samp_method_cd,
             basal_area, n_saplings, n_trees, n_seedlings)

# note: if we keep topo-interpolation downstream, may want to retain subplot
# resolution at this stage -- subplot attributes would attach to `stems` /
# `seeds` by (plt_cn, subp) before the summaries above.



# add topographic variables for nonforest plots ---------------------------------------------

# note: skipping topo for now

# # topography variables
# lf_slope_deg <- rast("/Volumes/T7/landfire/LF2020_SlpD_CONUS/Tif/LF2020_SlpD_CONUS.tif")
# lf_aspect_deg <- rast("/Volumes/T7/landfire/LF2020_Asp_CONUS/Tif/LF2020_Asp_CONUS.tif")
# NAflag(lf_slope_deg) <- -9999
# NAflag(lf_aspect_deg) <- -9999
#
# # plot locations for extraction
dp <- vect(select(d, lon, lat))
crs(dp) <- crs(rast("data/clim_historical.tif"))
# dp <- project(dp, crs(lf_slope_deg))
#
# # extract and munge
# topo <- c(lf_slope_deg, lf_aspect_deg) %>%
#       extract(dp, ID = FALSE) %>%
#       setNames(c("slope", "aspect")) %>%
#       mutate(slope = slope / 180 * pi, # degrees to radians
#              aspect = aspect / 180 * pi,
#              northness = cos(aspect) * sin(slope),
#              eastness = sin(aspect) * sin(slope)) %>%
#       select(northness_lf = northness,
#              eastness_lf = eastness)
#
# # merge and reconcile
# d <- d %>%
#       bind_cols(topo) %>%
#       mutate(northness = ifelse(is.na(northness), northness_lf, northness),
#              eastness = ifelse(is.na(eastness), eastness_lf, eastness)) %>%
#       select(-northness_lf, -eastness_lf)


# remove non-sampled plots with non-natural vegetation or with forest -----------------------

# load and extract landfire existing vegetation type
# veg <- rast("/Volumes/T7/landfire/LF2024_EVT_250_CONUS/Tif/LC24_EVT_250.tif") # load raster
veg <- rast("/Volumes/ultramarine/landfire/LF2020_BPS_220_CONUS/Tif/LC20_BPS_220.tif") # load raster
dp <- project(dp, crs(veg))
if(crs(veg) != crs(dp)) stop()
veg <- raster::extract(veg, dp)

# filter on veg
nonforest_veg <- c("Barren-Rock/Sand/Clay", "Grassland", "Perennial Ice/Snow", "Riparian",
                   "Savanna", "Shrubland", "Sparse")
nonforest <- foreign::read.dbf("/Volumes/ultramarine/landfire/LF2020_BPS_220_CONUS/Tif/LC20_BPS_220.tif.vat.dbf") %>%
      as_tibble() %>%
      filter(GROUPVEG %in% nonforest_veg)
d <- d %>%
      filter(samp_method_cd == 1 |
                   veg$BPS_CODE %in% unique(nonforest$BPS_CODE))


# merge with climate -------------------------------

# load climate data
climate <- rast("data/clim_historical.tif")

# filter to climate data bbox
bb <- ext(climate)
dc <- d %>%
      ungroup() %>%
      rename(x = lon, y = lat) %>%
      filter(between(x, bb$xmin, bb$xmax),
             between(y, bb$ymin, bb$ymax))

# filter to plots with non-missing climate data
dc <- cbind(dc, terra::extract(climate, dc[,1:2]))
dc <- filter(dc, is.finite(tmean))

# drop climate columns (they're re-extracted in app)
d <- dc %>% select(x:n_seedlings)


# export -------------------------------------------

# round to reduce size on disk
d_small <- d %>%
      # select(-ID) %>%
      mutate(across(
            where(is.double),
            ~ round(.x, 3)
      ))

write_csv(d_small, "data/fia.csv")
