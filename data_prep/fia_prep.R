library(data.table)
library(tidyverse)
library(janitor)
library(terra)
select <- dplyr::select


# load FIA tables -------------------------------

dir <- "/Volumes/T7/FIA/2025_06_27"

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

d <- plot %>%
      left_join(planted) %>%
      left_join(forest_type) %>%
      left_join(subplot) %>%
      left_join(tree) %>%
      full_join(seedling) %>%
      left_join(species) %>%
      mutate(PLOT_ID = paste(STATECD, UNITCD, COUNTYCD, PLOT),
             SUBPLOT_ID = paste(PLOT_ID, SUBP)) %>%

      filter(is.na(planted),
             PLOT_STATUS_CD %in% 1:2, # keep sampled plots (forest or nonforest); drop non-sampled (status 3)
             DESIGNCD %in% c(1, 111:118, 230:242, 311:323, 328, 501:506), # 4-subplot design
             SUBP <= 4) %>%

      clean_names() %>%
      mutate(lon = ifelse(lon>0, lon-360, lon)) %>%
      filter(lat > 24) %>%

      mutate(tree_id = paste(subplot_id, designcd, tree)) %>%
      filter(manual >= 2.0,
             is.finite(lon), is.finite(lat)) %>%

      mutate(yr = measyear) %>% # measurement year (not invyr) is the correct sampling date
      select(-invyr, -measyear)


# summarize to plot-species level -------------------

d <- d %>%
      filter(lat < 50) %>%
      group_by(plot_id, subplot_id) %>%
      filter(yr == max(yr),
             statuscd == 1 | is.na(statuscd) # live trees, or no trees
      ) %>%
      # Build a clean binomial; NA-out records that aren't identified to species
      # (missing genus/epithet, or genus-level / unknown codes like "spp.") so
      # they don't surface as fake species -- their stems still count in all-trees.
      mutate(species = ifelse(is.na(genus) | genus == "" |
                                    is.na(species) | species == "" |
                                    tolower(species) %in% c("spp.", "spp", "unknown"),
                              NA_character_, paste(genus, species))) %>%
      # Woodland species are measured at DRC, not DBH.
      # Convert DRC -> DBH (inches) via Chojnacky & Rogers 1999. Their
      # equations are fit for pinyon, Gambel oak, and Utah juniper; we extend
      # them by genus, and use the juniper equation for all other woodland genera.
      # Juniper is the fallback because the genera (mesquite, mountain-mahogany, etc.)
      # share its low-branching, multi-stemmed growth form mroe than pinyon or oak forms.
      # stm assumed multi-stem (0); floored at 0 where the linear model goes negative for
      # small stems. Used for BA only -- the sapling/tree split below stays on
      # the field-measured diameter.
      mutate(is_woodland = !is.na(woodland) & toupper(woodland) == "Y",
             dia_ba = case_when(
                   !is_woodland       ~ dia,
                   genus == "Pinus"   ~ pmax(-1.9374 + 0.9823 * dia, 0), # pinyons
                   genus == "Quercus" ~ pmax(-1.4599 + 0.9533 * dia, 0), # dryland/scrub oaks
                   TRUE               ~ pmax(-2.6843 + 1.0222 * dia, 0))) %>% # junipers + other woodland
      group_by(lon, lat, forest_type, species, samp_method_cd) %>%
      summarize(basal_area = sum(0.0054542 * dia_ba^2 * tpa_unadj, na.rm = T), # ft2/acre, all live stems (DBH basis; woodland DRC->DBH converted)
                n_saplings = sum(tpa_unadj[dia < 5], na.rm = T),  # stems/acre, saplings 1.0-4.9" DBH (microplot)
                n_trees = sum(tpa_unadj[dia >= 5], na.rm = T),    # stems/acre, trees >=5.0" DBH (subplot/macroplot)
                # seed_tpa is one value per subplot-species but the tree join
                # replicates it across each of that species' tree rows, so count
                # each subplot-species once.
                n_seedlings = sum(seed_tpa[!duplicated(paste(subp, spcd))], na.rm = T), # stems/acre, seedlings <1" DBH (microplot)
                # northness = mean(northness, na.rm = TRUE),
                # eastness = mean(eastness, na.rm = TRUE),
                .groups = "drop")

# note: if we keep topo-interpolation downstream, may want to retain subplot resolution at this stage



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
veg <- rast("/Volumes/T7/landfire/LF2020_BPS_220_CONUS/Tif/LC20_BPS_220.tif") # load raster
dp <- project(dp, crs(veg))
if(crs(veg) != crs(dp)) stop()
veg <- raster::extract(veg, dp)

# filter on veg
nonforest_veg <- c("Barren-Rock/Sand/Clay", "Grassland", "Perennial Ice/Snow", "Riparian",
                   "Savanna", "Shrubland", "Sparse")
nonforest <- foreign::read.dbf("/Volumes/T7/landfire/LF2020_BPS_220_CONUS/Tif/LC20_BPS_220.tif.vat.dbf") %>%
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
