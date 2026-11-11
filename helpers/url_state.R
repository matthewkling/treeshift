# URL state: the link that reproduces a set of results.
#
# The URL's query string records every model setting, plus what the map shows,
# so a link reopens the same results (for the same app and data versions). It
# is a curated schema with short, stable keys and a `v=` version token, because
# links get pasted into reports and emails and have to keep working:
#
#   - Keys are never renamed or reused. A change of meaning is a new key, and
#     url_decode() keeps reading the old one.
#   - Model settings are written whenever they apply, even at their default, so
#     a link keeps its meaning if a default changes. A key that is missing from
#     a link takes the default, so short hand-written links work too.
#   - What the map shows (raster and point layers, analog direction) is written
#     only when it differs from the default.
#   - Settings that don't apply are left out: a species link has no forest
#     type, and connectivity settings appear only when connectivity is on.
#
# Every value is checked against the app's own choices and bounds. A value that
# fails falls back to the default and is reported, never silently dropped, and
# nothing in a link is ever evaluated or rendered as HTML.
#
# State is a named list of app values keyed by input id (track, species,
# bbox_x, ...), plus `site` (lon, lat or NULL) and `polygon_unit` (a library
# id or ""). url_defaults() is the app's starting state, so the UI and the
# codec share one set of defaults. `ctx` holds the valid choices; see
# url_context().
#
# Schema v1 keys, in the order links write them:
#   v      schema version
#   e      entity: sp, ft, fg, all
#   es     single or community: single, multi
#   sp     species (spaces as underscores);  ft  forest type code;  fg  group code
#   st     life stage: seedling, sapling, adult, any
#   m      variable: occ, abund, ba, rich, pba, regen
#   t      timeframe: hist, fut, chg
#   cm     community measure;  n  categories shown;  sf  support floor
#   g      grain (km);  bx, by  bounding box (lon range, lat range)
#   poly   focal polygon (boundary library id);  site  focal site (lon,lat)
#   cv     climate variables;  ssp, era, gcm  scenario
#   tg, tc geographic and climate bandwidths
#   cs     connectivity on (1); csrc, ctgt, ccl, cll, ccs, cls its settings
#   r      raster layer;  rcs  connectivity layer;  pts  point layer
#   an     analog direction

url_schema_version <- 1L

url_keys <- c("v", "e", "es", "sp", "ft", "fg", "st", "m", "t", "cm", "n", "sf",
              "g", "bx", "by", "poly", "site", "cv", "ssp", "era", "gcm", "tg", "tc",
              "cs", "csrc", "ctgt", "ccl", "cll", "ccs", "cls", "r", "rcs", "pts", "an")

# Short codes for enumerated settings: code = app value
url_codes <- list(
      track = c(sp = "Species", ft = "Forest type", fg = "Forest group", all = "All trees"),
      metric = c(occ = "occurrence", abund = "abundance", ba = "basal area",
                 rich = "species richness", pba = "proportion basal area",
                 regen = "regeneration fraction"),
      time = c(hist = "historical", fut = "future", chg = "delta"),
      comm_measure = c(rich = "richness", div = "diversity", even = "evenness",
                       dom = "dominance", domid = "dominant entity", turn = "turnover",
                       gain = "expected gains", loss = "expected losses",
                       coldem = "colonization demand", lossfrac = "loss fraction",
                       win = "biggest winner", lose = "biggest loser"),
      electrode = c(cells = "grid cells", plots = "FIA plots", site = "selected site"),
      conductance = c(suit = "suitability", wind = "wind", suitwind = "suitability * wind",
                      uniform = "uniform"),
      loss = c(none = "none", uniform = "uniform", suit = "mean suitability"),
      stat = c(fia = "FIA variable", ess = "ESS", sim = "focal site similarity",
               conn = "connectivity"),
      cs_stat = stats::setNames(cs_stat_names, cs_stats),
      vector = c(occ = "occupied FIA plots", all = "all FIA plots", analogs = "focal site analogs",
                 anim = "current animation", field = "current field", none = "none"),
      analog = c(contemp = "Contemporary analogs", reverse = "Reverse analogs",
                 forward = "Forward analogs")
)

# Slider bounds, shared with the UI so a link can't hold a value the UI can't
url_bounds <- list(
      clim_agg = c(1, 10), theta_geog = c(5, 100), theta_clim = c(0.05, 0.5),
      comm_top_n = c(2, 25), cs_conductance_scale = c(-2, 2), cs_loss_scale = c(-3, 2)
)

# The app's valid choices, gathered once at startup (app.R). Choice vectors
# are named by their menu labels, which the settings record uses.
url_context <- function() {
      as_named <- function(x) stats::setNames(as.character(unlist(x)), names(x))
      list(species = spp,
           forest_types = as_named(forest_type_choices),
           forest_groups = as_named(forest_group_choices),
           ssps = as_named(ssp_choices), periods = as_named(period_choices),
           models = as_named(model_choices),
           default_ssp = default_ssp, default_period = default_period,
           default_model = default_model,
           units = boundary_units$id,
           unit_names = stats::setNames(boundary_units$name, boundary_units$id),
           bbox = clim_bbox)
}

url_defaults <- function(ctx) {
      list(track = "Species", entity_scope = "single", species = "Pinus edulis",
           forest_type = unname(ctx$forest_types[1]), forest_group = unname(ctx$forest_groups[1]),
           stage = "any", metric = "occurrence", fia_var = "presence probability",
           time = "historical", comm_measure = "richness", comm_top_n = 10,
           # Its default depends on the variable (support_floor_spec())
           support_floor = support_floor_spec(FALSE, "presence probability", "Species")$default,
           clim_agg = 10, bbox_x = c(-116, -100), bbox_y = c(29, 43),
           polygon_unit = "", site = NULL,
           clim_vars = url_menu_order(clim_vars_default, clim_vars()), scenario_ssp = ctx$default_ssp,
           scenario_period = ctx$default_period, scenario_model = ctx$default_model,
           theta_geog = 30, theta_clim = 0.1,
           cs_enable = FALSE, cs_sources = "grid cells", cs_targets = "grid cells",
           cs_conductance_layer = "suitability", cs_loss_layer = "none",
           cs_conductance_scale = 0, cs_loss_scale = -1,
           stat = "FIA variable", cs_stat = cs_stat_names[1],
           vector = "occupied FIA plots", analog_direction = "Contemporary analogs")
}

# Which settings apply, from the settings they depend on. Shared by the encoder
# and decoder, and mirroring the UI's conditionalPanels.
url_multi <- function(s) identical(s$entity_scope, "multi") && !identical(s$track, "All trees")
url_cs_on <- function(s) isTRUE(s$cs_enable) && !url_multi(s)
url_has_stage <- function(s) s$track %in% c("Species", "All trees")
url_has_floor <- function(s) {
      (url_multi(s) && !s$comm_measure %in% c("richness", "expected gains", "expected losses")) ||
            identical(s$metric, "regeneration fraction")
}

# Measure choices as the UI offers them: valid ones only, else the default
url_comm_measure <- function(s, want = s$comm_measure) {
      ch <- comm_measure_choices(s$track, s$fia_var, s$time)
      if (length(ch) == 0 || want %in% ch) want else unname(ch[1])
}


# Encoding -------------------------------------------------------------------

# Numbers to at most four decimals, without trailing zeros
url_num <- function(x) {
      paste(sub("\\.?0+$", "", formatC(round(x, 4), format = "f", digits = 4)), collapse = ",")
}

# Percent-encode, keeping the separators the schema uses readable
url_escape <- function(x) {
      x <- utils::URLencode(x, reserved = TRUE)
      gsub("%2C", ",", gsub("%3A", ":", x, fixed = TRUE), fixed = TRUE)
}

code_of <- function(map, value) names(map)[match(value, map)]

url_menu_order <- function(x, menu) x[order(match(x, menu))]

# State -> named character vector of keys, in schema order
url_fields <- function(s, ctx) {
      d <- url_defaults(ctx)
      f <- c(v = as.character(url_schema_version), e = code_of(url_codes$track, s$track))
      if (!identical(s$track, "All trees")) f["es"] <- s$entity_scope
      if (!url_multi(s)) {
            if (identical(s$track, "Species"))      f["sp"] <- gsub(" ", "_", s$species)
            if (identical(s$track, "Forest type"))  f["ft"] <- s$forest_type
            if (identical(s$track, "Forest group")) f["fg"] <- s$forest_group
      }
      if (url_has_stage(s)) f["st"] <- s$stage
      f["m"] <- code_of(url_codes$metric, s$metric)
      f["t"] <- code_of(url_codes$time, s$time)
      if (url_multi(s)) {
            f["cm"] <- code_of(url_codes$comm_measure, s$comm_measure)
            if (s$comm_measure %in% categorical_measures && isTRUE(is.finite(s$comm_top_n)))
                  f["n"] <- url_num(s$comm_top_n)
      }
      # A cleared numeric field reads as NA; it is left out, so the link takes
      # the default rather than carrying an invalid value
      if (url_has_floor(s) && isTRUE(is.finite(s$support_floor))) f["sf"] <- url_num(s$support_floor)
      f["g"] <- url_num(s$clim_agg)
      f["bx"] <- url_num(s$bbox_x)
      f["by"] <- url_num(s$bbox_y)
      if (nzchar(s$polygon_unit)) f["poly"] <- s$polygon_unit
      if (!is.null(s$site)) f["site"] <- url_num(s$site)
      # Multi-selects are written in menu order, so the same selection always
      # gives the same link. They can be emptied in the UI; an empty one is
      # left out.
      if (length(s$clim_vars)) f["cv"] <- paste(url_menu_order(s$clim_vars, clim_vars()), collapse = ",")
      f["ssp"] <- s$scenario_ssp
      f["era"] <- s$scenario_period
      if (length(s$scenario_model))
            f["gcm"] <- paste(url_menu_order(s$scenario_model, ctx$models), collapse = ",")
      f["tg"] <- url_num(s$theta_geog)
      f["tc"] <- url_num(s$theta_clim)
      if (url_cs_on(s)) {
            f["cs"] <- "1"
            f["csrc"] <- code_of(url_codes$electrode, s$cs_sources)
            f["ctgt"] <- code_of(url_codes$electrode, s$cs_targets)
            f["ccl"] <- code_of(url_codes$conductance, s$cs_conductance_layer)
            f["cll"] <- code_of(url_codes$loss, s$cs_loss_layer)
            f["ccs"] <- url_num(s$cs_conductance_scale)
            f["cls"] <- url_num(s$cs_loss_scale)
      }
      if (!identical(s$stat, d$stat)) f["r"] <- code_of(url_codes$stat, s$stat)
      if (identical(s$stat, "connectivity")) f["rcs"] <- code_of(url_codes$cs_stat, s$cs_stat)
      if (!identical(s$vector, d$vector)) f["pts"] <- code_of(url_codes$vector, s$vector)
      if (!identical(s$analog_direction, d$analog_direction))
            f["an"] <- code_of(url_codes$analog, s$analog_direction)
      f[!is.na(f)]
}

# State -> query string, with the leading "?"
url_encode <- function(s, ctx) {
      f <- url_fields(s, ctx)
      paste0("?", paste(names(f), url_escape(f), sep = "=", collapse = "&"))
}


# Decoding -------------------------------------------------------------------

# Query string -> list(state, problems). `problems` are plain-text notes on
# link values that couldn't be used, for showing to the user.
url_decode <- function(search, ctx) {
      s <- url_defaults(ctx)
      q <- shiny::parseQueryString(if (is.null(search)) "" else search)
      problems <- character(0)
      if (length(q) == 0) return(list(state = s, problems = problems))

      note <- function(label, raw, why = "isn't a recognized value") {
            raw <- if (nchar(raw) > 40) paste0(substr(raw, 1, 40), "...") else raw
            problems <<- c(problems, sprintf("%s: '%s' %s, so the default was used.", label, raw, why))
      }
      # A coded value, or NULL (missing, or reported as invalid)
      code <- function(key, map, label, allowed = map) {
            raw <- q[[key]]
            if (is.null(raw)) return(NULL)
            v <- unname(map[raw])
            if (is.na(v)) { note(label, raw); return(NULL) }
            if (!v %in% allowed) { note(label, raw, "isn't available with the link's other settings"); return(NULL) }
            v
      }
      # One of a set of plain values, or NULL
      one_of <- function(key, allowed, label, transform = identity) {
            raw <- q[[key]]
            if (is.null(raw)) return(NULL)
            v <- transform(raw)
            if (!v %in% allowed) { note(label, raw); return(NULL) }
            v
      }
      # n comma-separated numbers within [lo, hi], or NULL
      nums <- function(key, label, lo, hi, n = 1, integer = FALSE) {
            raw <- q[[key]]
            if (is.null(raw)) return(NULL)
            x <- suppressWarnings(as.numeric(strsplit(raw, ",", fixed = TRUE)[[1]]))
            ok <- length(x) == n && !anyNA(x) && all(x >= lo & x <= hi) &&
                  (!integer || all(x == round(x)))
            if (!ok) {
                  what <- if (n == 1) "a number" else paste(n, "numbers")
                  if (integer) what <- if (n == 1) "a whole number" else paste(n, "whole numbers")
                  note(label, raw, if (is.finite(hi)) sprintf("isn't %s between %s and %s", what, url_num(lo), url_num(hi))
                       else sprintf("isn't %s of at least %s", what, url_num(lo)))
                  return(NULL)
            }
            x
      }
      # A nonempty comma-separated subset of `allowed`, or NULL
      subset_of <- function(key, allowed, label) {
            raw <- q[[key]]
            if (is.null(raw)) return(NULL)
            x <- strsplit(raw, ",", fixed = TRUE)[[1]]
            if (length(x) == 0 || anyDuplicated(x) || !all(x %in% allowed)) { note(label, raw); return(NULL) }
            url_menu_order(x, allowed)
      }
      set <- function(field, value) if (!is.null(value)) s[[field]] <<- value

      if (!is.null(q$v)) {
            v <- suppressWarnings(as.integer(q$v))
            if (is.na(v)) {
                  problems <- c(problems, "The link's version number isn't valid; its settings were read as the current version.")
            } else if (v > url_schema_version) {
                  problems <- c(problems, "This link was made by a newer version of TreeShift, so some of its settings may not have been applied.")
            }
      }
      unknown <- setdiff(names(q), url_keys)
      if (length(unknown)) {
            problems <- c(problems, sprintf("Ignored settings this version doesn't recognize: %s.",
                                            paste(utils::head(unknown, 5), collapse = ", ")))
      }

      # In dependency order: what an entry may be depends on the entries above it
      set("track", code("e", url_codes$track, "Entity"))
      if (!identical(s$track, "All trees"))
            set("entity_scope", one_of("es", c("single", "multi"), "Single or community"))
      if (!url_multi(s)) {
            if (identical(s$track, "Species"))
                  set("species", one_of("sp", ctx$species, "Species", function(x) gsub("_", " ", x)))
            if (identical(s$track, "Forest type"))
                  set("forest_type", one_of("ft", ctx$forest_types, "Forest type"))
            if (identical(s$track, "Forest group"))
                  set("forest_group", one_of("fg", ctx$forest_groups, "Forest group"))
      }
      if (url_has_stage(s)) {
            set("stage", one_of("st", c("seedling", "sapling", "adult", "any"), "Life stage"))
      } else {
            s$stage <- "any"
      }
      mc <- metric_choices(s$track, s$stage)
      if (!s$metric %in% mc) s$metric <- mc[1]
      set("metric", code("m", url_codes$metric, "Variable", allowed = mc))
      s$fia_var <- resolve_var(s$track, s$metric, s$stage)
      set("time", code("t", url_codes$time, "Timeframe"))

      if (url_multi(s)) {
            s$comm_measure <- url_comm_measure(s, s$comm_measure)
            set("comm_measure", code("cm", url_codes$comm_measure, "Community measure",
                                     allowed = comm_measure_choices(s$track, s$fia_var, s$time)))
            if (s$comm_measure %in% categorical_measures)
                  set("comm_top_n", nums("n", "Categories shown", url_bounds$comm_top_n[1],
                                         url_bounds$comm_top_n[2], integer = TRUE))
      }
      s$support_floor <- support_floor_spec(url_multi(s), s$fia_var, s$track)$default
      if (url_has_floor(s)) set("support_floor", nums("sf", "Support floor", 0, Inf))

      b <- ctx$bbox
      set("clim_agg", nums("g", "Grain size", url_bounds$clim_agg[1], url_bounds$clim_agg[2], integer = TRUE))
      bx <- nums("bx", "Longitude range", b[["xmin"]], b[["xmax"]], n = 2)
      if (!is.null(bx) && bx[1] >= bx[2]) { note("Longitude range", q$bx, "isn't west to east"); bx <- NULL }
      set("bbox_x", bx)
      by <- nums("by", "Latitude range", b[["ymin"]], b[["ymax"]], n = 2)
      if (!is.null(by) && by[1] >= by[2]) { note("Latitude range", q$by, "isn't south to north"); by <- NULL }
      set("bbox_y", by)
      if (!is.null(q$poly)) {
            if (q$poly %in% ctx$units) {
                  s$polygon_unit <- q$poly
            } else {
                  problems <- c(problems, sprintf(
                        "Focal polygon '%s' isn't in this version's boundary library, so no polygon was set.",
                        substr(q$poly, 1, 40)))
            }
      }
      site <- q$site
      if (!is.null(site)) {
            lon <- suppressWarnings(as.numeric(strsplit(site, ",", fixed = TRUE)[[1]]))
            if (length(lon) == 2 && !anyNA(lon) && lon[1] >= b[["xmin"]] && lon[1] <= b[["xmax"]] &&
                lon[2] >= b[["ymin"]] && lon[2] <= b[["ymax"]]) {
                  s$site <- lon
            } else {
                  problems <- c(problems, sprintf(
                        "Focal site '%s' isn't a longitude,latitude within the map area, so no site was set.",
                        substr(site, 1, 40)))
            }
      }

      set("clim_vars", subset_of("cv", clim_vars(), "Climate variables"))
      set("scenario_ssp", one_of("ssp", ctx$ssps, "Emissions pathway"))
      set("scenario_period", one_of("era", ctx$periods, "Future time period"))
      set("scenario_model", subset_of("gcm", ctx$models, "Climate models"))
      set("theta_geog", nums("tg", "Geographic bandwidth", url_bounds$theta_geog[1], url_bounds$theta_geog[2]))
      set("theta_clim", nums("tc", "Climate bandwidth", url_bounds$theta_clim[1], url_bounds$theta_clim[2]))

      if (!is.null(q$cs)) {
            if (!identical(q$cs, "1")) {
                  note("Connectivity", q$cs)
            } else if (url_multi(s)) {
                  problems <- c(problems, "Connectivity doesn't apply to community composition, so it was left off.")
            } else {
                  s$cs_enable <- TRUE
            }
      }
      if (url_cs_on(s)) {
            set("cs_sources", code("csrc", url_codes$electrode, "Connectivity sources"))
            set("cs_targets", code("ctgt", url_codes$electrode, "Connectivity targets"))
            set("cs_conductance_layer", code("ccl", url_codes$conductance, "Horizontal conductance"))
            set("cs_loss_layer", code("cll", url_codes$loss, "Connectivity loss"))
            set("cs_conductance_scale", nums("ccs", "Conductance scale", url_bounds$cs_conductance_scale[1],
                                             url_bounds$cs_conductance_scale[2]))
            set("cs_loss_scale", nums("cls", "Loss scale", url_bounds$cs_loss_scale[1], url_bounds$cs_loss_scale[2]))
      }

      set("stat", code("r", url_codes$stat, "Raster layer",
                       allowed = if (url_cs_on(s)) stat_choices_cs else stat_choices_base))
      if (identical(s$stat, "connectivity")) {
            set("cs_stat", code("rcs", url_codes$cs_stat, "Connectivity layer"))
            # The app shows connectivity as change only
            s$time <- "delta"
      }
      set("vector", code("pts", url_codes$vector, "Point layer",
                         allowed = if (url_cs_on(s)) vector_choices_cs else vector_choices_base))
      set("analog_direction", code("an", url_codes$analog, "Analog direction"))

      list(state = s, problems = problems)
}
