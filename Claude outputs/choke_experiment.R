# Connectivity measures for the default pinyon case, outside Shiny, for
# experimenting with the choke point measure.
#
# Run from the project root. Uses the app's own code: the model run
# (ts_model) and server/connectivity.R's habitat, land use, conductance and
# solve, so results match the app at the same settings.

# App globals: packages, R/, helpers/, core_data
local({
      lines <- readLines("app.R")
      eval(parse(text = lines[seq_len(grep("^# UI ----", lines) - 1)]), envir = globalenv())
})

# ---- Settings ----------------------------------------------------------------
settings <- ts_settings()     # default: Pinus edulis occurrence, g = 10, the default bbox
reach    <- 100               # dispersal reach, km
landuse  <- 100               # weight on land use, %
era      <- "hist"            # "hist" or "fut"

# ---- Model run ---------------------------------------------------------------
result <- ts_model(settings, core_data)

# server/connectivity.R, with its reactives as plain functions of these values
cs <- new.env()
cs$aim       <- function() result$surfaces
cs$model     <- list(cs_landuse = landuse, cs_reach = reach)
cs$core_data <- core_data
cs$reactive  <- function(expr) { expr <- substitute(expr); function() eval(expr, cs) }
sys.source("server/connectivity.R", envir = cs)

habitat      <- cs$cs_habitat()
permeability <- cs$cs_permeability()
conductance  <- cs$cs_conductance(habitat[[era]], permeability, landuse)

# ---- Solve: every layer circuitscaper returns --------------------------------
# The same call as cs_compute(), keeping all layers
t0 <- Sys.time()
r <- circuitscaper::dp_omni(
      conductance, reach,
      source_strength = habitat[[era]],
      radius_units = "km", resistance_is = "conductances",
      source_threshold = cs$cs_threshold, seed = cs$cs_seed, power = TRUE)
message("solved in ", round(difftime(Sys.time(), t0, units = "secs")), " s; layers: ",
        paste(names(r), collapse = ", "))

power      <- r$power
even_power <- r$power / r$normalized_power   # power on an even landscape
even_power[!is.finite(even_power)] <- NA
mean_even  <- global(even_power, "mean", na.rm = TRUE)$mean

# ---- Choke point candidates --------------------------------------------------
# app:       power / mean even-landscape power (the app's former measure)
# per_cell:  power / the cell's own even-landscape power (r$normalized_power)
# shrunk:    power / (cell's even power + lambda * mean even power), which is
#            per_cell where even power is high and app where it is low
choke <- function(lambda) power / (even_power + lambda * mean_even)

candidates <- rast(list(app = power / mean_even,
                per_cell = r$normalized_power,
                shrunk_0.5 = choke(0.5),
                shrunk_1   = choke(1),
                shrunk_2   = choke(2)))

# Check: the solve matches the app's (which maps normalized current)
app_layers <- cs$cs_compute(era, habitat, reach, permeability, landuse)
stopifnot(isTRUE(all.equal(values(app_layers$normalized_current), values(r$normalized_current),
                           check.attributes = FALSE)))

# ---- Look --------------------------------------------------------------------
# Same scale as the app's legend: log10, clamped to 0.01 .. 10
show <- function(x) clamp(log10(x), -2, 1, values = TRUE)
plot(show(candidates), range = c(-2, 1), col = rev(hcl.colors(100, "Spectral")),
     nc = 3, mar = c(1, 1, 2, 4))

# Other layers to compare against: habitat, conductance, corridors, raw current
plot(rast(list(habitat = habitat[[era]], log10_conductance = log10(conductance),
       corridors = log10(r$normalized_current), log10_current = log10(r$current))),
     nc = 2)
