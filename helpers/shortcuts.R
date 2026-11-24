# Shortcuts: common questions, each a set of settings that answers it. A
# shortcut changes only what its question needs and keeps the rest -- the
# entity, area, focal site, climate data and bandwidths -- so it doesn't throw
# away the user's setup.
#
# Each shortcut is a function of the current state (the URL codec's state,
# helpers/url_state.R) returning the new state, or a string saying why it
# doesn't apply. Its label is a function of the entity type, so questions name
# what the user is looking at ("this forest type").

# The settings a shortcut can change without the model being rebuilt; any
# other change reopens the page at the new link (server/shortcuts.R)
shortcut_display_ids <- c("stat", "evidence_stat", "vector", "analog_direction")

# The entity type in a question: "this species", "forest types"
shortcut_entity <- function(track, plural = FALSE) {
      # All trees has no entity of its own, so its (unavailable) questions
      # are worded for species, as its single-or-community menu is
      if (identical(track, "All trees")) track <- "Species"
      if (plural) track_noun(track, plural = TRUE)
      else paste("this", track_noun(track))
}

# One focal entity, on the FIA variable map. Doesn't apply to All trees, which
# has no entity to focus on. Questions are about all life stages unless they
# say otherwise, so a stage left over from an earlier question doesn't carry
# into the next.
shortcut_single <- function(s, metric, time, stage = "any") {
      if (identical(s$track, "All trees"))
            return("Choose a species, forest type or forest group under Entity first.")
      s$entity_scope <- "single"
      s$stage <- stage
      s$metric <- metric
      s$time <- time
      s$stat <- "FIA variable"
      s
}

# Every entity at once, on a community measure, for all life stages. Measures
# of change by default.
shortcut_community <- function(s, measure, metric = "occurrence", time = "delta") {
      if (identical(s$track, "All trees"))
            return("Choose species, forest types or forest groups under Entity first.")
      s$entity_scope <- "multi"
      s$stage <- "any"
      s$metric <- metric
      s$time <- time
      s$comm_measure <- measure
      s$stat <- "FIA variable"
      s
}

shortcuts <- list(
      abundance = list(
            label = function(track) sprintf("Where is %s most abundant today?", shortcut_entity(track)),
            apply = function(s) {
                  # Forest types and groups are measured by occurrence only
                  if (s$track %in% c("Forest type", "Forest group"))
                        return(shortcut_single(s, "occurrence", "historical"))
                  shortcut_single(s, "basal area", "historical")
            }),
      habitat_change = list(
            label = function(track) sprintf("How might habitat for %s change?", shortcut_entity(track)),
            apply = function(s) shortcut_single(s, "occurrence", "delta")),
      regeneration = list(
            label = function(track) "Are young trees growing where adults do?",
            apply = function(s) {
                  if (!identical(s$track, "Species"))
                        return("Regeneration is modeled for species only.")
                  # Seedlings, unless saplings are already chosen
                  stage <- if (identical(s$stage, "sapling")) "sapling" else "seedling"
                  shortcut_single(s, "regeneration fraction", "historical", stage)
            }),
      dominant = list(
            label = function(track) sprintf("Which %s dominate here today?",
                                            shortcut_entity(track, plural = TRUE)),
            apply = function(s) {
                  # By basal area, as dominance usually means; forest types and
                  # groups are measured by occurrence only
                  metric <- if (s$track %in% c("Forest type", "Forest group")) "occurrence" else "basal area"
                  shortcut_community(s, "dominant entity", metric, "historical")
            }),
      biggest_loser = list(
            label = function(track) sprintf("Which %s stand to lose the most here?",
                                            shortcut_entity(track, plural = TRUE)),
            apply = function(s) shortcut_community(s, "biggest loser")),
      turnover = list(
            label = function(track) "How different might the future forest be?",
            apply = function(s) shortcut_community(s, "turnover")),
      connectivity = list(
            label = function(track) sprintf("What corridors link current and future habitat for %s?",
                                            shortcut_entity(track)),
            apply = function(s) {
                  s <- shortcut_single(s, "occurrence", "delta")
                  if (is.character(s)) return(s)
                  s$stat <- "connectivity"
                  s
            }),
      evidence = list(
            label = function(track) "How much evidence supports this map?",
            apply = function(s) {
                  s$stat <- "evidence"
                  s$evidence_stat <- "support"
                  s
            })
)

# A shortcut's new state, or why it doesn't apply. A new variable or community
# measure starts at its own support floor: floors are in the units of the
# layer they mask, so the old one wouldn't carry over.
shortcut_state <- function(id, s) {
      new <- shortcuts[[id]]$apply(s)
      if (is.character(new)) return(new)
      new$fia_var <- resolve_var(new$track, new$metric, new$stage)
      layer <- c("entity_scope", "fia_var", "comm_measure")
      if (!identical(new[layer], s[layer]))
            new$support_floor <- support_floor_spec(url_multi(new), new$fia_var, new$track)$default
      new
}

# The questions as buttons, for state s: one that doesn't apply is disabled,
# with the reason under it
shortcut_list <- function(s) {
      tags$ul(class = "ts-shortcuts",
              lapply(names(shortcuts), function(id) {
                    res <- shortcut_state(id, s)
                    unavailable <- is.character(res)
                    tags$li(
                          tags$button(type = "button", class = "ts-shortcut",
                                      disabled = if (unavailable) NA,
                                      onclick = sprintf("Shiny.setInputValue('shortcut', '%s', {priority: 'event'})", id),
                                      shortcuts[[id]]$label(s$track)),
                          if (unavailable) div(class = "ts-shortcut-note", res))
              }))
}
