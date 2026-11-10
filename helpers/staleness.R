## STALENESS -------------------------------------------
##
## Several outputs in this app are computed on an explicit button rather than
## reactively, which means they can outlive the settings that produced them.
## The mechanism is the same in each case: snapshot the inputs a result depends
## on at compute time ("stamp"), compare against the current inputs on every
## render, and mark the output rather than destroying it.
##
## Marking, not clearing, is deliberate. These results can take minutes to
## produce and stay readable after the settings move; blanking them throws away
## something the user may still want. Fading plus an explicit warning makes the
## result unmistakably provisional while leaving it legible.

# Which stamped fields differ between a stored stamp and the current one.
# `exempt` names fields that were deliberately varied by the run itself and so
# must not count as staleness -- a bandwidth sweep whose whole purpose is
# comparing bandwidths would otherwise mark itself stale the moment the user
# acts on its recommendation.
stamp_diff <- function(now, old, exempt = character()) {
      if (is.null(old)) return(character())
      keys <- setdiff(intersect(names(now), names(old)), exempt)
      keys[!vapply(keys, function(k) identical(now[[k]], old[[k]]), logical(1))]
}

# Human-readable phrase for a set of changed fields.
stale_phrase <- function(fields, labels) {
      nm <- unname(labels[fields])
      nm <- nm[!is.na(nm)]
      if (!length(nm)) return("settings")
      if (length(nm) == 1) return(nm)
      paste(paste(utils::head(nm, -1), collapse = ", "), "and", utils::tail(nm, 1))
}

# Inline warning banner. Returns NULL when nothing has changed, so it can be
# dropped straight into a renderUI.
stale_banner <- function(fields, labels,
                         action = "Re-run to update these results.") {
      if (!length(fields)) return(NULL)
      div(class = "ts-notice",
          sprintf("These results were produced before the %s changed. %s",
                  stale_phrase(fields, labels), action))
}

# Fade a ggplot and caption it when stale. Implemented as an in-plot overlay
# rather than CSS opacity so it needs no JS round-trip and survives into a
# downloaded image. Added last, so it sits above the plot's own layers.
gg_stale <- function(p, stale, label = "Stale \u2014 re-run to update") {
      if (!isTRUE(stale)) return(p)
      p +
            annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                     fill = "white", alpha = 0.5) +
            labs(caption = label) +
            # Same color as the .ts-notice rule, so an in-plot caption and a
            # sidebar notice read as the same kind of statement.
            theme(plot.caption = element_text(color = "#9a3412", face = "bold",
                                              hjust = 0, size = 11))
}
