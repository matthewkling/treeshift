# Summarize the usage log (R/log.R; events in server/usage_log.R).
#
#   Rscript scripts/usage_summary.R [days]
#
# Run from the project root, as a user who can read the app's log directory (or
# with TREESHIFT_LOG_DIR set to it). Covers the last `days` days (default 30).

source("scripts/load_core.R")

args <- commandArgs(trailingOnly = TRUE)
days <- if (length(args)) as.numeric(args[1]) else 30

log <- ts_log_read()
if (!nrow(log)) stop("No usage log entries found")
log <- dplyr::filter(log, time >= Sys.time() - days * 86400)
if (!nrow(log)) stop("No usage log entries in the last ", days, " days")

show <- function(title, x) {
      cat("\n", title, "\n", sep = "")
      print(as.data.frame(x), row.names = FALSE)
}
# Fields only some events have are absent from a log without those events
fields <- c("from_link", "user_agent", "kind", "outcome", "cache", "seconds",
            "message", "what", "query")
for (f in setdiff(fields, names(log))) log[[f]] <- NA

cat(sprintf("Usage, last %g days: %s to %s\n", days,
            format(min(log$time), "%Y-%m-%d"), format(max(log$time), "%Y-%m-%d")))

starts <- dplyr::filter(log, event == "session_start")
show("Sessions", dplyr::summarize(starts, sessions = dplyr::n(),
                                  opened_from_link = sum(from_link %in% TRUE)))

runs <- dplyr::filter(log, event == "model_run")
if (nrow(runs)) {
      show("Model runs", runs |>
                 dplyr::group_by(kind, outcome, cache) |>
                 dplyr::summarize(runs = dplyr::n(),
                                  median_seconds = stats::median(seconds, na.rm = TRUE),
                                  max_seconds = suppressWarnings(max(seconds, na.rm = TRUE)),
                                  .groups = "drop"))
      # The most requested settings: candidates for precomputing. Links
      # include display settings, so runs differing only in those count apart.
      show("Most run settings", runs |>
                 dplyr::filter(outcome == "ok") |>
                 dplyr::count(kind, query, sort = TRUE) |>
                 utils::head(10))
      problems <- dplyr::filter(runs, outcome != "ok")
      if (nrow(problems))
            show("Refused or failed runs", dplyr::count(problems, kind, outcome, message, sort = TRUE))
}

downloads <- dplyr::filter(log, event == "download")
if (nrow(downloads)) show("Downloads", dplyr::count(downloads, what, sort = TRUE))

# Phone browsers mark themselves "Mobi"; tablets mostly don't (iPads report
# themselves as desktop Safari), so they count as desktop or tablet
show("Sessions by device", starts |>
           dplyr::mutate(device = dplyr::case_when(is.na(user_agent) ~ "unknown",
                                                   grepl("Mobi", user_agent) ~ "phone",
                                                   TRUE ~ "desktop or tablet")) |>
           dplyr::count(device, sort = TRUE))

# Visitors' IP addresses aren't in this log; the web server's access log has them
