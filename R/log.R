# Usage log ----------------------------------------------------------------------
#
# What people do with the app: one JSON object per line (JSON Lines), one file
# per day per R process, so processes sharing the directory never interleave
# writes. Read a whole directory with scripts/usage_summary.R.
#
# Logging never gets in the user's way: a log that can't be opened or written is
# skipped with a warning.

#' Open the usage log
#'
#' Creates the directory if needed and deletes files older than `keep_days`.
#'
#' @param dir Log directory: `TREESHIFT_LOG_DIR` if set, else `logs/` in the
#'   working directory (the app directory, for the app).
#' @return A `ts_log`, for `ts_log_write()`.
ts_log_open <- function(dir = Sys.getenv("TREESHIFT_LOG_DIR", "logs"), keep_days = 365) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(dir)) stop("Can't create the log directory ", dir)
      ts_log_prune(dir, keep_days)
      structure(list(dir = normalizePath(dir)), class = "ts_log")
}

# Delete the log files dated more than keep_days ago (by the date that names
# them, not the file time, which copying can change)
ts_log_prune <- function(dir, keep_days) {
      files <- list.files(dir, pattern = "^\\d{4}-\\d{2}-\\d{2}-.*[.]jsonl$")
      dates <- as.Date(substr(files, 1, 10))
      old <- files[!is.na(dates) & dates < as.Date(format(Sys.time(), tz = "UTC")) - keep_days]
      unlink(file.path(dir, old))
}

#' Append an event to the usage log
#'
#' @param log `ts_log_open()`, or NULL for no logging.
#' @param event The event type.
#' @param fields A named list of the event's fields; NULL fields are written
#'   as null.
ts_log_write <- function(log, event, fields = list()) {
      if (is.null(log)) return(invisible())
      now <- Sys.time()
      line <- c(list(time = format(now, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"), event = event), fields)
      file <- file.path(log$dir, sprintf("%s-pid%d.jsonl", format(now, "%Y-%m-%d", tz = "UTC"),
                                         Sys.getpid()))
      tryCatch(cat(jsonlite::toJSON(line, auto_unbox = TRUE, null = "null", na = "null", digits = NA),
                   "\n", file = file, append = TRUE, sep = ""),
               error = function(e) warning("Couldn't write to the usage log: ", conditionMessage(e),
                                           call. = FALSE))
      invisible()
}

#' Read the usage log
#'
#' @param dir The log directory.
#' @return A tibble, one row per event, with the fields events share as
#'   columns. Fields an event doesn't have are NA.
ts_log_read <- function(dir = Sys.getenv("TREESHIFT_LOG_DIR", "logs")) {
      files <- sort(list.files(dir, pattern = "[.]jsonl$", full.names = TRUE))
      lines <- unlist(lapply(files, readLines, warn = FALSE))
      if (!length(lines)) return(tibble::tibble())
      dplyr::bind_rows(lapply(lines, function(l) {
            x <- jsonlite::fromJSON(l)
            # Null fields are NA; fields with several values (none so far)
            # are kept whole as text
            tibble::as_tibble(lapply(x, function(v) if (is.null(v)) NA else
                  if (length(v) == 1) v else as.character(jsonlite::toJSON(v))))
      })) |>
            dplyr::mutate(time = as.POSIXct(time, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"))
}
