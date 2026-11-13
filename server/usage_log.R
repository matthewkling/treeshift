# USAGE LOG -----------------------------------------------
#
# What this session does, in the usage log (R/log.R). Every event carries the
# session's id and the app version. The events:
#   session_start  the link the page opened on, the page that linked to it,
#                  and the browser
#   model_run      a model result: where it came from and how long it took
#   view           the settings link, each time it changes (url_state.R)
#   download       an export (downloads.R)
#   session_end    how long the session lasted, and how many runs it made
# Settings are recorded as the link's query string (helpers/url_state.R), which
# reopens them in the app.

# A random id: links this session's events without identifying anyone
log_session_id <- substr(rlang::hash(list(session$token, Sys.time())), 1, 12)
log_started <- Sys.time()
log_runs <- 0L

log_event <- function(event, ...) {
      ts_log_write(usage_log, event, c(list(session = log_session_id, version = app_info$app),
                                       list(...)))
}

# The current settings link, or NULL before the settings are complete
log_query <- function() {
      tryCatch(isolate(url_encode(current_state(), url_ctx)), error = function(e) NULL)
}

# Logged when the browser reports itself (ui/head.R), just after it connects:
# Shiny Server passes no request headers on to the app, so the browser can't
# be read from the request. The first view is usually logged a moment before.
# The referrer is null when there was none (a typed or bookmarked address);
# for another site's link, browsers usually send only the site, not the page. Visitors'
# IP addresses aren't logged here; the web server's access log has them.
observeEvent(input$ts_client, once = TRUE, {
      query <- isolate(session$clientData$url_search)
      referrer <- input$ts_client$referrer
      log_event("session_start", query = query, from_link = nzchar(query),
                referrer = if (nzchar(referrer %||% "")) referrer,
                user_agent = input$ts_client$user_agent)
})

session$onSessionEnded(function() {
      log_event("session_end",
                minutes = round(as.numeric(difftime(Sys.time(), log_started, units = "mins")), 2),
                runs = log_runs)
})

#' A model result through the result cache, logged as a model_run
#'
#' @param kind "surfaces", "community" or "connectivity".
#' @param key,compute As for `ts_cached()`.
#' @param describe Function of the result, returning extra fields to log
#'   (e.g. its cell count).
logged_run <- function(kind, key, compute, describe = function(value) list()) {
      start <- Sys.time()
      source <- NULL
      log_run <- function(outcome, ...) {
            log_runs <<- log_runs + 1L
            log_event("model_run", kind = kind, outcome = outcome, cache = source,
                      seconds = round(as.numeric(difftime(Sys.time(), start, units = "secs")), 2),
                      query = log_query(), ...)
      }
      value <- withCallingHandlers(
            ts_cached(result_cache, key, compute, report = function(s) source <<- s),
            # req() failures are Shiny's normal "not ready yet", not errors
            error = function(e) if (!inherits(e, "shiny.silent.error"))
                  log_run("error", message = conditionMessage(e)))
      do.call(log_run, c(list("ok"), describe(value)))
      value
}
