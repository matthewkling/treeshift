# The usage log: the writer and reader (R/log.R), and what a session records
# (server/usage_log.R)

withr::with_dir(app_dir, source("scripts/load_core.R"))

test_that("events round-trip through the log, one line each", {
      dir <- withr::local_tempdir()
      log <- ts_log_open(dir)
      ts_log_write(log, "a", list(query = "?e=sp&sp=Pinus_edulis", n = 3L, missing = NULL, na = NA))
      ts_log_write(log, "b", list(seconds = 1.25))
      files <- list.files(dir)
      expect_length(files, 1)
      expect_match(files, sprintf("^\\d{4}-\\d{2}-\\d{2}-pid%d[.]jsonl$", Sys.getpid()))
      expect_length(readLines(file.path(dir, files)), 2)

      x <- ts_log_read(dir)
      expect_identical(x$event, c("a", "b"))
      expect_identical(x$query[1], "?e=sp&sp=Pinus_edulis")
      expect_identical(x$n[1], 3L)
      expect_true(is.na(x$missing[1]) && is.na(x$na[1]))
      expect_identical(x$seconds[2], 1.25)
      expect_s3_class(x$time, "POSIXct")
      expect_lt(abs(as.numeric(difftime(x$time[1], Sys.time(), units = "secs"))), 60)
})

test_that("opening the log deletes files more than a year old", {
      dir <- withr::local_tempdir()
      today <- as.Date(format(Sys.time(), tz = "UTC"))
      old <- paste0(format(today - 366), "-pid1.jsonl")
      recent <- paste0(format(today - 364), "-pid1.jsonl")
      other <- "notes.txt"
      for (f in c(old, recent, other)) writeLines("{}", file.path(dir, f))
      ts_log_open(dir)
      expect_setequal(list.files(dir), c(recent, other))
})

test_that("a log that can't be written warns and carries on", {
      log <- structure(list(dir = file.path(withr::local_tempdir(), "gone")), class = "ts_log")
      expect_warning(ts_log_write(log, "a"), "Couldn't write to the usage log")
      expect_silent(ts_log_write(NULL, "a"))
})

test_that("a session logs its start, runs, views, downloads and end", {
      skip_on_cran()
      load_app_globals()
      dir <- withr::local_tempdir()
      withr::local_envvar(TREESHIFT_LOG_DIR = dir)

      app <- shinytest2::AppDriver$new(app_dir, width = 1400, height = 900,
                                       load_timeout = 600000, timeout = 600000)
      on.exit(app$stop(), add = TRUE)
      settle(app)
      base <- sub("\\?.*$", "", app$get_url())
      query <- "?e=sp&sp=Pinus_ponderosa&bx=-112,-106&by=37,42"

      s <- open_session(paste0(base, query))
      download_raster(s)
      s$stop()
      # The server logs the session's end once it notices the closed connection
      app$wait_for_idle(1000)
      Sys.sleep(2)

      x <- ts_log_read(dir)
      # As the browser sent it, commas encoded
      start <- dplyr::filter(x, event == "session_start", grepl("sp=Pinus_ponderosa", query))
      expect_equal(nrow(start), 1)
      expect_true(start$from_link)
      expect_true(nzchar(start$ip))
      expect_match(start$user_agent, "Chrome")
      ev <- dplyr::filter(x, session == start$session)

      run <- dplyr::filter(ev, event == "model_run")
      expect_equal(nrow(run), 1)
      expect_identical(run$kind, "surfaces")
      expect_identical(run$outcome, "ok")
      expect_identical(run$cache, "miss")
      expect_gt(run$seconds, 0)
      expect_gt(run$cells, 0)
      expect_match(run$query, "sp=Pinus_ponderosa", fixed = TRUE)

      expect_true(any(ev$event == "view" & grepl("sp=Pinus_ponderosa", ev$query, fixed = TRUE)))
      expect_identical(ev$what[ev$event == "download"], "raster")
      end <- dplyr::filter(ev, event == "session_end")
      expect_equal(nrow(end), 1)
      expect_equal(end$runs, 1)
      expect_identical(unique(ev$version), app_info$app)
})
