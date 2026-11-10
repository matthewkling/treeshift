# Julia / Circuitscape runtime for the connectivity module.

# Lazy Julia/Circuitscape initialization.
#
# cs_setup() warms up the Julia runtime and JIT-compiles Circuitscape, which
# is slow (tens of seconds). Doing it here at global scope makes EVERY new R
# process pay that cost at startup before serving its first page, even those who
# never open the connectivity panel.
#
# Instead, defer cs_setup() until the first time connectivity is actually
# computed (see connectivity()), and run it at most once per R process via
# this guard. Most users never touch connectivity, so they never pay the cost;
# the first connectivity user in a given process pays it once.
#
# LONGER-TERM: the better fix is to eliminate the warmup cost rather than
# defer it, by building a custom Julia sysimage with PackageCompiler.jl that
# bakes Circuitscape's compiled code in, and pointing JuliaCall at it. That
# makes the per-process warmup near-instant, at which point eager global init
# would be fine again. Lazy init is the stopgap until then.
.cs_initialized <- FALSE
ensure_cs_setup <- function() {
      if (.cs_initialized) return(invisible(FALSE))
      cs_setup()
      .cs_initialized <<- TRUE
      invisible(TRUE)
}
