# Thin wrappers over package functions, used throughout the app.

aggregate <- function(x, fact, ...){
      if(fact == 1) return(x)
      terra::aggregate(x, fact = fact, ...)
}

tooltip <- function(x, ...) bslib::tooltip(x, ..., options = list(delay = list(show = 2000, hide = 100)))
