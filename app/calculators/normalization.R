normalize.01 <- function(x, na.rm = FALSE) {
    return((x - min(x)) / (max(x) - min(x)))
}

normalize.z <- function(x, na.rm = FALSE) {
    return(scale(x))
}

normalize.min <- function(x, na.rm = FALSE) {
    return((x) / (min(x)))
}

normalize.median <- function(x, na.rm = FALSE) {
    return(x - median(x))
}

normalize.max <- function(x, na.rm = FALSE) {
    return((x) / (max(x)))
}