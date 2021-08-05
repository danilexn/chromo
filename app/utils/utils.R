null.if.empty <- function(x) {
    if(length(x) == 0) {
        return(NULL)
    }

    return(x)
}

is.empty.or.except <- function(x) {
    try(return(is.empty(x)))
    return(TRUE)
}

sourceFolder <- function(folder, recursive = FALSE, ...)
{
    files <- list.files(folder, pattern = "[.][rR]$",
                        full.names = TRUE, recursive = recursive)
    if (!length(files))
        stop(simpleError(sprintf('No R files in folder "%s"', folder)))
    src <- invisible(lapply(files, source, ...))
    message(sprintf('%s files sourced from folder "%s"', length(src), folder))
}

p.adjust.inv <- function (p, method = p.adjust.methods, n = length(p)) 
{
    # Only implemented for holm and fdr
    method <- match.arg(method)
    nm <- names(p)
    p <- as.numeric(p)
    p0 <- setNames(p, nm)
    if (all(nna <- !is.na(p)))
        nna <- TRUE
    p <- p[nna]
    lp <- length(p)
    stopifnot(n >= lp)
    if (n <= 1)
        return(p0)
    p0[nna] <- switch(method,
    holm = {
        i <- seq_len(lp)
        o <- order(p)
        ro <- order(o)
        pmin(1, cummin(p[o] / (n + 1L - i)))[ro]
    }, fdr = {
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummax(p[o] / (n/i)))[ro]
    }, none = p)
    p0
}