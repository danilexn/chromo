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