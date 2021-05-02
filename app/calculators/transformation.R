rotate.axes.maximum <- function(df, coords, particles) {
    df <- df %>% group_by(!!sym(particles)) %>%
        mutate(
            !!sym(coords[1]) := homogeneize.rotation(!!sym(coords[2]),!!sym(coords[1]))[, 2],!!sym(coords[2]) := homogeneize.rotation(!!sym(coords[2]),!!sym(coords[1]))[, 1]
        )
}

homogeneize.rotation <- function(x, y) {
    coords <- cbind(X = x, Y = y)
    rad <- seq(0, pi, l = 20)
    best.rotation <- c()
    for (i in rad) {
        coords.rot <- Rotation(coords, i)
        best.rotation <-
            c(best.rotation, abs(max(coords.rot[, 1]) - min(coords.rot[, 1])))
    }

    rot <- rad[which.min(best.rotation)]
    coords.rot <- Rotation(coords, rot)
    return(coords.rot)
}