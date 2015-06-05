#' plotCorrelation
#'
#' This function plot the raw correlation in a file: the file structure is - ID EXPERIMENTAL value  extimation value1   extimation value2   ..
#'
#' @param file the file to be analyzed
plotCorrelation <- function(file)
{
    library(ggplot2)
    data <- read.table(file, sep="\t",header = TRUE)
    dat <- data[-1]

    # Set color by cond
    ggplot(dat, aes(x=xvar, y=yvar, color=cond)) + geom_point(shape=1)

    # Same, but with different colors and add regression lines
    plot <- ggplot(dat, aes(x=xvar, y=yvar, color=cond)) +
        geom_point(shape=1) +
        scale_colour_hue(l=50) + # Use a slightly darker palette than normal
        geom_smooth(method=lm,   # Add linear regression lines
                    se=FALSE)    # Don't add shaded confidence region

    # Extend the regression lines beyond the domain of the data
    plot <- ggplot(dat, aes(x=xvar, y=yvar, color=cond)) + geom_point(shape=1) +
        scale_colour_hue(l=50) + # Use a slightly darker palette than normal
        geom_smooth(method=lm,   # Add linear regression lines
                    se=FALSE,    # Don't add shaded confidence region
                    fullrange=TRUE) # Extend regression lines


    # Set shape by cond
    plot <- ggplot(dat, aes(x=xvar, y=yvar, shape=cond)) + geom_point()

    # Same, but with different shapes
    plot <- ggplot(dat, aes(x=xvar, y=yvar, shape=cond)) + geom_point() +
        scale_shape_manual(values=c(1,2))  # Use a hollow circle and triangle


    return(plot)

}
