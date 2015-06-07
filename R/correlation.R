#' plotCorrelation
#'
#' This function plot the raw correlation in a file: the file structure is:
#   ID Object   EXPERIMENTAL-value  Predicted Value1    Predicted Value2    Predicted Value3   ..
#'
#' @param file the file to be analyzed
#'
rawCorrelation <- function(file)
{
    data <- read.table(file, sep="\t",header = TRUE)
    dat  <- removeFirstColumn(data)
    data <- processData(data)

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
