#' plotRawCorrelation
#'
#' This function plot the raw correlation between using the table variables
#' @param  file
#'
plotRawCorrelation <- function(file)
{
    data <- read.table(file, sep="\t",header = TRUE)
    dat  <- removeFirstColumn(data)
    dat  <- processData(dat)

    # Same, but with different colors and add regression lines
    plot <- ggplot(dat, aes(x=EXPERIMENTAL, y=SVM)) + geom_point(shape=1) +
        scale_colour_hue(l=50) + # Use a slightly darker palette than normal
        geom_smooth(method=lm,   # Add linear regression lines
                    se=FALSE)    # Don't add shaded confidence region
    return (plot)
}
