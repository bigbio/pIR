#' plotRawCorrelation
#'
#' This function plot the raw correlation between using the table variables
#'
#' @param  file
#'

plotRawCorrelation <- function(file)
{
    data <- read.table(file, sep="\t",header = TRUE)
    dat  <- removeFirstColumn(data)
    dat  <- processData(dat)

    for(i in names(dat)){
        df$paste(i,'length',sep="_") <- str_length(df$i)
    }
    # Same, but with different colors and add regression lines
    plot <- ggplot(dat, aes(x=EXPERIMENTAL, y=SVM)) + geom_point(shape=1) +
        scale_colour_hue(l=50) + # Use a slightly darker palette than normal
        geom_smooth(method=lm)     # Add shaded confidence region
    return (plot)
}
