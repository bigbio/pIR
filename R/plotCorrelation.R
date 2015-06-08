#' plotRawCorrelation
#'
#' This function plot the raw correlation between using the table variables
#'
#' @param  file the file that contains the expeted value and predicted values
#'

plotHistFrecuencyValues <- function(file)
{
    data <- read.table(file, sep="\t",header = TRUE)
    dat  <- removeFirstColumn(data)
    dat  <- processData(dat)

    plots <- plotHistFunc(dat) ## execute function

    return (plots)
}

#' multiplot
#'
#' This function allows to plot different ggplot charts in the same grids
#'
#' @param plotlists this param allow iterate thorugth a set of plots to
#' @param file if the user wants to plot to a file, it can use the file paramterer
#' @param cols Number of columns to be use.
#' @param layout the layout to be use

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {

    plots <- c(list(...), plotlist)
    numPlots = length(plots)

    if (is.null(layout)) {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots == 1) {
        print(plots[[1]])

    } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        for (i in 1:numPlots) {
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

#' plotHistFunc
#'
#' This function plots the histogram for each value of the datasets, it returns a set of plots
#' @param x dataset frame or table
#' @param na.rm ignore the null values
#'
plotHistFunc <- function(x, na.rm = TRUE, ...) {
    nm <- names(x)
    plots <- list()  # new empty list
    for (i in seq_along(nm)) {
        plot <- ggplot(x,aes_string(x = nm[i])) + geom_histogram(alpha = .5,fill = "mediumseagreen")
        plots[[i]] <- plot  # add each plot into plot list
    }
    return (plots)
}
