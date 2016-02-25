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
#' @param data dataset frame or table
#' @param na.rm ignore the null values
#'
plotHistFunc <- function(data, na.rm = TRUE, ...) {
    nm <- names(data)
    #    breaks <- seq(from = 0, to = 14, by = 0.5)
    #    print(breaks)
    plots <- list()  # new empty list
    for (i in seq_along(nm)) {
        plot <- ggplot(data, aes_string(x = nm[i])) + geom_histogram(alpha = .5,fill = "mediumseagreen") +
                scale_x_continuous(limits=c(2.5, 11)) +
                theme_bw() +
                ylab("Number of Observations") +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
        plots[[i]] <- plot  # add each plot into plot list
    }
    return (plots)
}


#' plotRawCorrelation
#'
#' This fucntion plots the raw correlation of all the predicted varaibles vs the expted variable (first column)
#' @param data frame
#'

plotRawCorrelation <- function(dat){
    nm <- names(dat)
    plots <- list()  # new empty list
    for (i in nm) {
        if(dat[1L] != dat[i]){
            newData <- data.frame(x = dat[1L], y = dat[i])
            colnames(newData) <- c("x", "y")
            #plot <- ggplot(x,aes_string(x = nm[i])) + geom_histogram(alpha = .5,fill = "mediumseagreen")
            plot <- ggplot(newData, aes(x=x, y=y)) +
                geom_point(shape=1, alpha = .5,size = 1, colour="mediumseagreen") +
                scale_colour_hue(l=50) + # Use a slightly darker palette than normal
                geom_smooth(method=lm, se=FALSE) + # Add shaded confidence region
                xlab(nm[1L]) +
                ylab(i) +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
            plots[[i]] <- plot  # add each plot into plot list
        }
    }
    return (plots)
}

#' plotBoxPlotCorrelation
#'
#' This fucntion plots the raw correlation of all the predicted varaibles vs the expted variable (first column)
#' @param data frame
#'

plotBoxPlotCorrelation <- function(dat){
    nm <- names(dat)
    plots <- list()  # new empty list
    for (i in nm) {
        if(dat[1L] != dat[i]){
            newData <- data.frame(x = dat[1L], y = dat[i])
            colnames(newData) <- c("x", "y")
            plot <- ggplot(newData, aes(x=x, y=y, group = x)) +
                geom_boxplot(aes(colour=factor(x)), fill=NA, outlier.shape = 1, outlier.size = 0.2, outlier.colour = "black", show_guide=FALSE) +
                xlab(nm[1L]) +
                ylab(i) +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
            plots[[i]] <- plot  # add each plot into plot list
        }

    }
    return (plots)
}

#' plotCorrelationMatrix
#'
#' This function plot the Correlation matrix between all the columns
#' @param data DataFRame
#'

plotCorrelationMatrix <- function(data, filter = NULL){
    datMy.scale<- scale(data); #scale all the features in the dataFRame
    corMatMy <- cor(datMy.scale) #compute the correlation matrix
    if(!is.null(filter)){
        highlyCor <- findCorrelation(corMatMy, filter) #Apply correlation filter
        datMyFiltered.scale <- datMy.scale[,-highlyCor] #then we remove all the variable correlated with more 0.7.
        corMatMy <- cor(datMyFiltered.scale)
        plot <- corrplot(corMatMy)
    }else{
        plot <- corrplot(corMatMy)  #visualize the matrix, clustering features by correlation index.
    }
    return (plot)
}

<<<<<<< HEAD
#' plotFractionCorrelation
#'
#' This fucntion plots the raw correlation of all the predicted varaibles (theoretical mean of the fraction) vs the expted variable (first column). It plots error bars and confidence intervals using a lineal model.
#' @param data frame
#'

plotFractionCorrelation <- function (dat){
nm <- names(dat)
plots <- list()  # new empty list
for (i in nm) {
        if(dat[1L] != dat[i]){
             newData <- data.frame(x = dat[1L], y = dat[i])
             colnames(newData) <- c("x", "y")
             fractions_stats <- ddply(newData, ~x, summarise, y_mean=mean(y), y_sd=sd(y)) # split data in fraction. compute mean and standard deviation
             # print(fractions_stats)
             plot <- ggplot(fractions_stats, aes(x=x, y=y_mean))+
             geom_errorbar(aes(ymin=y_mean-y_sd, ymax=y_mean+y_sd), width=.1) +
                           geom_point() +
                           geom_smooth(method="lm",formula=y~x) +
                           xlab (nm[1L]) +
                           ylab(i) +
                           theme_bw()+
                           theme(panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line = element_line(colour="black"))
             plots[[i]] <- plot # add each plot into plot list
        }
 
     }
    return (plots)
 }


=======
>>>>>>> ff7e566e64c270577ddf07143ee6d0751213acbe



