
data <- read.table("/Users/yperez/IdeaProjects/predictedAnalysis/data/Peptide-PI-Filtered.txt", sep="\t",header = TRUE)
dat  <- removeFirstColumn(data)
dat  <- processData(dat)


plotBoxPlotCorrelation <- function(dat){
    nm <- names(dat)
    plots <- list()  # new empty list
    for (i in nm) {
        newData <- data.frame(x = dat[1L], y = dat[i])
        colnames(newData) <- c("x", "y")
        #plot <- ggplot(x,aes_string(x = nm[i])) + geom_histogram(alpha = .5,fill = "mediumseagreen")
        plot <- ggplot(newData, aes(x=x, y=y, group = x, color = x, alpha = 0.2)) +
            geom_boxplot(outlier.shape = 1, outlier.size = 0.5, notch = FALSE, notchwidth = 0.5, outlier.colour = "black") +
            xlab(nm[1L]) +
            ylab(i) +
            theme_bw() +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
        plots[[i]] <- plot  # add each plot into plot list
    }
    return (plots)
}

plots <- plotBoxPlotCorrelation(dat)
multiplot(plotlist = plots, cols = 4)


round(dat$EXPERIMENTAL, digits = 1)

breaks <- seq(1, 10, by = .5)
plotData <- ggplot(dat, aes(x=EXPERIMENTAL, y=COFACTOR, group=EXPERIMENTAL))
plotData <- plotData + geom_boxplot()
#plotData <- plotData + scale_x_discrete(name= "Experimental", breaks=breaks)
plotData  # add each plot into plot list


