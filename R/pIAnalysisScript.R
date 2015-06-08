
data <- read.table("/Users/yperez/IdeaProjects/predictedAnalysis/data/Peptide-PI-Filtered.txt", sep="\t",header = TRUE)
dat  <- removeFirstColumn(data)
dat  <- processData(dat)


plotBoxPlotCorrelation <- function(dat){

}

plots <- plotBoxPlotCorrelation(dat)
multiplot(plotlist = plots, cols = 3)


round(dat$EXPERIMENTAL, digits = 1)

breaks <- seq(1, 10, by = .5)
plotData <- ggplot(dat, aes(x=EXPERIMENTAL, y=COFACTOR, group=EXPERIMENTAL))
plotData <- plotData + geom_boxplot()
#plotData <- plotData + scale_x_discrete(name= "Experimental", breaks=breaks)
plotData  # add each plot into plot list


