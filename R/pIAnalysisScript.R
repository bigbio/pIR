
data <- read.table("/Users/yperez/IdeaProjects/predictedAnalysis/data/Peptide-PI-Filtered.txt", sep="\t",header = TRUE)
dat  <- removeFirstColumn(data)
dat  <- processData(dat)


#plots <- plotBoxPlotCorrelation(dat)
#multiplot(plotlist = plots, cols = 3)

#plotData <- plotData + scale_x_discrete(name= "Experimental", breaks=breaks)
#plotData  # add each plot into plot list

corr <- computeCorrelationExperimental(dat, method = "pearson")
rmseValue <- computeRMSEExperimental(dat)



statsResults <- matrix(ncol = 3, nrow = ncol(rmseValue)+1)

for( i in colnames(rmseValue)){
}

