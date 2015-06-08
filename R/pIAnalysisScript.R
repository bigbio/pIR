
data <- read.table("/Users/yperez/IdeaProjects/predictedAnalysis/data/Protein-PI-Filtered.txt", sep="\t",header = TRUE)
dat  <- removeFirstColumn(data)
dat  <- processData(dat)
#plots <- plotRawCorrelation(dat)
#multiplot(plotlist = plots, cols = 4)




