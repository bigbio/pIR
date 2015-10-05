#' plotPTMCorrelation
#'
#' This function plot the PTMs correlation between to Files, one including the PTMs and another one without considering it
#'
#' @param pepFile1 with PTMs
#' @param pepFile2 without PTMs
#' @param height of the chart
#' @param width of the chart
#' @param cols number of colums
#'
plotPTMCorrelation <- function(pepFile1, pepFile2, height = 800, width = 800, cols = 3){

    # Read the data files. File1: excluding modification. File2: including modification
    pepPTM1 <- read.table(file=pepFile1,  header = TRUE, sep = "\t")
    pepPTM2 <- read.table(file=pepFile2,  header = TRUE, sep = "\t")

    # Remove the empty values and other data, also remove the string data Proteins and Peptides
    datPE1 <- removeFirstColumn(pepPTM1)
    datPE2 <- removeFirstColumn(pepPTM2)

    datPE1 <- processData(datPE1)
    datPE2 <- processData(datPE2)

    plotsPTM <- plotPTMShiftOnFractions(datPE1, datPE2)


    png("plotPtmOnFraction.png", width = 800, height = 800)
    multiplot(plotlist = plotsPTM, cols=3)
    dev.off()

    plotsPTM <- plotPTMShiftOverallData(datPE1, datPE2)

    png("plotPtmOverall.png", width = 800, height = 800)
    multiplot(plotlist = plotsPTM, cols=3)
    dev.off()

}
