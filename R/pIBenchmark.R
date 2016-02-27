#' Isoelectric poit of a set of peptides from a published paper
#'
#' This dataset contains more than 6000 peptides with the experimental isoelectric point and the predicted pI using different methods.
#' The variables are as follows:
#'
#' \itemize{
#'   \item SEQUENCE. The potein accession
#'   \item EXPERIMENTAL. The experimental isoelectric point
#'   \item COFACTOR. The isoelectric point using the method developed by Cargile and co-workers.
#'   \item ITERATIVE_THURLKILL
#'   \item ITERATIVE_TOSELAND
#'   \item ITERATIVE_GRIMSLEY
#'   \item ITERATIVE_SOLOMON
#'   \item ITERATIVE_EMBOSS
#'   \item ITERATIVE_LEHNINGER
#'   \item ITERATIVE_PATRICKIOS
#'   \item ITERATIVE_RICHARD
#'   \item ITERATIVE_RODWELL
#'   \item ITERATIVE_SILLERO
#'   \item BRANCA
#'   \item BJELL_DEFAULT
#'   \item BJELL_SKOOG
#'   \item BJELL_CALLIBRATED
#'   \item BJELL_EXPASY
#'   \item SVM
#' }
#'
#' @docType data
#' @keywords datasets
#' @name peptides
#' @usage data(peptides)
#' @format A data frame with 1041 rows and 19 variables
NULL

#' Isoelectric poit of a set of proteins from PID-DB
#'
#' This dataset contains more than 1000 proteins with the experimental isoelectric point and the predicted pI using different methods.
#' The variables are as follows:
#'
#' \itemize{
#'   \item PROTEIN. The potein accession
#'   \item EXPERIMENTAL. The experimental isoelectric point
#'   \item COFACTOR. The isoelectric point using the method developed by Cargile and co-workers.
#'   \item ITERATIVE_THURLKILL
#'   \item ITERATIVE_TOSELAND
#'   \item ITERATIVE_GRIMSLEY
#'   \item ITERATIVE_SOLOMON
#'   \item ITERATIVE_EMBOSS
#'   \item ITERATIVE_LEHNINGER
#'   \item ITERATIVE_PATRICKIOS
#'   \item ITERATIVE_RICHARD
#'   \item ITERATIVE_RODWELL
#'   \item ITERATIVE_SILLERO
#'   \item BRANCA
#'   \item BJELL_DEFAULT
#'   \item BJELL_SKOOG
#'   \item BJELL_CALLIBRATED
#'   \item BJELL_EXPASY
#'   \item SVM
#' }
#'
#' @docType data
#' @keywords datasets
#' @name proteins
#' @usage data(proteins)
#' @format A data frame with 1041 rows and 19 variables
NULL

#' pIBenchmark
#'
#' This function get the original file with peptides and proteins and do the complete benchmark. At the end it generates a file with all the comparison
#'
#' @param pepFile the peptide file
#' @param protFile the protein file
#' @param height
#' @param width
#' @param cols  number of columns by chart

pIBenchmark <- function(pepFile, protFile, height = 800, width = 800, cols = 3){
    # Read the data files
    peptides <- read.table(file=pepFile, header = TRUE, sep = "\t")
    proteins <- read.table(file=protFile, header = TRUE, sep = "\t")

    # Remove the empty values and other data, also remove the string data Proteins and Peptides
    datPE <- removeFirstColumn(peptides)
    datPR <- removeFirstColumn(proteins)

    datPE <- processData(datPE)
    datPR <- processData(datPR)

    distributionPeptides <- plotHistFunc(datPE, na.rm = TRUE)
    pdf("peptideDistributions.png", width = width, height = height)
    multiplot(plotlist = distributionPeptides, cols = cols)
    dev.off()

    distributionProteins <- plotHistFunc(datPR, na.rm = TRUE)
    png("proteinsDistributions.png", width = width, height = height)
    multiplot(plotlist = distributionProteins, cols=cols)
    dev.off()

    boxPlotCorrelationPeptides <- plotBoxPlotCorrelation(datPE)
    png("peptideBoxPlotCorrelation.png", width = width, height = height)
    multiplot(plotlist = boxPlotCorrelationPeptides, cols=cols)
    dev.off()

    rawCorrelationPeptides <- plotRawCorrelation(datPE)
    png("peptideRawCorrelation.png", width = width, height = height)
    multiplot(plotlist = rawCorrelationPeptides, cols=cols)
    dev.off()

    rawCorrelationProteins <- plotRawCorrelation(datPR)
    png("proteinRawCorrelation.png", width = width, height = height)
    multiplot(plotlist = rawCorrelationProteins, cols=cols)
    dev.off()

    resultStatsPE <- bindRMSECorrelationFrame(datPE, method = "pearson")
    write.csv(resultStatsPE, file = "peptideStats.csv")

    resultStatsPR <- bindRMSECorrelationFrame(datPR, method = "pearson")
    write.csv(resultStatsPR, file = "proteinStats.csv")

}



#' computePIvalues
#'
#' This function take a data frame in the way of: sequence pIExp, and return a new 
#' data frame with all pI values computed using multiple functions. This operation could be time-consuming.
#' 
#' @param originalData The original data frame
#'


computePIvalues <- function(originalData){
  
  data <- originalData
  
  colnames(data) <-c("sequence", "pIExp")
  
  #Add all the bjell methods and pk Sets
  data <- mdply(data, function(sequence, pIExp) { pIBjell(sequence = sequence, pkSetMethod = "bjell") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell")
  
  data <- mdply(data, function(sequence, pIExp, bjell) { pIBjell(sequence = sequence, pkSetMethod = "expasy") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy")
  
  data <- mdply(data, function(sequence, pIExp, bjell, expasy) { pIBjell(sequence = sequence, pkSetMethod = "skoog") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog")
  
  #Add aaindex attribute
  data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog) { aaIndex(sequence = sequence) })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex")
  
  #Add all the Iterative methods and pk Sets
  
  #add solomon
  data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, aaindex) { pIIterative(sequence = sequence, pkSetMethod = "solomon") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon")
  
  #add rodwell
  data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, aaindex, solomon) { pIIterative(sequence = sequence, pkSetMethod = "rodwell") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon", "rodwell")
  
  #add emboss
  data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, aaindex, solomon, rodwell) { pIIterative(sequence = sequence, pkSetMethod = "emboss") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon", "rodwell", "emboss")
  
  #add lehninger
  data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, aaindex, solomon, rodwell, emboss) { pIIterative(sequence = sequence, pkSetMethod = "lehninger") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon", "rodwell", "emboss", "lehninger")
  
  #add grimsley
  data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, aaindex, solomon, rodwell, emboss, lehninger) { pIIterative(sequence = sequence, pkSetMethod = "grimsley") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon", "rodwell", "emboss", "lehninger", "grimsley")
  
  #add patrickios
  data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, aaindex, solomon, rodwell, emboss, lehninger, grimsley) { pIIterative(sequence = sequence, pkSetMethod = "patrickios") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios")
  
  #add DtaSelect
  data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, aaindex, solomon, rodwell, emboss, lehninger, grimsley, patrickios) { pIIterative(sequence = sequence, pkSetMethod = "DtaSelect") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios", "DtaSelect")
  
  #add toseland
  data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, aaindex, solomon, rodwell, emboss, lehninger, grimsley, patrickios, DtaSelect) { pIIterative(sequence = sequence, pkSetMethod = "toseland") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios", "DtaSelect", "toseland")
  
  #add thurlkill
  data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, aaindex, solomon, rodwell, emboss, lehninger, grimsley, patrickios, DtaSelect, toseland) { pIIterative(sequence = sequence, pkSetMethod = "thurlkill") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios", "DtaSelect", "toseland", "thurlkill")
  
  #add nozaki_tanford
  data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, aaindex, solomon, rodwell, emboss, lehninger, grimsley, patrickios, DtaSelect, toseland, thurlkill) { pIIterative(sequence = sequence, pkSetMethod = "nozaki_tanford") })
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios", "DtaSelect", "toseland", "thurlkill", "nozaki_tanford")
  
  #add cofactor
  #data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, aaindex, solomon, rodwell, emboss, lehninger, grimsley, patrickios, DtaSelect, toseland, thurlkill, nozaki_tanford) { pICofactor(sequence = sequence) })
  
  #colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios", "DtaSelect", "toseland", "thurlkill", "nozaki_tanford", "cofactor")
  
  #Add SVM pI data
  #warning: the dataframe argument must containt the attributes: bjell, expasy and aaindex.
  data <- pISVMsequences(dataframe = data, defaultModel = FALSE) 
  
  colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios", "DtaSelect", "toseland", "thurlkill", "nozaki_tanford", "pISVM")
  
  #saving data with all variables
  write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")
  
  #print(data)
  
  return(data)
}