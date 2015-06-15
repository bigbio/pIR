#' pIIterativePK
#'
#' This data contains all the pk For the iterative method using different methods "solomon", "rodwell","emboss", "lehninger", "grimsley", "patrickios","DtaSelect"
#'
#' The variables are as follows:
#'
#' \itemize{
#'   \item AA Aminoacid
#'   \item SOLOMON
#'   \item SILLERO
#'   \item EMBOSS
#'   \item RODEWELL
#'   \item PATRICKIOS
#'   \item RICHARD
#'   \item LEHNINGER
#' }
#'
#' @docType data
#' @keywords datasets
#' @name pIIterativePK
#' @usage data(pIIterativePK)
#' @format A data frame with 8 rows and 9 variables
NULL

#' pIIterative
#'
#' This function compute the isoelectric point of a sequence using the Iterative method
#' Perez-Riverol et al. Isoelectric point optimization using peptide descriptors and support vector machines. J Proteomics. 2012 Apr 3;75(7):2269-74.
#'
#' @param sequence The amino acid sequence to be sue to compute the Isoelectric point
#' @param pkSetMethod The pK set to be used
#'
#' @examples
#' pIIterative(sequence="GLPRKILCAIAKKKGKCKGPLKLVCKC", pkSetMethod = "solomon")
#'

pIIterative <- function(sequence, pkSetMethod = "solomon"){

    pkSet <- loadPkSetIterative(pkSetMethod)

    pH  <- 6.5         # Starting point pI = 6.5 - theoretically it should be 7, but
                       # Average protein pI is 6.5 so we increase the probability.
    lastCharge <- 0
    gamma      <- 0.00001
    this.step = 3.5

    repeat {
        charge <- chargeAtPH(sequence, pH, pkSet)
        if( charge > 0){
            pH <- pH + this.step
        }else{
            pH <- pH - this.step
        }
        this.step = this.step/2
        error <- abs(charge-lastCharge)
        lastCharge <- charge
        if(error < gamma){
            break;
        }
    }
    pH <-specify_decimal(pH,3)
    return (pH)
}

#' computeAllIterativeValues
#'
#' This function compute the isoelectric point for all the pK sets
#' @param seq
#'

computeAllIterativeValues <- function(seq){
    solomon <- pIIterative(sequence = seq, pkSetMethod = "solomon")
    rodwell <- pIIterative(sequence = seq, pkSetMethod = "rodwell")
    emboss <- pIIterative(sequence = seq, pkSetMethod =  "emboss")
    lehninger <- pIIterative(sequence = seq, pkSetMethod = "lehninger")
    grimsley <- pIIterative(sequence = seq, pkSetMethod = "grimsley")
    patrickios <- pIIterative(sequence = seq, pkSetMethod = "patrickios")
    DtaSelect <- pIIterative(sequence = seq, pkSetMethod = "DtaSelect")

    values <- data.frame(method=c("solomon", "rodwell","emboss", "lehninger", "grimsley", "patrickios","DtaSelect"), values=c(solomon, rodwell,emboss, lehninger, grimsley, patrickios,DtaSelect))
    colnames(values) <- c("method", "values")
    return(values)

}
#' chargeAtPH
#'
#' This fucntion compute the charge of the peptide using a pkset
#'
#' @param sequence the aminoacid sequence
#' @param pH current PH
#' @param pKIterative the selected pK Set

chargeAtPH <- function(sequence, pH = 7, pKIterative){

    AAAcid  <- c("K","R", "H");
    AABasid <- c("D", "E", "C", "Y")

    charge <- 0.0

    NTerm_Pk <- retrievePKValue("NTerm", pKIterative)
    CTerm_Pk <- retrievePKValue("CTerm", pKIterative)

    charge <- charge + pcharge(pH, NTerm_Pk)
    charge <- charge - pcharge(CTerm_Pk,pH)
    sequence <- toupper(sequence)
    aaV <- strsplit(sequence, "", fixed = TRUE)
    for (i in 1:nchar(sequence)){
        aa <- aaV[[1]][i]
        if(any(AAAcid == aa) && !is.na(retrievePKValue(aa, pKIterative))){
            charge <- charge + pcharge(pH, retrievePKValue(aa, pKIterative))
        }
        if(any(AABasid == aa) && !is.na(retrievePKValue(aa, pKIterative))){

            charge = charge - pcharge(retrievePKValue(aa, pKIterative),pH)
        }
    }
    return (charge)
}

#' pcharge
#'
#' This pcharge function computhe the charge
#' @param pH current pH
#' @param pk current pk

pcharge <- function(pH, pk){
    val <- 10^(pH - pk)
    val <- 1/(1 + val)
    return (val)
}

#' specify_decimal
#'
#' This function round a value to an specific decimal places
#' @param x the double value
#' @param k the decimal places

specify_decimal <- function(x, k){
    return(format(round(x, k), nsmall=k))
}


