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

    pkSet <- loadPkSet(pkSetMethod)

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

#' computeAllPKValues
#'
#' This function compute the isoelectric point for all the pK sets
#' @param seq
#'

computeAllPKValues <- function(seq){
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

#' loadPkSet
#'
#' This function retrieve a pk dataset using the name of the pKSet
#' \itemize{
#'   \item sillero : Sillero, A., Maldonado, A. (2006) Isoelectric point determination of proteins and other macromolecules: oscillating method. Comput Biol Med., 36:157-166.
#'   \item rodwell    : Rodwell, J. Heterogeneity of component bands in isoelectric focusing patterns. Analytical Biochemistry, 1982, 119 (2), 440-449.
#'   \item solomon    : Solomon, T.W.G. (1998) Fundamentals of Organic Chemistry, 5th edition. Published by Wiley.
#'   \item emboss     : EMBOSS data are from http://emboss.sourceforge.net/apps/release/5.0/emboss/apps/iep.html.
#'   \item lehninger  : Nelson, David L., Albert L. Lehninger, and Michael M. Cox. Lehninger principles of biochemistry. Macmillan, 2008.
#'   \item patrickios : Patrickios, Costas S., and Edna N. Yamasaki. "Polypeptide amino acid composition and isoelectric point ii. comparison between experiment and theory." Analytical biochemistry 231.1 (1995): 82-91.
#'   \item grimsley   : Gerald R Grimsley, J Martin Scholtz and C Nick Pace. A summary of the measured pK values of the ionizable groups in folded proteins. Protein Sci. 2009 Jan; 18(1): 247–251.
#'   \item wikipedia  : http://en.wikipedia.org/wiki/List_of_standard_amino_acids
#'   \item DtaSelect  : http://fields.scripps.edu/DTASelect/20010710-pI-Algorithm.pdf
#' }
#'
#' @param pkSetMethod name of the pk
#'

loadPkSet <- function(pkSetMethod = "solomon"){
    pKValues <- c()
    if(pkSetMethod == "solomon"){
        pkValues <- data.frame(key=c("NTerm","CTerm","C","D","E","H","K","R","Y"), c(9.6,2.4,8.3,3.9,4.3,6.0,10.5,12.5,10.1))
        colnames(pkValues) <- c("key", "value")

    }
    if(pkSetMethod == "rodwell"){
        pkValues <- data.frame(key=c("NTerm","CTerm","C","D","E","H","K","R","Y"), c(8.0,3.1,8.33,3.68,4.25,6.0,11.5,11.5,10.07))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSetMethod == "emboss"){
        pkValues <- data.frame(key=c("NTerm","CTerm","C","D","E","H","K","R","Y"), c(8.6,3.6,8.5,3.9,4.1,6.5,10.8,12.5,10.1))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSetMethod == "lehninger"){
        pkValues <- data.frame(key=c("NTerm","CTerm","C","D","E","H","K","R","Y"), c(9.69,2.34,8.33,3.86,4.25,6.0,10.5,12.4,10.0))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSetMethod == "grimsley"){
        pkValues <- data.frame(key=c("NTerm","CTerm","C","D","E","H","K","R","Y"), c(7.7,3.3,6.8,3.5,4.2,6.6,10.5,12.04,10.3))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSetMethod == "patrickios"){
        pkValues <- data.frame(key=c("NTerm","CTerm","C","D","E","H","K","R","Y"), c(11.2,4.2,0.0,4.2,4.2,0.0,11.2,11.2,0.0))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSetMethod == "DtaSelect"){
        pkValues <- data.frame(key=c("NTerm","CTerm","C","D","E","H","K","R","Y"), c(8.0,3.1,8.5,4.4,4.4,6.5,10.0,12.0,10.0))
        colnames(pkValues) <- c("key", "value")
    }
    return (pkValues)
}

#' retrievePKValue
#'
#' This function retrieve the pK values for an specific aminoacid
#'
#' @param aminoacid
#' @param pKIterative

retrievePKValue <- function(aa, pKIterative){
    pkValue <- NA
    for(i in 1:nrow(pKIterative)){
        if(pKIterative[i,"key"] == aa){
            pkValue <- pKIterative[i,"value"]
        }
    }
    return (pkValue)
}

#' specify_decimal
#'
#' This function round a value to an specific decimal places
#' @param x the double value
#' @param k the decimal places

specify_decimal <- function(x, k){
    return(format(round(x, k), nsmall=k))
}


