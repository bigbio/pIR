#' pIIterativePK
#'
#' This data contains all the pk For the iterative method using different methods
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
#' J Proteomics. 2012 Apr 3;75(7):2269-74. doi: 10.1016/j.jprot.2012.01.029. Epub 2012 Feb 3.
#'
#' @param sequence the amino acid sequence
#' @param paramter
#'
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
        if(error > gamma){
            break;
        }
    }

    return (pH)
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

    for (i in 0:nchar(sequence)){
        aa<-c(factor(prot<-strsplit(toupper(sequence),"")[[1]], levels = LETTERS))
        if(any(AAAcid == aa) && !is.na(retrievePKValue(aa, pKIterative))){
            charge <- charge + pcharge(pH, retrievePKValue(aa, pKIterative))
        }
        if(any(AABasid == aa) && !is.na(retrievePKValue(aa, pKIterative))){
            charge = charge - pcharge(retrievePKValue(aa, pKIterative),pH)
        }
    }
    return (charge)

}

pcharge <- function(pH, pk){
    val <- 10^ (pH - pk)
    val <- 1/(1 + val)
    return (val)
}

loadPkSet <- function(pkSetMethod = "solomon"){

    pKValues <- c()
    if(pkSetMethod == "solomon"){
        pkValues <- data.frame(key=c("CTerm", "NTerm", "D", "E", "K", "R", "H", "C", "Y"), c(2.4, 9.6,3.9,4.3,10.5,12.5,6.0,8.3,10.1))
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
#'
retrievePKValue <- function(aa, pKIterative){
    pkValue <- NA
    for(i in 1:nrow(pKIterative)){
        print(pKIterative[i,"key"])
        if(pKIterative[i,"key"] == aa){
            pkValue <- pKIterative[i,"value"]
            print(pKIterative[i,"value"])
        }
    }
    print(pkValue)
    return (pkValue)
}

#pi <- pIIterative(sequence = "GLPRKILCAIAKKKGKCKGPLKLVCKC", pkSetMethod = "solomon")


