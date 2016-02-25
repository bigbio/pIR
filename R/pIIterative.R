#' pIIterativePK
#'
<<<<<<< HEAD
#' This data contains all the pk For the iterative method using different methods "solomon", "rodwell","emboss", "lehninger", "grimsley", "patrickios","DtaSelect","toseland","thurlkill","nozaki_tanford" 
=======
#' This data contains all the pk For the iterative method using different methods "solomon", "rodwell","emboss", "lehninger", "grimsley", "patrickios","DtaSelect"
>>>>>>> ff7e566e64c270577ddf07143ee6d0751213acbe
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
<<<<<<< HEAD
#'   \item TOSELAND
#'   \item THURLKILL
#'   \item NOZAKI_TANFORD
=======
<<<<<<< HEAD
#'   \item TOSELAND
#'   \item THURLKILL
#'   \item NOZAKI_TANFORD
=======
>>>>>>> ff7e566e64c270577ddf07143ee6d0751213acbe
>>>>>>> 0b2156cbfde9fe3414f7daa1f85b32c9d1e12bf0
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
    sequence <- reformat(seq= sequence)
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
    pH <-specify_decimal(pH,4)
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
<<<<<<< HEAD
    toseland <- pIIterative(sequence = seq, pkSetMethod = "toseland")
    thurlkill <- pIIterative(sequence = seq, pkSetMethod = "thurlkill")
    nozaki_tanford <- pIIterative(sequence = seq, pkSetMethod = "nozaki_tanford")
    values <- data.frame(method=c("solomon", "rodwell","emboss", "lehninger", "grimsley", "patrickios","DtaSelect","toseland","thurlkill","nozaki_tanford"), values=c(solomon, rodwell,emboss, lehninger, grimsley, patrickios,DtaSelect,toseland,thurlkill,nozaki_tanford))
=======
<<<<<<< HEAD
    toseland <- pIIterative(sequence = seq, pkSetMethod = "toseland")
    thurlkill <- pIIterative(sequence = seq, pkSetMethod = "thurlkill")
    nozaki_tanford <- pIIterative(sequence = seq, pkSetMethod = "nozaki_tanford")

    values <- data.frame(method=c("solomon", "rodwell","emboss", "lehninger", "grimsley", "patrickios","DtaSelect","toseland","thurlkill","nozaki_tanford"), values=c(solomon, rodwell,emboss, lehninger, grimsley, patrickios,DtaSelect,toseland,thurlkill,nozaki_tanford))
=======

    values <- data.frame(method=c("solomon", "rodwell","emboss", "lehninger", "grimsley", "patrickios","DtaSelect"), values=c(solomon, rodwell,emboss, lehninger, grimsley, patrickios,DtaSelect))
>>>>>>> ff7e566e64c270577ddf07143ee6d0751213acbe
>>>>>>> 0b2156cbfde9fe3414f7daa1f85b32c9d1e12bf0
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

    aaV <- strsplit(sequence, "", fixed = TRUE)

    #To evaluate N-terminal contribution
    if(aaV[[1]][1]=="n" || aaV[[1]][1]=="m"){       #n:Acetylation, m:Acetylation+Oxidation
        NTerm_Pk <- 0.0
        charge <- charge + pcharge(pH, NTerm_Pk)
    } else {
        NTerm_Pk <- retrievePKValue("NTerm", pKIterative)      #otherwise
        charge <- charge + pcharge(pH, NTerm_Pk)
    }

    CTerm_Pk <- retrievePKValue("CTerm", pKIterative)
    charge <- charge - pcharge(CTerm_Pk,pH)

    #sequence <- toupper(sequence) #convert any char to uppercase

    for (i in 1:nchar(sequence)){
        aa <- aaV[[1]][i]
        if(any(AAAcid == aa) && !is.na(retrievePKValue(aa, pKIterative))){
            charge <- charge + pcharge(pH, retrievePKValue(aa, pKIterative))
        }
        if(any(AABasid == aa) && !is.na(retrievePKValue(aa, pKIterative))){

            charge = charge - pcharge(retrievePKValue(aa, pKIterative),pH)
        }
        if(aa == "p"){              #computing phosphorylation contribution
            aa <- aaV[[1]][i+1]     #getting the next amino acid in the sequence...
            charge = charge + pchargePhosphorylation(aa, pH)
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

#' pchargePhosphorylation
#'
#' This pcharge function computhe the charge contribution of phosphorylation for S, T and Y amino acid
#' @param aa amino acid phosphorylated
#' @param pH current pH
#'

pchargePhosphorylation <- function(aa, pH){

    val <- 0

    if(!is.na(aa)){

        if(aa == "S" || aa == "T"){

            STpKa1 <- 1.2       #pk values for phospho-amino S and T (from ProMoST tool)
            STpKa2 <- 6.5

            val_1 <- 10^(STpKa1 - pH)
            val_1 <- -1/(1 + val_1)

            val_2 <- 10^(STpKa2 - pH)
            val_2 <- -1/(1 + val_2)

            val <- val_1 + val_2

            return (val)
        }

        if(aa == "Y"){

            STpKa1 <- 1.2       #!!!!put here the specific pKa and pKb for Y
            STpKa2 <- 6.5

            val_1 <- 10^(STpKa1 - pH)
            val_1 <- -1/(1 + val_1)

            val_2 <- 10^(STpKa2 - pH)
            val_2 <- -1/(1 + val_2)

            val <- val_1 + val_2

            return (val)
        }

    }
    return (val)
}

#' specify_decimal
#'
#' This function round a value to an specific decimal places
#' @param x the double value
#' @param k the decimal places

specify_decimal <- function(x, k){
    value <- format(round(x, k), nsmall=k)
    return(as.numeric(value))
}


