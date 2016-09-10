#' pIIterativePK
#'
#' This data contains all the pk For the iterative method using different methods "solomon", "rodwell","emboss", "lehninger", "grimsley", "patrickios","DtaSelect","toseland","thurlkill","nozaki_tanford"
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
#'   \item TOSELAND
#'   \item THURLKILL
#'   \item NOZAKI_TANFORD
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
#' @param gamma The parameter to fix precision
#'
#' @examples
#' pIIterative(sequence="GLPRKILCAIAKKKGKCKGPLKLVCKC", pkSetMethod = "solomon")
#'

pIIterative <- function(sequence, pkSetMethod = "solomon", gamma = 0.001){
  #reformat sequence
  sequence <- reformat(seq= sequence)
  
  #load requeried pK set
  pkSet <- loadPkSetIterative(pkSetMethod)
  
  pH = 7.0
  
  #compute pI at reducied pH range 0.0-7.0. It avoid to
  #use the whole pH range.
  if(chargeAtPH(sequence, pH, pkSet) < 0){
    pHs <- seq(0, 7, gamma)
    charges <- chargeAtPH(sequence, pHs, pkSet)
    pH <- pHs[which.min(abs(charges))]
    return(specify_decimal(pH, 4))
  }
  
  #compute pI at reducied pH range 7.0-14.0. It avoid to
  #use the whole pH range.
  if(chargeAtPH(sequence, pH, pkSet) > 0){
    pHs <- seq(7, 14, gamma)
    charges <- chargeAtPH(sequence, pHs, pkSet)
    pH <- pHs[which.min(abs(charges))]
    return(specify_decimal(pH, 4))
  }
  
    #them return current pH value
    #pH <-specify_decimal(pH,4)
    return (pH)
}


#' pIIterativeMultipleSequences
#'
#' This function compute the isoelectric point of sequences contained into dataframe using the Iterative method.
#' At lest, a column named "sequence" is expected.
#' It return a dataframe with predicted pI value binded.
#' 
#' 
#' @param df The dataframe with sequences
#' @param pKSetMethod The pK set method
#'

pIIterativeMultipleSequences <- function(sequences = df, pkSetMethod = "solomon"){
  
  #getting column with sequences
  data <- subset(sequences, select=c("sequence"))
  
  #compute pI
  data <- mdply(data, function(sequence) { pIIterative(sequence = sequence, pkSetMethod) })
  
  #rename column added with pK set name
  names(data)[names(data)=="V1"] <- c(pkSetMethod)
  
  return(data)
  
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
    toseland <- pIIterative(sequence = seq, pkSetMethod = "toseland")
    thurlkill <- pIIterative(sequence = seq, pkSetMethod = "thurlkill")
    nozaki_tanford <- pIIterative(sequence = seq, pkSetMethod = "nozaki_tanford")
    
    values <- data.frame(method=c("solomon", "rodwell","emboss", "lehninger", "grimsley", "patrickios","DtaSelect","toseland","thurlkill","nozaki_tanford"), values=c(solomon, rodwell,emboss, lehninger, grimsley, patrickios,DtaSelect,toseland,thurlkill,nozaki_tanford))

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

     sequence <- gsub("[\r\n ]", "", sequence)
  
     lev <- c("n","m","p","R","H","K","D","E","C","Y") #here the entries that will be considered (PTMs, aa basic, aa acid)
  
     aa <- table(factor(prot <- strsplit(sequence, "")[[1]], levels = lev))
  
     #To evaluate N-terminal contribution
     if(prot[1]=="n" || prot[1]=="m"){       #n:Acetylation, m:Acetylation+Oxidation
        nterm <- 0.0
      } else {
        nterm <- (1/(1 +  10^(1 *  (pH - retrievePKValue("NTerm", pKIterative)))))      #otherwise
      }
  
     #To evaluate C-terminal contribution
     cterm <- (-1/(1 + 10^(-1 * (pH - retrievePKValue("CTerm", pKIterative)))))
  
     carg <- aa["R"] * (1/(1 + 10^(1 * (pH - retrievePKValue("R", pKIterative)))))
     chis <- aa["H"] * (1/(1 + 10^(1 * (pH - retrievePKValue("H", pKIterative)))))
     clys <- aa["K"] * (1/(1 + 10^(1 * (pH - retrievePKValue("K", pKIterative)))))
     casp <- aa["D"] * (-1/(1 + 10^(-1 * (pH - retrievePKValue("D", pKIterative)))))
     cglu <- aa["E"] * (-1/(1 + 10^(-1 * (pH - retrievePKValue("E", pKIterative)))))
     ccys <- aa["C"] * (-1/(1 + 10^(-1 * (pH - retrievePKValue("C", pKIterative)))))
     ctyr <- aa["Y"] * (-1/(1 + 10^(-1 * (pH - retrievePKValue("Y", pKIterative)))))
  
     #computing phosphorylation contribution
     pcharge <- 0
      if(grepl(pattern = "p", x=sequence, fixed = "TRUE")){  
        for (i in 1:length(prot)){
          if(prot[i] == "p"){              
            phosphoAA <- prot[i+1]         #getting the next aminoacid (Phospho) in the sequence...
            pcharge = pcharge + pchargePhosphorylation(phosphoAA, pH)
          }
         }
       }
  
      charge <- as.numeric(nterm + carg + clys + chis + casp + cglu + ctyr + ccys + cterm + pcharge)
  
      return(charge)
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


