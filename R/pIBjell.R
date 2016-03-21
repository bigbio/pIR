#' computeAllBjellValues
#'
#' This function compute the isoelectric point for all the pK sets
#' @param seq
#'

computeAllBjellValues <- function(seq){
    expasy     <- pIBjell(sequence = seq, pkSetMethod = "expasy")
    skoog      <- pIBjell(sequence = seq, pkSetMethod = "skoog")
    calibrated <- pIBjell(sequence = seq, pkSetMethod = "calibrated")
    values <- data.frame(method=c("expasy", "skoog","calibrated"), values=c(expasy, skoog, calibrated))
    colnames(values) <- c("method", "values")
    return(values)

}

#' pIBjell
#'
#' This function compute the isoelectric point of ami anocid sequences using the Bjell method
#' @param sequence The amino acid sequence
#' @param pKSetMethod The pK set method
#' @param gamma The parameter to fix precision
#'
pIBjell <- function(sequence, pkSetMethod = "expasy", gamma=0.001){
    sequence <- reformat(seq= sequence)
    NtermPK <- loadNTermPK(pkSet = pkSetMethod)
    CtermPK <- loadCTermPK(pkSet= pkSetMethod)
    GroupPK <- loadGroupPK(pkSet = pkSetMethod)

    pH = 7.0
    
    #compute pI at reducied pH range 0.0-7.0. It avoid to
    #use the whole pH range.
    if(getcharge(sequence, NtermPK, CtermPK, GroupPK, pH) < 0){
      pHs <- seq(0, 7, gamma)
      charges <- getcharge(sequence, NtermPK, CtermPK, GroupPK, pHs)
      pH <- pHs[which.min(abs(charges))]
      return(specify_decimal(pH, 4))
    }
    
    #compute pI at reducied pH range 7.0-14.0. It avoid to
    #use the whole pH range.
    if(getcharge(sequence,NtermPK, CtermPK, GroupPK, pH) > 0){
      pHs <- seq(7, 14, gamma)
      charges <- getcharge(sequence,NtermPK, CtermPK, GroupPK, pHs)
      pH <- pHs[which.min(abs(charges))]
      return(specify_decimal(pH, 4))
    }
    
    #them return current pH value
    return(pH) 
}

#' pIBjellMultipleSequences
#'
#' This function compute the isoelectric point of sequences contained into dataframe using the Bjell method.
#' At lest, a column named "sequence" is expected.
#' It return a dataframe with predicted pI value binded.
#' 
#' 
#' @param df The dataframe with sequences
#' @param pKSetMethod The pK set method
#'

pIBjellMultipleSequences <- function(sequences = df, pkSetMethod = "expasy"){
  
  #getting column with sequences
  data <- subset(sequences, select=c("sequence"))
  
  #compute pI
  data <- mdply(data, function(sequence) { pIBjell(sequence = sequence, pkSetMethod) })
  
  #rename column added with pK set name
  names(data)[names(data)=="V1"] <- c(pkSetMethod)
    
    return(data)
  
}


#' getcharge
#'
#' This function compute the charge to a SequenceAA
#'
#' @param sequence
#' @param NTermPK the N-term pk Set
#' @param CTermPK the C-term pk Set
#' @param GroupPK the Side Group pk Set
#' @param pH The actual pH
#'

getcharge <- function (sequence, NTermPK, CTermPK, GroupPK, pH){

    #sequence <- toupper(sequence)
    
    lev <- c("n","m","p","R","H","K","D","E","C","Y")  #here the entries that will be considered (PTMs, aa basic, aa acid)
   
    aaTable <- table(factor(prot <- strsplit(sequence, "")[[1]], levels = lev))
    
    aaNTerm <- prot[1]
    aaCTerm <- prot[length(prot)]
    
    #To evaluate N-terminal contribution
    if(aaNTerm=="n" || aaNTerm=="m"){       #n:Acetylation, m:Acetylation+Oxidation
      nterm <- 0.0
    } else if(aaNTerm=="o"){ 
      nterm <- (1/(1 +  10^(1 *  (pH - retrievePKValue(prot[2], NTermPK)))))   #if Met oxidated in N-terminal position avoid "o".
    }else{
      nterm <- (1/(1 +  10^(1 *  (pH - retrievePKValue(aaNTerm, NTermPK)))))   #otherwise
    }

    #To evaluate C-terminal contribution
    cterm <- (-1/(1 + 10^(-1 * (pH - retrievePKValue(aaCTerm, CTermPK)))))
    
    carg <- aaTable["R"] * (1/(1 + 10^(1 * (pH - retrievePKValue("R", GroupPK)))))
    chis <- aaTable["H"] * (1/(1 + 10^(1 * (pH - retrievePKValue("H", GroupPK)))))
    clys <- aaTable["K"] * (1/(1 + 10^(1 * (pH - retrievePKValue("K", GroupPK)))))
    casp <- aaTable["D"] * (-1/(1 + 10^(-1 * (pH - retrievePKValue("D", GroupPK)))))
    cglu <- aaTable["E"] * (-1/(1 + 10^(-1 * (pH - retrievePKValue("E", GroupPK)))))
    ccys <- aaTable["C"] * (-1/(1 + 10^(-1 * (pH - retrievePKValue("C", GroupPK)))))
    ctyr <- aaTable["Y"] * (-1/(1 + 10^(-1 * (pH - retrievePKValue("Y", GroupPK)))))
    
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

