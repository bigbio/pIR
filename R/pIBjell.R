#' computeAllBjellValues
#'
#' This function compute the isoelectric point for all the pK sets
#' @param seq
#'

computeAllBjellValues <- function(seq){
    expasy     <- pIBjell(sequence = seq, pkSetMethod = "expasy")
    skoog      <- pIBjell(sequence = seq, pkSetMethod = "skoog")
    calibrated <- pIBjell(sequence = seq, pkSetMethod =  "calibrated")
    bjell      <- pIBjell(sequence = seq, pkSetMethod = "bjell")
    values <- data.frame(method=c("expasy", "skoog","calibrated", "bjell"), values=c(expasy, skoog, calibrated, bjell))
    colnames(values) <- c("method", "values")
    return(values)

}

#' pIBjell
#'
#' This function compute the isoelectric point of ami anocid sequences using the Bjell method
#' @param sequence The amino acid sequence
#' @param pKSetMethod The pK set method
#'
pIBjell <- function(sequence, pkSetMethod = "expasy"){
    sequence <- reformat(seq= sequence)
    NtermPK <- loadNTermPK(pkSet = pkSetMethod)
    CtermPK <- loadCTermPK(pkSet= pkSetMethod)
    GroupPK <- loadGroupPK(pkSet = pkSetMethod)

    pH  <- 6.5         # Starting point pI = 6.5 - theoretically it should be 7, but
    # Average protein pI is 6.5 so we increase the probability.
    lastCharge <- 0
    gamma      <- 0.00001
    this.step = 3.5

    repeat {
        charge <- getcharge(sequence = sequence, NTermPK = NtermPK, CTermPK = CtermPK, GroupPK = GroupPK, pH = pH)
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
    #print(pH)
    return (pH)
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

    aaV <- strsplit(sequence, "", fixed = TRUE)
    aaNTerm <- aaV[[1]][1]
    aaCTerm <- aaV[[1]][nchar(sequence)]


    pHpK <- 0.0
    this.FoRmU <- 0.0

    #To evaluate N-terminal contribution
    if(aaNTerm=="n" || aaNTerm=="m"){       #n:Acetylation, m:Acetylation+Oxidation
        pHpK <- pH - 0
        this.FoRmU <- this.FoRmU + (1.0 / (1.0 + (10.0^pHpK)))

    } else if(aaNTerm=="o"){       #o: Oxidation, it no modifies N-terminal contribution
        pHpK <- pH - retrievePKValue(aaV[[1]][2], NTermPK) #to eliminate PTM (neutral) to compute N-terminal contribution
        this.FoRmU <- this.FoRmU + (1.0 / (1.0 + (10.0^pHpK)))

    } else {
        pHpK <- pH - retrievePKValue(aaNTerm, NTermPK)          #otherwise
        this.FoRmU <- this.FoRmU + (1.0 / (1.0 + (10.0^pHpK)))
    }



    pHpK <- retrievePKValue(aaCTerm, CTermPK) - pH
    this.FoRmU = this.FoRmU + (-1.0 / (1.0 + (10.0^pHpK)))


    for(i in 1:nchar(sequence)) {

        aa <- aaV[[1]][i]

        if (!is.na(retrievePKValue(aa, GroupPK))){
            valuepK <- retrievePKValue(aa, GroupPK)

            if (valuepK < 0.0) {
                pHpK <- pH + valuepK;
                this.FoRmU <- this.FoRmU + (1.0 / (1.0 + (10.0^pHpK)))
            } else {
                pHpK <- valuepK - pH;
                this.FoRmU <- this.FoRmU + (-1.0 / (1.0 + (10.0^pHpK)))
            }
        }

        if(aa == "p"){            #computing phosphorylation contribution
            aa <- aaV[[1]][i+1]     #getting the next amino acid in the sequence...
            this.FoRmU = this.FoRmU + pchargePhosphorylation(aa, pH)
        }

    }
    return (this.FoRmU)
}

