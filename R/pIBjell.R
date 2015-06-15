
#' pIBjell
#'
#' This function compute the isoelectric point of ami anocid sequences using the Bjell method
#' @param sequence The amino acid sequence
#' @param pKSetMethod The pK set method
#'
pIBjell <- function(sequence, pkSetMethod = "expasy"){

    piInt = 0.5
    this.pH = -1.0

    NtermPK <- loadNTermPK(pkSet = pkSetMethod)
    print(NtermPK)

    CtermPK <- loadCTermPK(pkSet= pkSetMethod)
    print(CtermPK)

    GroupPK <- loadGroupPK(pkSet = pkSetMethod)
    print(GroupPK)

    # This algorithm used this strategy: Take an of step = 0.5 and make a loop while
    # the charge of the SequenceAA was >= 0.0 then take the last value of the charge and the last
    # value of ph to make an exaustive step with and step of 0.0001.

    repeat{
        if (getcharge(sequence = sequence, NTermPK = NtermPK, CTermPK = CtermPK, GroupPK = GroupPK, pH = this.pH) < 0.0){
            break
            }
        this.pH = this.pH + piInt
        if (this.pH > 14.0){
                break
        }
        }

    this.pH = this.pH - piInt;
    piInt = 0.001

    # Take a short step to compute the pI in this range.
    # This algorithm is very exaustive because the error
    # in the algorithm is in the 3t decimal of the number.

    pIActual = 0.0;

    repeat{
        pIActual = getcharge(sequence = sequence, NTermPK = NtermPK, CTermPK = CtermPK, GroupPK = GroupPK, pH = this.pH);
        if (pIActual < 0.0){
            this.pH = this.pH - piInt
            break
            }
        this.pH = this.pH + piInt

        if(this.pH > 14.0){
            break
        }
        }

    if (pIActual >= 0.0) this.pH = 100.001;

    # this is an unreal value.

    pHround = round(this.pH * 100.0)

    this.pH = -1.0

    return (pHround / 100.0)
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

    sequence <- toupper(sequence)

    aaV <- strsplit(sequence, "", fixed = TRUE)
    aaNTerm <- aaV[[1]][1]
    aaCTerm <- aaV[[1]][nchar(sequence)]


    pHpK <- 0.0
    this.FoRmU <- 0.0

    pHpK <- pH - retrievePKValue(aaNTerm, NTermPK)

    this.FoRmU <- this.FoRmU + (1.0 / (1.0 + (10.0^pHpK)))

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
        }
    return (this.FoRmU)
}


#pi < - pIBjell(sequence = seq, pkSetMethod = "expasy")

