#' pICofactor
#'
#' This function compute the isoelctric following the Cofactor method proposed by
#' @param sequence

pICofactor <- function(sequence)
    {

    pIAdjust_E_Nterm <- pICofactorAdjust_E_Nterm()
    pIAdjust_E_Cterm <- pICofactorAdjust_E_Cterm()
    pIAdjust_D_Nterm <- pICofactorAdjust_D_Nterm()
    pIAdjust_D_Cterm <- pICofactorAdjust_D_Cterm()
    pIAdjust_Cterminus <- pICofactorAdjust_Cterm()

    OldPI = 10.0
    HiPI  = 14.0
    LowPI = 0.0
    C_PI  = 7.0
    PItemp = 0

    while (abs(C_PI - OldPI) > 0.001){

        PItemp = 0
        sequence <- toupper(sequence)
        aaV <- strsplit(sequence, "", fixed = TRUE)
        for (i in 1:nchar(sequence)){
            aa <- aaV[[1]][i]
            ctermAA <- aaV[[1]][nchar(sequence)]
            ntermAA <- aaV[[1]][1]
            PItemp <- PItemp + GetPIValuesForAABjellvist(aa, C_PI, i, nchar(sequence),sequence)
            PItemp <- PItemp + GetPIValuesForTermBjellvist(ctermAA, C_PI, 1,sequence)
            PItemp <- PItemp + GetPIValuesForTermBjellvist(ntermAA, C_PI,0,sequence)
            OldPI = C_PI
            if (PItemp > 0){
                C_PI = (C_PI + LowPI) / 2.0
                HiPI = OldPI
            }else{
                C_PI = (C_PI + HiPI) / 2.0
                LowPI = OldPI
            }
        }
        }
    return (C_PI)
}


GetPIValuesForAABjellvist <- function(LL, Old_PI, position, last, sequence){
    result <- 0

    if (LL == "C"){
        if (position == 1){
            result <- 1.0 / (1.0 + (10^(8.0 - Old_PI)))
        }else if (position == last){
            result <- 1.0 / (1.0 + (10^(9.0 - Old_PI)))
        }else{
            result <- 1.0 / (1.0 + (10^(8.28 -Old_PI)))
        }
    }

    if (LL == "D"){
        return 1.0 / (1.0 + (10^(3.945 + GetComplexPIValue(LL, position, last,sequence) - Old_PI)));
    }

    if (LL == "E"){
        return 1.0 / (1.0 + (10^(4.38 + GetComplexPIValue(LL, position, last, sequence) - Old_PI)));
    }

    if (LL == "Y"){
        if (position == 1){
            return 1.0 / (1.0 + (10^ (9.84 - Old_PI)));
        }else if (position == last){
            return 1.0 / (1.0 + (10^(10.34 - Old_PI)));
        }else{
            return 1.0 / (1.0 + (10^(9.84 - Old_PI)));
        }
    }

    if (LL == "H"){
        if (position == 1){
            return -1.0 / (1.0 + (10^(Old_PI - 4.96)));
        }else if (position == last){
            return -1.0 / (1.0 + (10^ (Old_PI - 6.89)));
        }else{
            return -1.0 / (1.0 + (10^(Old_PI - pIAdjust_H)));
        }
    }

    if (LL == "K"){
        if (position == 1){
            return -1.0 / (1.0 + (10^ (Old_PI - 10.3)));
        }else if (position == last){
            return -1.0 / (1.0 + (10^(Old_PI - 10.3)));
        }else{
            return -1.0 / (1.0 + (10^(Old_PI - 9.8)));
        }
    }

    if (LL == "R"){
        if (position == 1){
            result <- ((-1.0)/(1.0 + (10^(Old_PI - 10.8))))
        }else if (position == last){
            result <- ((-1.0)/(1.0 + (10^ (Old_PI - 10.8))))
        }else{
            result <- (-1.0/(1.0 + (10^(Old_PI - 12.0))))
        }
    }

    return (result)
}

