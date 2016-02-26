#' loadPkSetIterative
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
#'   \item toseland   : Toseland, C.P., et al. (2006). "PPD v1.0--an integrated, web-accessible database of experimentally determined protein pKa values." Nucleic Acids Res 34(Database issue): D199-203.
#'   \item thurlkill  : Thurlkill, R.L., et al. (2006). "pK values of the ionizable groups of proteins." Protein Sci 15(5): 1214-1218.
#'   \item nozaki_tanford  : Nozaki, Y., and C. Tanford. 1967. J. Am. Chem. Soc. 89:736–742.  
#' }
#'
#' @param pkSetMethod name of the pk
#'

loadPkSetIterative <- function(pkSetMethod = "solomon"){
    pkValues <- c()
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
    if(pkSetMethod == "toseland"){
        pkValues <- data.frame(key=c("NTerm","CTerm","C","D","E","H","K","R","Y"), c(8.71,3.19,6.87,3.6,4.29,6.33,10.45,12.0,9.61))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSetMethod == "thurlkill"){
        pkValues <- data.frame(key=c("NTerm","CTerm","C","D","E","H","K","R","Y"), c(8.0,3.67,8.55,3.67,4.25,6.54,10.4,12.0,9.84))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSetMethod == "nozaki_tanford"){
        pkValues <- data.frame(key=c("NTerm","CTerm","C","D","E","H","K","R","Y"), c(7.50,3.67,9.5,4.0,4.25,6.3,10.4,12.0,9.6))
        colnames(pkValues) <- c("key", "value")
    }

    return (pkValues)
}

#' loadNTermPK
#'
#' This function set NTerm contribution fo the Bjell method
#'
#' @param pkSet The pk Set to be use
#'

loadNTermPK <- function(pkSet = "expasy"){
    pkValues <- c()
    if(pkSet == "expasy"){
        pkValues <- data.frame(key=c("A",  "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P","S",  "T",  "W", "Y", "V"),c( 7.59, 7.5, 7.5, 7.5, 7.5, 7.7, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.0, 7.5, 8.359999999999999, 6.93, 6.82, 7.5, 7.5, 7.44))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSet == "skoog"){
        pkValues <- data.frame(key=c("A",  "R",    "N",  "D",  "C",     "E",  "Q",  "G", "H", "I",  "L", "K",    "M",  "F",  "P",  "S",  "T",   "W",  "Y",    "V"),c( 9.69, 9.0399, 8.80, 9.82, 10.7799, 9.76, 9.13, 9.6, 9.17,9.68, 9.6, 8.9499, 9.21, 9.13, 10.6, 9.15, 10.43, 9.39, 9.1099, 9.6199))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSet == "bjell"){
        pkValues <- data.frame(key=c("A", "R",  "N",  "D", "C",    "E",  "Q",  "G", "H",  "I",  "L",  "K",  "M",  "F",  "P",    "S",  "T",  "W",  "Y",  "V"),c( 7.5, 6.76, 7.22, 7.7, 8.1199, 7.19, 6.73, 7.5, 7.18, 7.48, 7.46, 6.67, 6.98, 6.96, 8.3599, 6.86, 7.02, 7.11, 6.83, 7.44))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSet == "calibrated"){
        pkValues <- data.frame(key=c("A", "R",  "N",  "D", "C",    "E",  "Q",  "G", "H",  "I",  "L",  "K",  "M",  "F",  "P",    "S",  "T",  "W",  "Y",  "V"),c(7.59, 7.5,  6.7,  7.5, 6.5,    7.7,  7.5,  7.5, 7.5,  7.5,  7.5,  7.5,  7.0,  7.5,  8.3599, 6.93, 6.82, 7.5,  7.5,  7.44))
        colnames(pkValues) <- c("key", "value")
    }
    return(pkValues)
}

#' loadCTermPK
#'
#' This function set CTerm contribution fo the Bjell method
#'
#' @param pkSet The pk Set to be use
#'
#'
loadCTermPK <- function(pkSet = "expasy"){
    pkValues <- c()
    if(pkSet == "expasy"){
        pkValues <- data.frame(key=c("A",  "R",  "N",  "D",  "C",  "E",  "Q",  "G",  "H",  "I",  "L",  "K",  "M",  "F",  "P",  "S",  "T",  "W",  "Y",  "V"),
                                   c( 3.55, 3.55, 3.55, 4.55, 3.55, 4.75, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55 ,3.55, 3.55, 3.55))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSet == "skoog"){
        pkValues <- data.frame(key=c("A",  "R",    "N",  "D",  "C",  "E",  "Q",  "G",  "H",  "I",  "L",  "K",  "M",  "F",  "P",  "S",  "T",  "W",  "Y", "V"),
                                   c(2.35,  2.17,   2.02, 2.09, 1.71, 2.19, 2.17, 2.34, 1.82, 2.36, 2.36, 2.18, 2.28, 1.83, 1.99, 2.21, 2.63, 2.38, 2.2, 2.32))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSet == "bjell"){
        pkValues <- data.frame(key=c("A",  "R",  "N",  "D",  "C",  "E",  "Q",  "G",  "H",  "I",  "L",  "K",  "M",  "F",  "P",    "S",  "T",  "W",  "Y",  "V"),
                                   c( 2.35, 2.17, 2.02, 2.09, 1.71, 2.19, 2.17, 2.34, 1.82, 2.36, 2.36, 2.18, 2.28, 1.83, 1.99,   2.21, 2.63, 2.38, 2.2,  2.32))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSet == "calibrated"){
        pkValues <- data.frame(key=c("A",  "R",  "N",  "D",  "C",  "E",  "Q",  "G",  "H",  "I",  "L",  "K",  "M",  "F",  "P",  "S",  "T",  "W",  "Y",  "V"),
                                   c( 3.55, 3.55, 3.55, 4.55, 3.55, 4.75, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55, 3.55))
        colnames(pkValues) <- c("key", "value")
    }
    return(pkValues)

}

#' setGroupPK
#'
#' This function set the group contribution fo the Bjell method
#'
#' @param pkSetMethod The pk Set to be use
#'

loadGroupPK <- function(pkSet = "expasy"){
    pkValues <- c()
    if(pkSet == "expasy"){
        pkValues <- data.frame(key=c("R", "D", "C", "E", "H", "K", "Y"),
                                   c(-12.0, 4.05, 9.0, 4.45, -5.98, -10.0, 10.0))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSet == "skoog"){
        pkValues <- data.frame(key=c("R",   "D",  "C",  "E",   "H",   "K",    "Y"),
                                   c(-12.48, 3.86, 8.33, 4.25, -6.0, -10.5299, 10.07))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSet == "bjell"){
        pkValues <- data.frame(key=c("R",  "D",  "C",    "E",   "H",   "K",  "Y"),
                                   c(-12.5, 4.07, 8.2799, 4.45, -6.08, -9.80, 9.84))
        colnames(pkValues) <- c("key", "value")
    }
    if(pkSet == "calibrated"){
        pkValues <- data.frame(key=c("R", "D", "C", "E", "H", "K", "Y"),
                                   c(-12.0, 4.05, 9.0, 4.45, -5.98, -10.0, 10.0))
        colnames(pkValues) <- c("key", "value")
    }
    return(pkValues)

}

#' retrievePKValue
#'
#' This function retrieve the pK values for an specific aminoacid
#'
#' @param aminoacid
#' @param pKIterative

retrievePKValue <- function(aa, pKSet){
    pkValue <- NA
    for(i in 1:nrow(pKSet)){
        if(pKSet[i,"key"] == aa){
            pkValue <- pKSet[i,"value"]
        }
    }
    return (pkValue)
}


pICofactorAdjust_E_Nterm <- function(){
    pKValue <- data.frame(key1 <- c(0.058000, -0.100000, 0.000000, 0.012000, 0.004000, 0.068000, 0.002000, 0.002000, 0.066000, 0.000000, -0.064000, 0.044000, 0.042000, 0.068000, 0.000000, 0.036000, 0.050000, -0.324000, 0.048000, 0.034000, 0.000000, 0.050000, 0.084000, 0.000000, 0.078000, -0.414000),
                          key2 <- c(-0.002000, 5.504000, -0.112000, 0.000000, 0.076000, -0.006000, -0.038000, -0.108000, 0.020000, 0.000000, -0.072000, 0.020000, 0.022000, 0.008000, 0.000000, -0.010000, 0.004000, -0.084000, -0.032000, -0.020000, 0.000000, 0.004000, -0.036000, 0.000000, -0.038000, -0.230000),
                          key3 <- c(0.000000, 5.448000, 0.054000, 0.040000, 0.012000, 0.000000, -0.010000, -0.060000, -0.020000, 0.000000, -0.118000, -0.004000, -0.040000, -0.018000, 0.000000, 0.010000, 0.000000, -0.128000, -0.020000, -0.008000, 0.000000, 0.000000, 0.000000, 0.000000, -0.018000, -0.158000))

    colnames(pKValue) <- c(1,2,3)

    return (pKValue)
}

pICofactorAdjust_E_Cterm <- function(){
    pKValue <- data.frame(key2 <- c(0.022000, 0.014000, 0.036000, -0.020000, 0.120000, 0.006000, 0.000000, -0.092000, -0.002000, 0.000000, -0.038000, 0.002000, -0.008000, 0.006000, 0.000000, 0.072000, 0.000000, -0.022000, 0.020000, 0.022000, 0.000000, 0.012000, 0.006000, 0.000000, 0.050000, 5.496000),
                          key2 <- c(0.000000, 0.012000, -0.050000, -0.040000, -0.028000, 0.020000, 0.012000, -0.132000, 0.008000, 0.000000, -0.086000, 0.008000, 0.012000, 0.022000, 0.000000, 0.006000, 0.008000, -0.088000, -0.002000, -0.012000, 0.000000, 0.000000, 0.006000, 0.000000, 0.022000, 5.024000),
                          key3 <- c(0.000000, 0.000000, 0.016000, 0.000000, 0.026000, 0.000000, 0.000000, -0.070000, 0.000000, 0.000000, -0.032000, 0.000000, -0.008000, 0.004000, 0.000000, 0.004000, -0.010000, -0.028000, 0.000000, 0.006000, 0.000000, 0.000000, 0.032000, 0.000000, -0.022000, 4.400000))

    colnames(pKValue) <- c(1,2,3)
    return (pKValue)
}


pICofactorAdjust_D_Nterm <- function(){
    pKValue <- data.frame(key1 <- c(0.012000, -0.108000, 0.004000, 0.132000, 0.122000, 0.040000, -0.016000, -0.056000, 0.078000, 0.000000, -0.090000, 0.040000, 0.028000, 0.044000, 0.000000, 0.008000, 0.020000, -0.144000, 0.006000, 0.000000, 0.000000, 0.034000, 0.058000, 0.000000, 0.054000, -0.722000),
                          key2 <- c(-0.012000, 5.312000, -0.084000, 0.034000, 0.100000, -0.008000, -0.040000, -0.084000, 0.002000, 0.000000, -0.104000, 0.000000, 0.006000, -0.022000, 0.000000, -0.014000, -0.002000, -0.038000, -0.042000, -0.020000, 0.000000, 0.000000, -0.022000, 0.000000, -0.010000, -0.252000),
                          key3 <- c(0.000000, 4.544000, 0.016000, 0.070000, 0.040000, 0.008000, -0.020000, -0.082000, 0.000000, 0.000000, -0.084000, -0.012000, -0.024000, -0.012000, 0.000000, 0.010000, 0.000000, -0.136000, -0.012000, -0.020000, 0.000000, 0.000000, -0.002000, 0.000000, -0.008000, -0.140000))

    colnames(pKValue) <- c(1,2,3)
    return (pKValue)
}

pICofactorAdjust_D_Cterm <- function(){
    pKValue <- data.frame(key1 <- c(-0.014000, 0.558000, 0.008000, 0.028000, 0.108000, 0.012000, -0.020000, -0.136000, 0.012000, 0.000000, -0.058000, 0.020000, -0.010000, -0.006000, 0.000000, -0.044000, -0.012000, -0.058000, -0.008000, 0.008000, 0.000000, 0.014000, 0.004000, 0.000000, 0.014000, 5.272000),
                          key2<- c(0.000000, 0.052000, 0.002000, 0.088000, 0.060000, 0.010000, 0.000000, -0.136000, -0.002000, 0.000000, -0.012000, 0.000000, 0.006000, -0.050000, 0.000000, 0.028000, 0.026000, -0.054000, -0.010000, -0.062000, 0.000000, 0.020000, 0.046000, 0.000000, 0.018000, 6.312000),
                          key3 <- c(0.000000, 0.000000, 0.014000, -0.020000, 0.000000, 0.000000, -0.020000, -0.024000, 0.020000, 0.000000, -0.096000, 0.000000, 0.000000, 0.020000, 0.000000, 0.000000, -0.012000, -0.104000, -0.058000, 0.000000, 0.000000, 0.000000, -0.010000, 0.000000, 0.000000, 4.368000))

    colnames(pKValue) <- c(1,2,3)
    return (pKValue)

};

pICofactorAdjust_Cterm <- function(){
    pKValue <- data.frame(key2 <- c(0.290000, 5.704000, 5.296000, 0.028000, 0.310000, 0.536000, 5.896000, -0.290000, 0.212000, 0.000000, 0.000000, 0.164000, 0.196000, 0.076000, 0.000000, 0.632000, 0.060000, -0.080000, 0.026000, 0.108000, 0.000000, 0.338000, 0.088000, 0.000000, 0.192000, 5.856000),
                          key2 <- c(-0.040000, 4.000000, -0.112000, 0.120000, 0.156000, 0.000000, -0.066000, 0.012000, 0.030000, 0.000000, -0.158000, -0.022000, -0.056000, -0.074000, 0.000000, -0.016000, -0.020000, -0.184000, -0.098000, -0.066000, 0.000000, 0.020000, 0.008000, 0.000000, -0.008000, 5.800000),
                          key3 <- c(-0.020000, 3.944000, -0.076000, 0.076000, 0.096000, -0.004000, 0.010000, 0.180000, -0.026000, 0.000000, -0.044000, -0.008000, -0.012000, 0.016000, 0.000000, 0.002000, -0.008000, -0.276000, -0.008000, -0.048000, 0.000000, -0.020000, 0.002000, 0.000000, 0.010000, 4.488000))
    colnames(pKValue) <- c(1,2,3)
    return (pKValue)
}


