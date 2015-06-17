#' svmPI
#'
#' This function compute the isoelectric point usig support vector machines
#'
#' @param df data frame where the first column is the peptde sequence and the second column the experimental
#'
svmPI <- function(seq){

}

# dafaultTrainData
# This fucntion acces to the data to trace train the SVM method

defaultTrainData <- function(){
    filePath <- system.file("extdata", "svmDataShort.csv", package = "pIR")
    data <- read.table(file = filePath, sep = ",", header = TRUE)
    return(data)
}

# svmPIDefault
#
# This function train the original dataset from the SVM dataset a get the model

svmPIDefault <- function(){
    data <- defaultTrainData()
    data <- mdply(data, function(peptide, pIExp) { pIBjell(sequence = peptide, pkSetMethod = "bjell") })
    colnames(data) <-c("peptide", "pIExp", "bjell")
    data <- mdply(data, function(peptide, pIExp, bjell){

        sequence <- toupper(peptide)
        aaV <- strsplit(sequence, "", fixed = TRUE)
        aaNTerm <- aaV[[1]][1]
        aaCTerm <- aaV[[1]][nchar(sequence)]

        pKNValues <- loadNTermPK(pkSet = "bjell")
        pKCValues <- loadCTermPK(pkSet = "bjell")

        pKNTerm <- retrievePKValue(aaNTerm, pKNValues)
        pKCTerm <- retrievePKValue(aaCTerm, pKCValues)

        zimmerman <- pKNTerm + pKCTerm

        aaZimmermanDes <- data.frame(key <-c("A",  "L",  "R", "K", "N", "M", "D", "F", "C","P", "Q", "S", "E", "T", "G", "W", "H", "Y", "I", "V"),
                                     value  <- c( 6.00, 5.98, 10.76, 9.74, 5.41, 5.74, 2.77, 5.48, 5.05, 6.30, 5.65, 5.68,3.22, 5.66,5.97, 5.89, 7.59, 5.66, 6.02, 5.96))
        colnames(aaZimmermanDes) <- c("key", "value")
        count <- 2;

        for(i in 1:nchar(sequence)) {

            aa <- aaV[[1]][i]
            if(retrievePKValue(aa, aaZimmermanDes)){
                zimmerman <-  zimmerman + retrievePKValue(aa, aaZimmermanDes)
                count <- count + 1
            }
        }

        zimmerman = zimmerman/count
        return (zimmerman)
    })
    colnames(data) <-c("peptide", "pIExp", "bjell", "aaIndex")
    print(data)
    return(data)
}

#'aaIndex
#'
#' This function compute the isoelectric point using the aa Index
#' @param sequence The peptide sequence
#'
aaIndex <- function(sequence){

    sequence <- toupper(sequence)
    aaV <- strsplit(sequence, "", fixed = TRUE)
    aaNTerm <- aaV[[1]][1]
    aaCTerm <- aaV[[1]][nchar(sequence)]

    pKNValues <- loadNTermPK(pkSet = "bjell")
    pKCValues <- loadCTermPK(pkSet = "bjell")

    pKNTerm <- retrievePKValue(aaNTerm, pKNValues)
    pKCTerm <- retrievePKValue(aaCTerm, pKCValues)

    zimmerman <- pKNTerm + pKCTerm

    aaZimmermanDes <- data.frame(key <-c("A",  "L",  "R", "K", "N", "M", "D", "F", "C","P", "Q", "S", "E", "T", "G", "W", "H", "Y", "I", "V"),
    value  <- c( 6.00, 5.98, 10.76, 9.74, 5.41, 5.74, 2.77, 5.48, 5.05, 6.30, 5.65, 5.68,3.22, 5.66,5.97, 5.89, 7.59, 5.66, 6.02, 5.96))
    colnames(aaZimmermanDes) <- c("key", "value")
    count <- 2;

    for(i in 1:nchar(sequence)) {

        aa <- aaV[[1]][i]
        if(retrievePKValue(aa, aaZimmermanDes)){
            zimmerman <-  zimmerman + retrievePKValue(aa, aaZimmermanDes)
            count <- count + 1
    }
    }

    zimmerman = zimmerman/count
    return (zimmerman)
}


