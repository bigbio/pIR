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
    filePath <- system.file("extdata", "svmData.csv", package = "pIR")
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

    peptides_propeties <- subset(data, select=c("bjell", "aaIndex"))
    peptides_experimental <- subset(data, select=c("pIExp"))
    svmProfileValue <- svmProfile(dfExp = peptides_experimental, dfProp = peptides_propeties)

    return(svmProfile)
}

svmProfile <- function(dfExp, dfProp){

    #load Data
    # This is the data file with the descriptors:
    peptides_desc <- as.matrix(dfProp);
    # This is the Data File with the Experimental Isoelectric Point
    peptides_class<-as.matrix(dfExp);

    #Scale and center data
    peptides_desc<- scale(peptides_desc,center=TRUE,scale=TRUE);

    #Divide the dataset in train and test sets

    # Create an index of the number to train
    inTrain <- createDataPartition(peptides_class, p = 3/4, list = FALSE);

    #Create the Training Dataset for Descriptors
    trainDescr <- peptides_desc[inTrain,];

    # Create the Testing dataset for Descriptors
    testDescr <- peptides_desc[-inTrain,];

    trainClass <- peptides_class[inTrain];
    testClass <- peptides_class[-inTrain];

    #Support Vector Machine Object
    svmProfileValue <- rfe(x=trainDescr, y = trainClass, sizes = c(1:5),rfeControl = rfeControl(functions = caretFuncs,number = 2),method = "svmRadial",fit = FALSE);

    return (svmProfileValue)
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


