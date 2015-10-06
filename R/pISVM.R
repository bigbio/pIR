# dafaultTrainData
# This function acces to the data to trace train the SVM method

defaultPeptideTrainData <- function(){
    filePath <- system.file("extdata", "svmPeptideData.csv", package = "pIR")
    data <- read.table(file = filePath, sep = ",", header = TRUE)
    return(data)
}

# dafaultTrainData
# This function acces to the data to trace train the SVM method

defaultProteinTrainData <- function(){
    filePath <- system.file("extdata", "svmProteinData.csv", package = "pIR")
    data <- read.table(file = filePath, sep = ",", header = TRUE)
    return(data)
}

# svmPIDefault
#
# This function train the original dataset from the SVM dataset a get the model

svmBuildPeptideData <- function(loadData = FALSE, method = "svmRadial", numberIter = 2){
    data <- defaultPeptideTrainData()
    if(loadData){
        data <- svmPIBuildSVM(originalData = data)
        save(data,file="data/svmPeptideData.rda")
    }else{
        load("data/svmPeptideData.rda")
    }

    peptides_properties <- subset(data, select=c("bjell", "expasy", "skoog", "calibrated", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios", "DtaSelect", "aaindex"))
    peptides_experimental <- subset(data, select=c("pIExp"))
    svmProfileValue <- svmProfile(dfExp = peptides_experimental, dfProp = peptides_properties, method = method, numberIter = numberIter)
    return(svmProfileValue)
}

svmBuildProteinData <- function(loadData = FALSE, method = "svmRadial", numberIter = 2){
    data <- defaultProteinTrainData()
    if(loadData){
        data <- svmPIBuildSVM(originalData = data)
        save(data,file="data/svmProteinData.rda")
    }else{
        load("data/svmProteinData.rda")
    }
    peptides_properties <- subset(data, select=c("bjell", "expasy", "skoog", "calibrated", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios", "DtaSelect", "aaindex"))
    peptides_experimental <- subset(data, select=c("pIExp"))
    svmProfileValue <- svmProfile(dfExp = peptides_experimental, dfProp = peptides_properties, method = method, numberIter = numberIter)
    return(svmProfileValue)
}


#' svmPIBuildSVMData
#'
#' This function take a data frame in the way of: sequence pIExp
#' @param originalData The original dra frame
#'
svmPIBuildSVM <- function(originalData){
    data <- originalData

    #Add all the bjell methods and pk Sets
    data <- mdply(data, function(sequence, pIExp) { pIBjell(sequence = sequence, pkSetMethod = "bjell") })
    colnames(data) <-c("sequence", "pIExp", "bjell")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    data <- mdply(data, function(sequence, pIExp, bjell) { pIBjell(sequence = sequence, pkSetMethod = "expasy") })
    colnames(data) <-c("sequence", "pIExp","bjell", "expasy")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    data <- mdply(data, function(sequence, pIExp, bjell, expasy) { pIBjell(sequence = sequence, pkSetMethod = "skoog") })
    colnames(data) <-c("sequence", "pIExp", "bjell","expasy", "skoog")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog) { pIBjell(sequence = sequence, pkSetMethod = "calibrated") })
    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "calibrated")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    # Add iterative values
    data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, calibrated) { pIIterative(sequence = sequence, pkSetMethod = "solomon") })
    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "calibrated", "solomon")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, calibrated, solomon) { pIIterative(sequence = sequence, pkSetMethod = "rodwell") })
    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "calibrated", "solomon", "rodwell")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, calibrated, solomon, rodwell) { pIIterative(sequence = sequence, pkSetMethod = "emboss") })
    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "calibrated", "solomon", "rodwell", "emboss")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, calibrated, solomon, rodwell, emboss) { pIIterative(sequence = sequence, pkSetMethod = "lehninger") })
    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "calibrated", "solomon", "rodwell", "emboss", "lehninger")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, calibrated, solomon, rodwell, emboss, lehninger) { pIIterative(sequence = sequence, pkSetMethod = "grimsley") })
    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "calibrated", "solomon", "rodwell", "emboss", "lehninger", "grimsley")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, calibrated, solomon, rodwell, emboss, lehninger, grimsley) { pIIterative(sequence = sequence, pkSetMethod = "patrickios") })
    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "calibrated", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog, calibrated, solomon, rodwell, emboss, lehninger, grimsley, patrickios) { pIIterative(sequence = sequence, pkSetMethod = "DtaSelect") })
    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "calibrated", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios", "DtaSelect")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    data <- read.table(file="data.csv", header = TRUE, sep = ",")
    data <- mdply(data, function(X, sequence, pIExp, bjell, expasy, skoog, calibrated, solomon, rodwell, emboss, lehninger, grimsley, patrickios, DtaSelect){
        sequence <- toupper(sequence)
        sequence <- reformat(seq = sequence)
        sequence <- removePTM(seq = sequence) #to avoid PTMs character in sequence.

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
        print(zimmerman)
        return (zimmerman)
    })
    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "calibrated", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios", "DtaSelect", "aaindex")
    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    return(data)
}

#' svmPIData
#'
#' This function takes a data frame an return a model with the best variables
#' @param data The df
#'

svmPIData <- function(data){
    peptides_propeties <- subset(data, select=c("bjell", "aaIndex"))
    peptides_experimental <- subset(data, select=c("pIExp"))
    svmProfileValue <- svmProfile(dfExp = peptides_experimental, dfProp = peptides_propeties, method = "svmRadial", numberIter = 5)
    return(svmProfileValue)
}

svmProfile <- function(dfExp, dfProp, method = "svmRadial", numberIter = 2){

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
    svmProfileValue <- rfe(x=trainDescr, y = trainClass, sizes = c(1:ncol(peptides_desc)),rfeControl = rfeControl(functions = caretFuncs,number = numberIter),method = method,fit = FALSE);

    return (svmProfileValue)
}

#'aaIndex
#'
#' This function compute the isoelectric point using the aa Index
#' @param sequence The peptide sequence
#'
aaIndex <- function(sequence){

    sequence <- reformat(seq = sequence)
    sequence <- removePTM(seq = sequence) #to avoid PTMs character in sequence.

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
        if(is.na(retrievePKValue(aa, aaZimmermanDes))){
            zimmerman <-  zimmerman + retrievePKValue(aa, aaZimmermanDes)
            count <- count + 1
        }
    }

    zimmerman = zimmerman/count
    return (zimmerman)
}
