# dafaultTrainData
# This function acces to the data to trace train the SVM method

defaultPeptideTrainData <- function(){
    filePath <- system.file("extdata", "svmDataDefault.csv", package = "pIR")
    data <- read.table(file = filePath, sep = ",", header = TRUE)
    return(data)
}

loadDefaultModel <- function(){
    svmModel <- load("data/svmModelDefault.rda")
    return(svmModel)
}


# This function transform a new data (dataframe) using center and scale transformation.
# This transformation is necesary previous to use svm method.

transformData <- function(instance = instance){

    #load default dataset
    svmDataDefault <- load("data/svmDataDefault.rda")

    #retrive attributes from default dataset. It must be equal to new instance attribute.
    peptides_propeties <- subset(data, select=c("bjell", "expasy", "aaindex"))

    #create preProcess object
    preObject <- preProcess(peptides_propeties, method = c("center", "scale"))

    #process new instance using default(training) setting
    preData <- predict(preObject, newdata = instance)

    return(preData)
}


# This function predict the pI from any sequence using SVM model.
# Use as > pI <- pISVMpeptide("NSENATDANQVAHMYQSQDQVK")
# value  > pI = 3.766294

pISVMpeptide <- function(sequence){

    aaindex    <- aaIndex(sequence = sequence)
    bjell      <- pIBjell(sequence = sequence, pkSetMethod = "bjell")
    expasy     <- pIBjell(sequence = sequence, pkSetMethod = "expasy")

    svmModel   <- loadDefaultModel()

    #Building sequence as dataframe
    seq <- data.frame(bjell, expasy, aaindex)

    #pre-processing new instance using same training data processing
    seq <- transformData(seq)

    pI <- predict(svmModelDefault, newdata = seq)
    return(pI)
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

svmBuildPeptideData <- function(defaultModel = FALSE, method = "svmRadial", numberIter = 2){
    data <- defaultPeptideTrainData()
    if(defaultModel){
        loadDefaultModel()
        svmModel <- svmModel
    }else{
        load("data/svmDataDefault.rda")
        peptides_properties <- subset(data, select=c("bjell", "expasy", "aaindex"))
        peptides_experimental <- subset(data, select=c("pIExp"))
        svmModel <- svmProfile(dfExp = peptides_experimental, dfProp = peptides_properties, method = method, numberIter = numberIter)

    }
    return(svmModel)
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

    #Add all the bjell methods and pk Sets
    data <- mdply(data, function(sequence, pIExp, bjell) { pIBjell(sequence = sequence, pkSetMethod = "expasy") })

    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy")

    #Add all the bjell methods and pk Sets
    data <- mdply(data, function(sequence, pIExp, bjell, expasy) { pIBjell(sequence = sequence, pkSetMethod = "skoog") })

    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog")

    #Add all the bjell methods and pk Sets
    data <- mdply(data, function(sequence, pIExp, bjell, expasy, skoog) { aaIndex(sequence = sequence) })

    colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "skoog", "aaindex")

    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    #print(data)

    return(data)
}

#' svmPIData
#'
#' This function takes a data frame an return a model with the best variables
#' @param data The df
#'

svmPIData <- function(data){
    peptides_propeties <- subset(data, select=c("bjell", "expasy", "aaIndex"))
    peptides_experimental <- subset(data, select=c("pIExp"))
    svmProfileValue <- svmProfile(dfExp = peptides_experimental, dfProp = peptides_propeties, method = "svmRadial", numberIter = 300)
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
    inTrain <- createDataPartition(peptides_class, p = 3/4, list = FALSE)[,1];

    #Create the Training Dataset for Descriptors
    trainDescr <- peptides_desc[inTrain,];

    # Create the Testing dataset for Descriptors
    testDescr <- peptides_desc[-inTrain,];

    trainClass <- peptides_class[inTrain];
    testClass <- peptides_class[-inTrain];

    mod <- getModelInfo("svmRadial", regex = FALSE)[[1]]

    mod$predict <- function(modelFit, newdata, submodels = NULL) {
        svmPred <- function(obj, x) {
            hasPM <- !is.null(unlist(obj@prob.model))
            if(hasPM) {
                pred <- lev(obj)[apply(predict(obj, x, type = "probabilities"), 1, which.max)]
            } else pred <- predict(obj, x)
            pred
        }
        out <- try(svmPred(modelFit, newdata), silent = TRUE)
        if(is.character(lev(modelFit))) {
            if(class(out)[1] == "try-error") {
                warning("kernlab class prediction calculations failed; returning NAs")
                out <- rep("", nrow(newdata))
                out[seq(along = out)] <- NA
            }
        } else {
            if(class(out)[1] == "try-error") {
                warning("kernlab prediction calculations failed; returning NAs")
                out <- rep(NA, nrow(newdata))
            }
        }
        if(is.matrix(out)) out <- out[,1]
        out
    }

    #Support Vector Machine Object
    svmProfileValue <- rfe(trainDescr,trainClass, sizes = (1:3),rfeControl = rfeControl(functions = caretFuncs,number = numberIter, verbose = TRUE),method = mod);

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
#   print(zimmerman)
    return (zimmerman)
}

