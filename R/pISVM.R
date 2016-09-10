#' defaultPeptideTrainData
#' This function acces to the data to trace train the SVM method

defaultPeptideTrainData <- function(){
    filePath <- system.file("extdata", "svmDataDefault.csv", package = "pIR")
    data <- read.table(file = filePath, sep = ",", header = TRUE)
    return(data)
}


#' loadSVMModel
#'
#' This function acces to a SVM model trained
#' @param model The model to be loaded
#'
#'
loadSVMModel <- function(model = "default"){

  if(model=="default"){
    filePath <- system.file("extdata", "svmModelDefault.rda", package = "pIR")
    svm <- load(file = filePath, .GlobalEnv)
    return(svm)
  }
  if(model=="heller"){
    filePath <- system.file("extdata", "svmModelHeller.rda", package = "pIR")
    svm <- load(file = filePath, .GlobalEnv)
    return(svm)
  }
  if(model=="branca"){
    filePath <- system.file("extdata", "svmModelBranca.rda", package = "pIR")
    svm <- load(file = filePath, .GlobalEnv)
    return(svm)
  }
  #otherwise return default trained model
  filePath <- system.file("extdata", "svmModelDefault.rda", package = "pIR")
  svm <- load(file = filePath, .GlobalEnv)
  return(svm)
}



#' transformData
#'
#' This function transform a new data (dataframe) using center and scale transformation.
#' This transformation is necesary previous to use svm method.
#'
#' @param instances The data to be transformed
#' @param defaultTrainingSet A flag to get the training set used in the transformation (if TRUE, default training set will be used)
#' @param trainingSet The training set used to apply transformation.
#'
#' @details If trainingSet flag is not "default", "heller" or "branca", the current data will be used to apply the transformation requeried
#' @details before use a SVM model, the set supply must contains the following variables: calibrated, expasy and aaindex.
#'
transformData <- function(instances = instance, trainingSet = "default"){

  if(trainingSet=="default"){
    #load default dataset
    filePath <- system.file("extdata", "svmDataSetDefault.rda", package = "pIR")
    load(filePath)

  }else if(trainingSet=="heller"){
    #load heller dataset
    filePath <- system.file("extdata", "svmHellerDataSet.rda", package = "pIR")
    load(filePath)

  }else if(trainingSet=="branca"){
    #load branca dataset
    filePath <- system.file("extdata", "svmBrancaDataSet.rda", package = "pIR")
    load(filePath)

  }else if(trainingSet=="current"){
    #otherwise, use the provided dataset to apply transformation...
    data <- instances
  }

   #retrive attributes from dataset. It must be equal to new instance attributes.
   peptides_propeties <- subset(data, select=c("calibrated", "expasy", "aaindex"))

   #create preProcess object
   preObject <- preProcess(peptides_propeties, method = c("center", "scale"))

   #process new instances using (trainingSet) option
   preData <- predict(preObject, newdata = instances)

   return(preData)
}




#' pISVMsequences
#'
#' This function predict the pI from multiple sequences contained into dataframe.
#'
#' @param df The dataset with sequences. It must contains the variables: "calibrated", "expasy" and "aaindex".
#' @param model The SVM-based model to be used in the prediction (use "default", "heller" or "branca" options)
#' @param newModel A flag enabling the posibility to choose a new model to be used.
#'
#' @details By default, this method use a svm-model from the setting of the "model" parameter and keeping
#' @details the parameter "newModel" = FALSE. However, it is possible to build a new svm model from the current dataset
#' @details setting "newModel"=TRUE. To do it, The input dataframe must contains the variables requeried to train
#' @details a svm model: "calibrated", "expasy" and "aaindex".
#'
pISVMsequences <- function(df = dataframe, model = "default", newModel = FALSE){

   if(!newModel){

    #loading default svm model
    svm <- loadSVMModel(model)

    #processing data...requeried previous to use svm model
    processedData <- transformData(instances = df, trainingSet = model)

    #predicting pI with svm model
     pI <- predict(svmModel, newdata = processedData)

    #build dataframe with new pI values predicted
    df$pIsvm <- as.vector(pI, mode = "any")

  }else{

    #build new svm model from current data
    svmModel <- svmPIData(data = df)

    #read varibles from dataframe
    sequences_propeties <- subset(df, select=c("calibrated", "expasy", "aaindex"))

    #processing data...requeried previous to use svm model
    processedData <- transformData(instances = sequences_propeties, trainingSet = "current")

    #predicting pI with svm model
    pI <- predict(svmModel, newdata = processedData)

    #build dataframe with new pI values predicted
    df$pIsvm <- as.vector(pI, mode = "any")

  }

  return(df)
}


#' pISVMpeptide
#'
#' This function predict the pI of a single sequence.
#'
#' @param sequence The sequence to be used
#' @param model The SVM-based model to be used in the prediction (use "default", "heller" or "branca" options)
#'
#'
pISVMpeptide <- function(sequence, model = "default"){
    #computing the features requeried to vectorize the peptide/protein
    aaindex      <- aaIndex(sequence = sequence)
    calibrated   <- pIBjell(sequence = sequence, pkSetMethod = "calibrated")
    expasy       <- pIBjell(sequence = sequence, pkSetMethod = "expasy")

    #loading svm-model
    svm  <- loadSVMModel(model = model)

    #Building sequence as a dataframe
    seq <- data.frame(calibrated, expasy, aaindex)

    #pre-processing new instance using the same model training-dataset processing
    seq <- transformData(instances = seq, trainingSet = model)

    pI <- predict(svmModel, newdata = seq)
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
        load("data/svmDataSetDefault.rda")
        peptides_properties <- subset(data, select=c("calibrated", "expasy", "aaindex"))
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
    peptides_properties <- subset(data, select=c("expasy", "skoog", "calibrated", "solomon", "rodwell", "emboss", "lehninger", "grimsley", "patrickios", "DtaSelect", "aaindex"))
    peptides_experimental <- subset(data, select=c("pIExp"))
    svmProfileValue <- svmProfile(dfExp = peptides_experimental, dfProp = peptides_properties, method = method, numberIter = numberIter)
    return(svmProfileValue)
}


#' svmPIBuildSVMData
#'
#' This function take a data frame in the way of: sequence pIExp. It return a dataframe with features to train a svm model
#' @param originalData The original dra frame
#'
svmPIBuildData <- function(originalData){

    data <- originalData

    colnames(data) <-c("sequence", "pIExp")

    #Add all the bjell methods and pk Sets
    data <- mdply(data, function(sequence, pIExp) { pIBjell(sequence = sequence, pkSetMethod = "calibrated") })

    colnames(data) <-c("sequence", "pIExp", "calibrated")

    #Add all the bjell methods and pk Sets
    data <- mdply(data, function(sequence, pIExp, calibrated) { pIBjell(sequence = sequence, pkSetMethod = "expasy") })

    colnames(data) <-c("sequence", "pIExp", "calibrated", "expasy")

    #Add all the bjell methods and pk Sets
    data <- mdply(data, function(sequence, pIExp, calibrated, expasy) { pIBjell(sequence = sequence, pkSetMethod = "skoog") })

    colnames(data) <-c("sequence", "pIExp", "calibrated", "expasy", "skoog")

    #Add all the bjell methods and pk Sets
    data <- mdply(data, function(sequence, pIExp, calibrated, expasy, skoog) { aaIndex(sequence = sequence) })

    colnames(data) <-c("sequence", "pIExp", "calibrated", "expasy", "skoog", "aaindex")

    #write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

    #print(data)

    return(data)
}

#' svmPIData
#'
#' This function takes a data frame an return a model with the best variables.
#' The training process could be time-consuming.
#'
#' @param data The df
#'

svmPIData <- function(data){
    peptides_propeties <- subset(data, select=c("calibrated", "expasy", "aaindex"))
    peptides_experimental <- subset(data, select=c("pIExp"))
    svmProfileValue <- svmProfile(dfExp = peptides_experimental, dfProp = peptides_propeties, method = "svmRadial", numberIter = 2)
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
    sequence <- toupper(sequence)

    #aaV <- strsplit(sequence, "", fixed = TRUE)
    #aaNTerm <- aaV[[1]][1]
    #aaCTerm <- aaV[[1]][nchar(sequence)]

    aaZimmermanDes <- data.frame(key    <- c("A",  "L",  "R",  "K",   "N",   "M",   "D",  "F", "C",  "P",  "Q",  "S", "E",  "T",  "G",  "W",  "H",  "Y",  "I",  "V"),
                                 value  <- c( 6.00, 5.98, 10.76, 9.74, 5.41, 5.74, 2.77, 5.48, 5.05, 6.30, 5.65, 5.68,3.22, 5.66, 5.97, 5.89, 7.59, 5.66, 6.02, 5.96))
    colnames(aaZimmermanDes) <- c("key", "value")
    #count <- 2;

    lev <- c("A",  "L",  "R", "K", "N", "M", "D", "F", "C","P", "Q", "S", "E", "T", "G", "W", "H", "Y", "I", "V")

    aaTable <- table(factor(prot <- strsplit(sequence, "")[[1]], levels = lev))

    pKNValues <- loadNTermPK(pkSet = "calibrated")
    pKCValues <- loadCTermPK(pkSet = "calibrated")

    pKNTerm <- retrievePKValue(prot[1], pKNValues)
    pKCTerm <- retrievePKValue(prot[length(prot)], pKCValues)

    oldw <- getOption("warn")
    options(warn = -1)

    temp <- aaTable*aaZimmermanDes

    options(warn = oldw)
    zimm <- pKNTerm + pKCTerm + sum(temp$value)

    zimm = zimm/(length(prot)+2)

    return (as.numeric(zimm))
}

