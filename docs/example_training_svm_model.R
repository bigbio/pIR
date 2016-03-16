#This script train a new svm model that can be used to predict the isoelectric point from new instances.

#reading peptide/protein dataset from any source. It must contain two columns: sequences and experimental pI.

data <- read.table(file = "data/svmDataDefault.csv", header = FALSE, sep = ",")

colnames(data) <-c("sequence", "pIExp")

#Add the bjell isoelectric point using calibrated pk Set

data <- mdply(data, function(sequence, pIExp) { pIBjell(sequence = sequence, pkSetMethod = "calibrated") })

colnames(data) <-c ("sequence", "pIExp", "calibrated")

#Add the bjell isoelectric point using expasy pK Set

data <- mdply(data, function(sequence, pIExp, calibrated) { pIBjell(sequence = sequence, pkSetMethod = "expasy") })

colnames(data) <-c("sequence", "pIExp", "calibrated", "expasy")

#Add the aaindex property

data <- mdply(data, function(sequence, pIExp, calibrated, expasy) { aaIndex(sequence = sequence) })

colnames(data) <-c("sequence", "pIExp", "calibrated", "expasy", "aaindex")

#write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

#Save data with attributes (pI values and aaindex property)

save(data, file = "data/svmDataSetDefault.rda")

#retrieve attributes subset (predictors)

peptides_properties <- subset(data, select=c("calibrated", "expasy","aaindex"))

#getting class variable (pI experimental)

peptides_experimental <- subset(data, select=c("pIExp"))

#training svm model using the svmProfile function

svmModel <- svmProfile(dfExp = peptides_experimental, dfProp = peptides_properties, method = "svmRadial", numberIter = 2)

#Save the new model trained
save(svmModel, file="data/svmModelDefault.rda")
