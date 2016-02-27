#This script train a new svm model that can be used to predict the isoelectric point from new instances.

#reading peptide/protein dataset from any source. It must contain two columns: sequences and experimental pI.

data <- read.table(file = "pIR/data/svmDataDefault.csv", header = FALSE, sep = ",")

colnames(data) <-c("sequence", "pIExp")

#Add the bjell isoelectric point using default bjell pk Set

data <- mdply(data, function(sequence, pIExp) { pIBjell(sequence = sequence, pkSetMethod = "bjell") })

colnames(data) <-c ("sequence", "pIExp", "bjell")

#Add the bjell isoelectric point using expasy pK Set

data <- mdply(data, function(sequence, pIExp, bjell) { pIBjell(sequence = sequence, pkSetMethod = "expasy") })

colnames(data) <-c("sequence", "pIExp", "bjell", "expasy")

#Add the aaindex property

data <- mdply(data, function(sequence, pIExp, bjell, expasy) { aaIndex(sequence = sequence) })

colnames(data) <-c("sequence", "pIExp", "bjell", "expasy", "aaindex")

#write.table(data, file = "data.csv", sep = ",", col.names = NA, qmethod = "double")

#Save data with attributes (pI values and aaindex property)

save(data, "svmDataDefault.rda")

#retrieve attributes subset (predictors)

peptides_properties <- subset(data, select=c("bjell", "expasy","aaindex"))

#getting class variable (pI experimental)

peptides_experimental <- subset(data, select=c("pIExp"))

#training svm model using the svmProfile function

svmModel <- svmProfile(dfExp = peptides_experimental, dfProp = peptides_properties, method = "svmRadial", numberIter = 2)

#Save the new model trained
save(svmModel, "svmModel.rda")
