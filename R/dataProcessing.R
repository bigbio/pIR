#' processData
#'
#' This function remove NULL and empty values from the data
#' @param data frame or table

processData <- function (data)
{
    row.has.na <- apply(data, 1, function(x){any(is.na(x))})
    sum(row.has.na)
    return (data[!row.has.na,])
}

#' removeFirstColumn
#'
#' This function remove the first in the table or frame
#'
#' @param data frame or table

removeFirstColumn <- function(data)
{
    return (data[-1])
}

#' computeCorrelationMatrix
#'
#' This function compute the correlation for the frame
#' @param data frame
#' @param method

computeCorrelationMatrix <- function(data, method = "pearson"){
    matrixData <- as.matrix(data)
    class(matrixData)
    return (rcorr(matrixData, type = method))
}

#' getCorrelationExperimental
#'
#' This furntion return a matrix where each colum is a correlation respected the experimental value
#' @param data
#' @param method

computeCorrelationExperimental <- function(data, method = "pearson"){
    matrixValues <- computeCorrelationMatrix(data, method)
    values <- matrixValues$r[1,]
    values <- values[-1]
    return (tFrame(data.frame(values)))
}

#' rmse
#'
#' This function return the rmse value between to vectors
#' @param expValue the value expected
#' @param predictValue the value predicted

rmse <- function(expValue, predictValue){
    RMSE <- sqrt(mean((expValue-predictValue)^2))
    return (RMSE)
}

#' computeRMSEExperimental
#'
#' This function plot the RMSE experimental for all the predicted values
#' @param dat

computeRMSEExperimental <- function(dat){
    nm <- names(dat)
    rmseValue <- data.frame(matrix(ncol = ncol(dat), nrow = 1))
    colnames(rmseValue) <- nm
    rmseValue <- removeFirstColumn(rmseValue)
    for (i in nm) {
        if(dat[1L] != dat[i]){
            newData <- data.frame(x = dat[1L], y = dat[i])
            colnames(newData) <- c("x", "y")
            value <- rmse(newData$x, newData$y)
            rmseValue[i] <- value
        }
    }
    return (rmseValue)
}

#' bindRMSECorrelation
#'
#' This function bind the RMSE and correlation in one frame
#'
#' @param correlation
#' @param RMSE

bindRMSECorrelation <- function(rmse, corr){
    total <- rbind(corr, rmse)
    row.names(total) <- c("correlation", "RMSE")
    return (total)
}

#' bindRMSECorrelationFrame
#'
#' This function bind the RMSE and correlation in one frame
#'
#' @param data
#' @param method

bindRMSECorrelationFrame <- function(dat, method = "pearson"){
    corr <- computeCorrelationExperimental(dat, method)
    rmseValue <- computeRMSEExperimental(dat)
    return (bindRMSECorrelation(rmseValue, corr))
}

#' reformat
#'
#' This function reformat the sequence to remove inconsistencies
#'
#' @param seq sequence
#'

reformat <- function(seq){
    seq <- gsub("X", "", seq)
    seq <- gsub("B", "", seq)
    seq <- gsub("J", "", seq)
    seq <- gsub("Z", "", seq)
    seq <- gsub("U", "", seq)
    return (seq)
}

#' removePTM
#'
#' This function reformat the sequence to remove all PTM labers.
#' Useful for example to compute AAindex Isoelectric point
#'
#' @param seq sequence
#'

removePTM <- function(seq){
    seq <- gsub("o", "", seq)
    seq <- gsub("m", "", seq)
    seq <- gsub("n", "", seq)
    seq <- gsub("p", "", seq)
    return (seq)
}


#' processTerminalSequence
#'
#' This function reformat the sequence removing Search Engine notation.
#' 
#' Example:
#' 
#' sequence <- "K.SDFGHQASSR.L"
#' s <- processTerminalSequence(seq = sequence)
#' 
#' The result will be:
#' s = SDFGHQASSR
#' 
#' @param seq sequence
#' 
processTerminalSequence <- function (seq){
  
  before_dot <- 3
  after_dot   <- nchar(seq) - 2
  seq <- substr(seq, start = before_dot, stop = after_dot)
  return (seq)
}

