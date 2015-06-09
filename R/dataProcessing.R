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
    return (values)
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

computeRMSEExperimental <- function(data){
    nm <- names(dat)
    rmseValue <- data.frame(matrix(ncol = ncol(dat), nrow = 1))
    colnames(rmseValue) <- nm
    rmseValue <- removeFirstColumn(rmseValue)
    for (i in nm) {
        if(dat[1L] != data[i]){
            newData <- data.frame(x = dat[1L], y = dat[i])
            colnames(newData) <- c("x", "y")
            value <- rmse(newData$x, newData$y)
            rmseValue[i] <- value
        }
    }
    return (rmseValue)
}

