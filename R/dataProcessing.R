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
