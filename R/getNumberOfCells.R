#' \code{getNumberOfCells} return the total number of cells at time
#'
#' @param model A CellModel
#' @param time the time of interest
#' @return the size of the cell population
#' @examples
#' getNumberOfCells(runModel(10,10),10)
#' @export
#' 
setGeneric("getNumberOfCells", function(model, time)
    standardGeneric("getNumberOfCells"))

#' \code{getNumberOfCells} Plots the total number of cells vs time
#'
#' @param model A CellModel
#' @param time the time of interest
#' @return the size of the cell population over time
#' @examples
#' getNumberOfCells(runModel(10,10),10)
#' @export

setMethod("getNumberOfCells", "CellModel",

    function(model, time) {

        radii = seq(3, length(model@cells[[timeToRow(model,time)]]), 6)
        return(sum(model@cells[[timeToRow(model,time)]][radii] > 0))
    
    }

)
