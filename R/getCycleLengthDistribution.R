#' \code{getCycleLengthDistribution} Plots a distribution of the growth
#'      rates of all the cells
#'
#' @param model A CellModel
#' @param time The time of interest for the growth rates
#' @return returns a vector of growth rates for cells alive at time
#' @examples
#' getCycleLengthDistribution(runModel(1,10),0)
#' @export

setGeneric("getCycleLengthDistribution", function(model, time)
    standardGeneric("getCycleLengthDistribution"))

#' \code{getCycleLengthDistribution} Plots a distribution of the growth
#'      rates of all the cells
#'
#' @param model A CellModel
#' @param time The time of interest for the growth rates
#' @return returns a vector of growth rates for cells alive at time
#' @examples
#' getCycleLengthDistribution(runModel(1,10),0)
#' @export

setMethod("getCycleLengthDistribution", "CellModel",

    function(model,time) {
          
        nG <- model@parameters[8]
        timeIncrement <- model@parameters[9]
        row <- timeToRow(model,time)
        gr_rates <- seq(6,length(model@cells[[row]]),6)
        gr <- model@cells[[row]][gr_rates]
        ret_val <- 1 + 2 * (sqrt(2) - 1) * timeIncrement * nG / gr
        return(ret_val)

    }
          
)
