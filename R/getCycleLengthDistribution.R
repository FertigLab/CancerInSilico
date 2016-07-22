#' \code{getCycleLengthDistribution} Plots a distribution of the growth
#'      rates of all the cells
#'
#' @param model A CellModel
#' @param time The time of interest for the growth rates
#' @return returns a vector of growth rates for cells alive at time
#' @examples
#' getCycleLengthDistribution(runModel(100,10),10)
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
#' getCycleLengthDistribution(runModel(100,10),10)
#' @export

setMethod("getCycleLengthDistribution", "CellModel",

    function(model, time) {
    
        row <- timeToRow(model,time)
        gr_rates <- seq(6,length(model@cells[[row]]),6)
        gr <- model@cells[[row]][gr_rates]
        nG <- 10 ## must be up-to-date with line 36 in runModel
        t <- model@parameters[6]
        ret_val <- 3 * (sqrt(2) - 1) * t * nG / gr
        return(ret_val)

    }
          
)
