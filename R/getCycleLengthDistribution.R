#' \code{getCycleLengthDistribution} Plots a distribution of the growth
#'      rates of all the cells
#'
#' @param model A CellModel
#' @param time The time of interest for the growth rates
#' @return returns a vector of growth rates for cells alive at time
#' @examples
#' getCycleLengthDistribution(runCancerSim(1,10),0)
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
#' getCycleLengthDistribution(runCancerSim(1,10),0)
#' @rdname CellModel-class
#' @export

setMethod("getCycleLengthDistribution", "CellModel",

    function(model,time) {
          
        nG <- .nG(model)
        timeIncrement <- .timeIncrement(model)
        row <- timeToRow(model,time)
        gr_rates <- seq(6,length(model@m_cells[[row]]),6)
        gr <- model@m_cells[[row]][gr_rates]
        ret_val <- 1 + 2 * (sqrt(2) - 1) * timeIncrement * nG / gr
        return(ret_val)

    }
          
)
