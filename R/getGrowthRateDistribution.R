#' \code{getGrowthRateDistribution} Plots a distribution of the growth
#'      rates of all the cells
#'
#' @param model A CellModel
#' @param time The time of interest for the growth rates
#' @return returns a vector of growth rates for cells alive at time
#' @examples
#' getGrowthRateDistribution(runModel(100,10),10)
#' @export

setGeneric("getGrowthRateDistribution", function(model, time)
    standardGeneric("getGrowthRateDistribution"))

#' \code{getGrowthRateDistribution} Plots a distribution of the growth
#'      rates of all the cells
#'
#' @param model A CellModel
#' @param time The time of interest for the growth rates
#' @return returns a vector of growth rates for cells alive at time
#' @examples
#' getGrowthRateDistribution(runModel(100,10),10)
#' @export

setMethod("getGrowthRateDistribution", "CellModel",

    function(model, time) {

        radii = seq(3,ncol(model),6)
        gr_rates = radii[model[time,radii]>0] + 3
        return(model[time,gr_rates])

    }
          
)
