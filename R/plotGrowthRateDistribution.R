#' \code{plotGrowthRateDistribution} Plots a distribution of the growth
#'      rates of all the cells
#'
#' @param mat A CellModel
#' @param time The time of interest for the growth rates
#' @return Plots the distribution of growth rates at time
#' @examples
#' plotGrowthRateDistribution(runModel(100,10),10)

setGeneric("plotGrowthRateDistribution", function(mat, time)
    standardGeneric("plotGrowthRateDistribution"))

#' \code{plotGrowthRateDistribution} Plots a distribution of the growth
#'      rates of all the cells
#'
#' @param mat A CellModel
#' @param time The time of interest for the growth rates
#' @return Plots the distribution of growth rates at time
#' @examples
#' plotGrowthRateDistribution(runModel(100,10),10)
#' @export

setMethod("plotGrowthRateDistribution", "CellMatrix",

    function(mat, time) {

        radii = seq(3,ncol(mat),6)
        gr_rates = radii[mat[time,radii]>0] + 3
        plot(density(mat[time,gr_rates]),main=paste("Plot of Growth Rate Distribution at Time",time),xlab = "Growth Rate",ylab = "Density")

    }
          
)
