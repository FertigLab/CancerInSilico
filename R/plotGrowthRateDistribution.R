#' \code{plotGrowthRateDistribution} Plots a distribution of the growth
#'      rates of all the cells
#'
<<<<<<< HEAD
#'@param mat A Cell Matrix
#'@param time The time of interest for the growth rates
#'@export
=======
>>>>>>> 1f571cc86bdbc85a7bd5f5e6d24cf291a11c89a3
#' @param mat A Cell Matrix
#' @param time The time of interest for the growth rates
#' @return Plots the distribution of growth rates at time
#' @examples
#' plotGrowthRateDistribution(runModel(100,10),10)
#' @export
<<<<<<< HEAD

=======
>>>>>>> 1f571cc86bdbc85a7bd5f5e6d24cf291a11c89a3

setGeneric("plotGrowthRateDistribution", function(mat, time)
    standardGeneric("plotGrowthRateDistribution"))

<<<<<<< HEAD
=======
#' \code{plotGrowthRateDistribution} Plots a distribution of the growth
#'      rates of all the cells
#'
#' @param mat A Cell Matrix
#' @param time The time of interest for the growth rates
#' @return Plots the distribution of growth rates at time
#' @export
>>>>>>> 1f571cc86bdbc85a7bd5f5e6d24cf291a11c89a3

setMethod("plotGrowthRateDistribution", "CellMatrix",

    function(mat, time) {

        radii = seq(3,ncol(mat),6)
        gr_rates = radii[mat[time,radii]>0] + 3
        plot(density(mat[time,gr_rates]))

    }
          
)
