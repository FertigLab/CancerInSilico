#' \code{getDensity} Gets the density of cells at a given time
#'
#' @param mat A Cell Model object
#' @param time The time of interest
#' @return The density of cells at that time
#' @examples
#' getDensity(runModel(100,10),5)
#' @export

setGeneric("getDensity", function(mat,time)
  standardGeneric("getDensity"))

#' \code{getDensity} Gets the density of cells at a given time
#'
#' @param mat A Cell Model object
#' @param time The time of interest
#' @return The density of cells at that time
#' @export

setMethod("getDensity", "CellMatrix",

    function(mat,time) {

        return(1)

    }
    
)
