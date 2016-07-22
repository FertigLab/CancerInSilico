#' \code{getDensity} Gets the density of cells at a given time
#'
#' @param model A Cell Model
#' @param time The time of interest
#' @return The density of cells at that time (not quite the same as confluency)
#' @examples
#' getDensity(runModel(100,10),5)
#' @export

setGeneric("getDensity", function(model,time)
  standardGeneric("getDensity"))

#' \code{getDensity} Gets the density of cells at a given time
#'
#' @param model A Cell Model
#' @param time The time of interest
#' @return The density of cells at that time (not quite the same as confluency)
#' @examples
#' getDensity(runModel(100,10),5)
#' @export

setMethod("getDensity", "CellModel",

    function(model,time) {
      
      radii <- seq(3,length(model@cells[[time/model@parameters[6] + 1]]),6)
      xcoords <- seq(1,length(model@cells[[time/model@parameters[6] + 1]]),6)
      ycoords <- seq(2,length(model@cells[[time/model@parameters[6] + 1]]),6)
      
      #farthest distance from (0,0) of cell
      d <- max(sqrt(model@cells[[time/model@parameters[6] + 1]][xcoords]**2 + model@cells[[time/model@parameters[6] +1]][ycoords]**2))
               
      return(sum(model@cells[[(time/model@parameters[6]) + 1]][radii] ** 2) / (d ^ 2))
      
    }
    
)
