#' \code{getDensity} Gets the density of cells at a given time
#'
#' @param model A Cell Model
#' @param time The time of interest
#' @return The density of cells at that time (not quite the same as confluency)
#' @examples
#' getDensity(runCancerSim(1,10),10)
#' @export

setGeneric("getDensity", function(model,time)
  standardGeneric("getDensity"))

#' \code{getDensity} Gets the density of cells at a given time
#'
#' @param model A Cell Model
#' @param time The time of interest
#' @return The density of cells at that time (not quite the same as confluency)
#' @examples
#' getDensity(runCancerSim(1,10),10)
#' @rdname CellModel-class
#' @export

setMethod("getDensity", "CellModel",

    function(model,time) {

    row <- timeToRow(model,time)
    xcoords <- seq(1,length(model@m_cells[[row]]),6)
    ycoords <- xcoords + 1
    radii <- xcoords + 2
    

    #farthest distance from (0,0) of cell
    d <- max(sqrt(model@m_cells[[row]][xcoords]**2 + model@m_cells[[row]][ycoords]**2))
       
    return(sum(model@m_cells[[row]][radii] ** 2) / (d ^ 2))

}
    
)
