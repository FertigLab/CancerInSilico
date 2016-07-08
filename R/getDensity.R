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

setMethod("getDensity", "CellMatrix",

    function(mat,time) {
      radii <- seq(3,ncol(mat),6)
      xcoords <- seq(1,ncol(mat),6)
      ycoords <- seq(2,ncol(mat),6)
      #farthest distance from (0,0) of cell
      d <- max(mat[time,xcoords])**2 + max(mat[time,ycoords])**2
      border <- 1;
      
      while(d > border^2){
        border <- border + 1
      }
      dens <- (sum(pi * mat[time,radii] ** 2))/(pi*border**2)
      return(dens)
    }
    
)
