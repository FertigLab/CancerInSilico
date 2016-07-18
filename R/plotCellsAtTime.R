#' \code{plotCellsAtTime} Plots a CellModel at a certain point in time
#'
#' @param model A CellModel
#' @param time The timestep at which to plot the modelrix.
#' @return Plot a visual representation of cells at time
#' @examples
#' plotCellsAtTime(runModel(10,100),60)
#' @export

setGeneric("plotCellsAtTime", function(model,time)
    standardGeneric("plotCellsAtTime"))

#' \code{plotCellsAtTime} Plots a CellModel at a certain point in time
#'
#' @param model A CellModel
#' @param time The timestep at which to plot the modelrix.
#' @return Plot a visual representation of cells at time
#' @examples
#' plotCellsAtTime(runModel(10,100),60)
#' @export

setMethod("plotCellsAtTime", "CellModel",

	function(model,time)  {

    radii <- seq(3,ncol(model),6)
    
    numCells <- sum(model[time,radii]>0)

    xcoords <- seq(1,numCells * 6,6)

    mn <- min(min(model[,xcoords]),min(model[,xcoords+1])) - 2
    mx <- max(max(model[,xcoords]),max(model[,xcoords+1])) + 2

    plot(c(mn,mx),c(mn,mx),main=paste("Plot of CellModel At Time",time),type="n",asp=1)
          
    x_1 <- model[time,xcoords] + (0.5 * model[time,xcoords+3] - model[time,xcoords+2]) * cos(model[time,xcoords+4])
    x_2 <- model[time,xcoords] - (0.5 * model[time,xcoords+3] - model[time,xcoords+2]) * cos(model[time,xcoords+4])
    y_1 <- model[time,xcoords+1] + (0.5 * model[time,xcoords+3] - model[time,xcoords+2]) * sin(model[time,xcoords+4])
    y_2 <- model[time,xcoords+1] - (0.5 * model[time,xcoords+3] - model[time,xcoords+2]) * sin(model[time,xcoords+4])
    
    x <- c(x_1,x_2)
    y <- c(y_1,y_2)
    rad <- c(model[time,xcoords+2], model[time,xcoords+2])
    
    symbols(x,y, circles=rad, inches=FALSE, add=TRUE, bg="bisque4", fg="bisque4")
   
  }

)
