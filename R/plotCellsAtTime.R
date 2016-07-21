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
	    modtime <- (time / model@parameters[5]) + 1
        
	    radii <- seq(1,ncol(model@cells),6)
        
        numCells <- sum(model@cells[modtime,radii]>0)
    
        xcoords <- seq(1,numCells * 6,6)
    
        mn <- min(min(model@cells[,xcoords]),min(model@cells[,xcoords+1])) - 2
        mx <- max(max(model@cells[,xcoords]),max(model@cells[,xcoords+1])) + 2
        
        plot(c(mn,mx),c(mn,mx),main=paste("Plot of CellModel At Time",time),type="n",asp=1)
              
        x_1 <- model@cells[modtime,xcoords] + (0.5 * model@cells[modtime,xcoords+3] - model@cells[modtime,xcoords+2]) * cos(model@cells[modtime,xcoords+4])
        x_2 <- model@cells[modtime,xcoords] - (0.5 * model@cells[modtime,xcoords+3] - model@cells[modtime,xcoords+2]) * cos(model@cells[modtime,xcoords+4])
        y_1 <- model@cells[modtime,xcoords+1] + (0.5 * model@cells[modtime,xcoords+3] - model@cells[modtime,xcoords+2]) * sin(model@cells[modtime,xcoords+4])
        y_2 <- model@cells[modtime,xcoords+1] - (0.5 * model@cells[modtime,xcoords+3] - model@cells[modtime,xcoords+2]) * sin(model@cells[modtime,xcoords+4])
        
        x <- c(x_1,x_2)
        y <- c(y_1,y_2)
        rad <- c(model@cells[modtime,xcoords+2], model@cells[modtime,xcoords+2])
        
        symbols(x,y, circles=rad, inches=FALSE, add=TRUE, bg="bisque4", fg="bisque4")
   
  }

)
