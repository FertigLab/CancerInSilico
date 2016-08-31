#' \code{plotCellsAtTime} Plots a CellModel at a certain point in time
#'
#' @param model A CellModel
#' @param time The timestep at which to plot the modelrix.
#' @return Plot a visual representation of cells at time
#' @examples
#' plotCellsAtTime(runCancerSim(1,10),10)
#' @export

setGeneric("plotCellsAtTime", function(model,time)
    standardGeneric("plotCellsAtTime"))

#' \code{plotCellsAtTime} Plots a CellModel at a certain point in time
#'
#' @param model A CellModel
#' @param time The timestep at which to plot the modelrix.
#' @return Plot a visual representation of cells at time
#' @examples
#' plotCellsAtTime(runCancerSim(1,10),10)
#' @rdname CellModel-class
#' @export

setMethod("plotCellsAtTime", "CellModel",

	function(model,time)  {

    	row <- timeToRow(model,time)
        
        numCells <- length(model@m_cells[[row]]) / 6
    
        xcoords <- seq(1,numCells * 6,6)
        ycoords <- xcoords + 1
        radii <- xcoords + 2
        axis_len <- xcoords + 3
        axis_ang <- xcoords + 4
    
        mn <- min(min(model@m_cells[[row]][xcoords]),min(model@m_cells[[row]][ycoords])) - 2
        mx <- max(max(model@m_cells[[row]][xcoords]),max(model@m_cells[[row]][ycoords])) + 2
        
        plot(c(mn,mx),c(mn,mx),main=paste("Plot of CellModel At Time",time),xlab = "",ylab="",type="n",asp=1)
              
        x_1 <- model@m_cells[[row]][xcoords] + (0.5 * model@m_cells[[row]][axis_len] - model@m_cells[[row]][radii]) * cos(model@m_cells[[row]][axis_ang])
        x_2 <- model@m_cells[[row]][xcoords] - (0.5 * model@m_cells[[row]][axis_len] - model@m_cells[[row]][radii]) * cos(model@m_cells[[row]][axis_ang])
        y_1 <- model@m_cells[[row]][ycoords] + (0.5 * model@m_cells[[row]][axis_len] - model@m_cells[[row]][radii]) * sin(model@m_cells[[row]][axis_ang])
        y_2 <- model@m_cells[[row]][ycoords] - (0.5 * model@m_cells[[row]][axis_len] - model@m_cells[[row]][radii]) * sin(model@m_cells[[row]][axis_ang])
        
        x <- c(x_1,x_2)
        y <- c(y_1,y_2)
        rad <- c(model@m_cells[[row]][radii], model@m_cells[[row]][radii])
        
        symbols(x,y, circles=rad, inches=FALSE, add=TRUE, bg="bisque4", fg="bisque4")
   
  }

)
