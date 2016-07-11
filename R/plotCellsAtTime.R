#' \code{plotCellsAtTime} Plots a CellModel at a certain point in time
#'
#' @param mat A CellModel object
#' @param time The timestep at which to plot the matrix. Must be below
#'      the specified max amount of timesteps
#' @return Plot a visual representation of cells at time
#' @export

setGeneric("plotCellsAtTime", function(mat,time)
    standardGeneric("plotCellsAtTime"))

setMethod("plotCellsAtTime", "CellMatrix",

    function(mat,time)  {

        radii = seq(3,ncol(mat),6)
        #Currently Assuming All Cells are Alive (No Cell Death)
        numCells = sum(mat[time,radii]>0)

        xcoords = seq(1,numCells * 6,6)
        ycoords = xcoords + 1
        radii = ycoords + 1
        axis_len = radii + 1
        axis_ang = axis_len + 1
        
        mn = min(min(mat[,xcoords]),min(mat[,ycoords])) - 2
        mx = max(max(mat[,xcoords]),max(mat[,ycoords])) + 2

        #Opens new device to put plot into if no plot exists already
        
        if(identical(dev.list(),NULL) == TRUE){
          dev.new()
          dev.set(which = 1)
        }
        dev.set(which = 4)
        
        theta = seq(0,2*pi,length=20)
        plot(c(mn,mx),c(mn,mx),main=paste("Plot of CellModel At Time",time),type="n")
          
        for (n in xcoords) {
  
            x_1 =  mat[time,n] + (0.5 * mat[time,n+3] - mat[time,n+2]) *
                    cos(mat[time,n+4])
            y_1 =  mat[time,n+1] + (0.5 * mat[time,n+3] - mat[time,n+2]) *
                    sin(mat[time,n+4])
            x_2 =  mat[time,n] - (0.5 * mat[time,n+3] - mat[time,n+2]) *
                    cos(mat[time,n+4])
            y_2 =  mat[time,n+1] - (0.5 * mat[time,n+3] - mat[time,n+2]) *
                    sin(mat[time,n+4])
            
            lines(x_1 + mat[time,n+2] * cos(theta),y_1 + mat[time,n+2] * sin(theta),type="l",new=FALSE)
            lines(x_2 + mat[time,n+2] * cos(theta),y_2 + mat[time,n+2] * sin(theta),type="l",new=FALSE)

        }
        
      }

)
