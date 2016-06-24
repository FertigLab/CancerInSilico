#'\code{plotCellsAtTime} Plots a CellMatrix at a certain point in time
#'
#'
#'@param mat A Cell Matrix
#'@param time The timestep at which to plot the matrix. Must be below
#' the specified max amount of timesteps
#'@example
#'plotCellsAtTime(runModel(24, 23), 21)
#'@export

setGeneric("plotCellsAtTime", function(mat,time)
  standardGeneric("plotCellsAtTime"))

setMethod("plotCellsAtTime", "CellMatrix",
  function(mat,time)
    {
    radii = seq(3,ncol(mat),6)
    numCells = sum(mat[time,radii]>0)
    
    #Information of cells based on Matrix Values (Used in plotCellsAtTime)
    
    xcoords = seq(1,(numCells-1)*7,6)
    ycoords = xcoords + 1
    radii = ycoords + 1
    axis_len = radii + 1
    axis_ang = axis_len + 1
    
    mn = min(min(mat[,xcoords]),min(mat[,ycoords])) - 2
    mx = max(max(mat[,xcoords]),max(mat[,ycoords])) + 2
    
    #Opens new device to put plot into
    dev.new()
    dev.set(which = 1)
    
    plot(c(mn,mx),c(mn,mx),type="n")
    text(20,20,labels = "test")
    
    #How Many Cells are Alive
    for (n in xcoords) {
      #Currently Assuming All Cells are Alive (No Cell Death)
      x_1 =  mat[time,n] + (- 0.5 * mat[time,n+3] + mat[time,n+2]) * cos(mat[time,n+4])
      y_1 =  mat[time,n+1] + (- 0.5 * mat[time,n+3] + mat[time,n+2]) * sin(mat[time,n+4])
      x_2 =  mat[time,n] + (0.5 * mat[time,n+3] - mat[time,n+2]) * cos(mat[time,n+4])
      y_2 =  mat[time,n+1] + (0.5 * mat[time,n+3] - mat[time,n+2]) * sin(mat[time,n+4])
      AddCircle(x_1,y_1,mat[time,n+2])
      AddCircle(x_2,y_2,mat[time,n+2])
    }
    
  }
)