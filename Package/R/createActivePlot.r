setGeneric("createActivePlot", function(mat,time)
  standardGeneric("createActivePlot"))




setMethod("createActivePlot", "cellMatrix",
  function(mat,time)
    {

    dev.new()
    dev.set(which = 1)
    dev.next(which = dev.cur())
    plot(c(mn,mx),c(mn,mx),type="n")
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