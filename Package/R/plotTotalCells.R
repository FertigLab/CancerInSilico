setGeneric("plotTotalCells", function(mat)
            standardGeneric("plotTotalCells"))

setMethod("plotTotalCells", "cellMatrix",

  function(mat) {
  
    total_cells = c()
    radii = seq(3,ncol(mat),6)
  
    for (t in 1:nrow(mat)) {
      total_cells[t] = sum(mat[t,radii]>0)
    }
  
    plot(total_cells,type="l")
  }

)