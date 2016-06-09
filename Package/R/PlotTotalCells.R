plotTotalCells <-function(dataMatrix) {
  
  total_cells = c()
  radii = seq(3,ncol(dataMatrix),6)
  
  for (t in 1:nrow(dataMatrix)) {
    total_cells[t] = sum(dataMatrix[t,radii]>0)
  }
  plot(total_cells)
  
}