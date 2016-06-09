runModel <- function(initialNum, runTime, density) {
  
  data <- CellModel(initialNum, runTime, density)
  cellMat <- new("cellMatrix", data)
  
  return(cellMat)
  
}
