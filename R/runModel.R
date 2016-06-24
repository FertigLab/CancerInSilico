#' Creates a Cell Matrix based on certain simulated variables
#' 
#' @param initialNum A number greater than 0
#' @param runTime A number greater than 0
#' @param density 
#' @param meanGrowth 
#' @param varGrowth 
#' @param maxMigration 
#' @param maxDeform 
#' @param maxRotate 
#' @param epsilon 
#' @param delta
#' @param outIncrement
#' @param randseed 

runModel <- function(initialNum, runTime, density = 0.05,
                     meanGrowth = 0.15, varGrowth = 0.0, maxMigration = 0.5,
                     maxDeform = 0.075, maxRotate = 0.3, epsilon = 0.05,
                     delta = 5.0, outIncrement = 10, randSeed = 0)
{
  
 output = tryCatch({
 
   CellModel(initialNum, runTime, density, meanGrowth,
     varGrowth, maxMigration, maxDeform, maxRotate,
     epsilon, delta, outIncrement, randSeed)
 
 }, error = function(cond) {

   message(cond,'\n')
   stop()
    
 })
  
  cellMat <- new("cellMatrix", output)

  return(cellMat)
  
}
