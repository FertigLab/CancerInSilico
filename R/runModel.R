#' \code{runModel} Creates a Cell Matrix based on certain simulated variables
#' 
#' 
#' @param initialNum A number representing initial number of cells
#' @param runTime A number representing max number of timesteps
#' @param density A decimal for density of cells
#' @param meanGrowth A decimal representing mean growth rate of cells
#' @param varGrowth 
#' @param maxMigration 
#' @param maxDeform 
#' @param maxRotate 
#' @param epsilon 
#' @param delta 
#' @param outIncrement A number representing increments 
#'  between timesteps
#' @param randseed A seed for the random number generator
#' @examples 
#' runModel(100,100)
#' @export

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

   message(cond, '\n')
   stop()
    
 })
  
  cellMat <- new("cellMatrix", output)

  return(cellMat)
  
}
