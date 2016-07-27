#'\code{runModel} Runs the model
#'
#'
#' @param initialNum how many cells initially
#' @param runTime how long the simulation represents in realtime
#' @param density the density the cells are seeded at
#' @param cycleLengthDist cycle time distribution
#' @param inheritGrowth whether or not daughter cells have the same cycle-length as parents
#' @param outputIncrement time increment to print status at
#' @param randSeed seed for the model
#' @param epsilon epsilon model specific parameter
#' @param nG model specific parameter
#' @return A CellModel containing all info from the model run
#' @examples
#' runModel(100,10)
#' @export

runModel <- function(initialNum,
                     runTime,
                     density = 0.01,
                     cycleLengthDist = 12,
                     inheritGrowth = F,
                     outputIncrement = 6,
                     randSeed = 0,
                     epsilon = 4,
                     nG = 10)
    
{
    
    if (density > 0.1) {
        
        message("density too high to seed efficiently\n")
        stop()
        
    }
    
    delta <- 0.2 ## must be less than 4 or calculations break
    
    timeIncrement = delta / (4 * nG * (4 - sqrt(2)))
    if (timeIncrement > delta * (min(cycleLengthDist) - 1) / (8 * nG * (sqrt(2) - 1))) {
        timeIncrement = delta * (min(cycleLengthDist) - 1) / (8 * nG * (sqrt(2) - 1))
    }
    maxDeform <- 2 * timeIncrement * nG * (4 - sqrt(2))
    grRates <- 2 * (sqrt(2) - 1) * timeIncrement * nG / (cycleLengthDist - 1)
    mcSteps <- ceiling(runTime / timeIncrement)
    maxTranslation <- delta / 2
    maxRotate <- acos((16 + delta ^ 2 - 4 * delta) / 16)
    outputIncrement <- floor(outputIncrement / timeIncrement)
    
    output <- tryCatch({
        
        CellModel(initialNum, mcSteps, density, maxTranslation,
                  maxDeform, maxRotate, epsilon, delta, outputIncrement,
                  randSeed, grRates, inheritGrowth, nG, timeIncrement)
        
    }, error = function(cond) {
        
        message(cond, '\n')
        stop()
        
    })
    
    cellMat <- new("CellModel",cells = output,parameters = c(initialNum,runTime,density,mean(cycleLengthDist),inheritGrowth,timeIncrement,outputIncrement,randSeed,epsilon))
    
    return(cellMat)
    
}