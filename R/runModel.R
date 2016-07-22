#'\code{runModel} Runs the model
#'
#'
#' @param initialNum how many cells initially
#' @param runTime how long the simulation represents in realtime
#' @param density the density the cells are seeded at
#' @param cycleLengthDist cycle time distribution
#' @param inheritGrowth whether or not daughter cells have the same cycle-length as parents
#' @param timeIncrement time length of a step in the model, based on the hour
#' @param outputIncrement time increment to print status at
#' @param randSeed seed for the model
#' @param epsilon epsilon model specific parameter
#' @return A CellModel containing all info from the model run
#' @examples
#' runModel(100,10)
#' @export

runModel <- function(initialNum,
                     runTime,
                     density = 0.01,
                     cycleLengthDist = 12,
                     inheritGrowth = F,
                     timeIncrement = 0.25,
                     outputIncrement = 40,
                     randSeed = 0,
                     epsilon = 4)
                         
{

    if (density > 0.1) {
        
        message("density too high to seed efficiently\n")
        stop()

    }

    nG <- 10
    delta <- 0.5
    grRates <- 3 * (sqrt(2) - 1) * timeIncrement * nG / cycleLengthDist
    mcSteps <- ceiling(runTime / timeIncrement)
    maxTranslation <- mean(grRates)
    maxRotate <- 3 * mean(grRates)
    maxDeform <- 2 * max(grRates)

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
