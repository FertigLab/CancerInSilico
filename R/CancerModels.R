#'\code{runCancerSim} runs a cell-based model of cancer
#'
#' @details This function provides a centralized R interface to run c++ code for cell-based models implemented in this package. Standard parameters, as well as model-specific parameters, are passed in to this function along with a model name. This function then runs the model and returns a CellModel object containing all the information from the model. This object can then be accessed with various functions designed to interact with the class. To see a list of available functions, there is a show() command implemented for CellModel objects.
#' @param initialNum how many cells initially (integer)
#' @param runTime how long the simulation runs (model hours)
#' @param density the density the cells are seeded at, must be in (0,0.1]
#' @param cycleLengthDist cycle time distribution
#' @param inheritGrowth whether or not daughter cells have the same cycle-length as parents
#' @param outputIncrement time increment to print status at
#' @param randSeed seed for the model
#' @param modelType the name of the cell-based model to use
#' @param ... model specific parameters (depends on modelType)
#' @return A CellModel containing all info from the model run
#' @examples
#' runCancerSim(1,4)
#' @export

runCancerSim <- function(initialNum,
                     runTime,
                     density = 0.01,
                     cycleLengthDist = 12,
                     inheritGrowth = FALSE,
                     outputIncrement = 6,
                     randSeed = 0,
                     modelType = "DrasdoHohme2003",
                     ...)

{

    # check if model type is valid
    if (modelType != "DrasdoHohme2003") {

      stop("invalid model type")

    } else {

        # call function that implements specific model
        return (runDrasdoHohme(initialNum, runTime, density, cycleLengthDist,
                inheritGrowth, outputIncrement, randSeed, ...))

    }

}

#'\code{runDrasdoHohme} runs the model based on Drasdo and Hohme (2003)
#'
#' @details This function calls the C++ implementation of the Drasdo and Hohme (2003) model.
#' @param initialNum how many cells initially
#' @param runTime how long the simulation represents in realtime
#' @param density the density the cells are seeded at
#' @param cycleLengthDist cycle time distribution
#' @param inheritGrowth whether or not daughter cells have the same cycle-length as parents
#' @param outputIncrement time increment to print status at
#' @param randSeed seed for the model
#' @param ... nG, epsilon parameters (specific to this model)
#' @return A CellModel containing all info from the model run

runDrasdoHohme <- function(initialNum,
                         runTime,
                         density,
                         cycleLengthDist,
                         inheritGrowth,
                         outputIncrement,
                         randSeed,
                         ...)
  
{
  
    # process model specific parameters
    nG <- list(...)$nG
    epsilon <- list(...)$epsilon

    # set to default value if  missing
    if (is.null(nG)) {nG = 24}
    if (is.null(epsilon)) {epsilon = 10}

    # check if density is valid  
    if (density > 0.1) {

        message("density too high to seed efficiently\n")
        stop()

    }
  
    # model specific, must be less than 4 or calculations break
    delta <- 0.2 

    # calculate model time increment based on parameters
    timeIncrement = delta / (4 * nG * (4 - sqrt(2)))

    # calculate the maximum time increment that still allows for accurate model
    max_incr = delta * (min(cycleLengthDist) - 1) / (8 * nG * (sqrt(2) - 1))

    # enforce the maximum on the time increment
    if (timeIncrement > max_incr) {

        timeIncrement = max_incr

    }

    # calculate parameters internal to the model
    maxDeform <- 2 * timeIncrement * nG * (4 - sqrt(2))
    grRates <- 2 * (sqrt(2) - 1) * timeIncrement * nG / (cycleLengthDist - 1)
    mcSteps <- ceiling(runTime / timeIncrement)
    maxTranslation <- delta / 2
    maxRotate <- acos((16 + delta ^ 2 - 4 * delta) / 16)
    outputIncrement2 <- floor(outputIncrement / timeIncrement)

    # call the model, this function implemented in C++ code    
    output <- CellModel(initialNum, mcSteps, density, maxTranslation,
          maxDeform, maxRotate, epsilon, delta, outputIncrement2,
          randSeed, grRates, inheritGrowth, nG, timeIncrement)

    # create CellModel object with parameters and data from the simulation
    cellMat <- new("CellModel",
                mCells = output,
                mInitialNumCells = initialNum,
                mRunTime = runTime,
                mInitialDensity = density,
                mInheritGrowth = inheritGrowth,
                mOutputIncrement = outputIncrement,
                mRandSeed = randSeed,
                mEpsilon = epsilon,
                mNG = nG,
                mTimeIncrement = timeIncrement,
                mCycleLengthDist = cycleLengthDist)

    # return CellModel object
    return(cellMat)

}
