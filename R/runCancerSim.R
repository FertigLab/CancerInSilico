#'\code{runCancerSim} Runs the model
#'
#' @details This function provides a centralized R interface to run c++ code for cell-based
#' models implemented in this package. Standard parameters, as well as model-specific parameters,
#' are passed in to this function along with a model name. This function then runs the model 
#' and returns a CellModel object containing all the information from the model. This object
#' can then be accessed with various functions designed to interact with the class. To see a list
#' of available functions, there is a show() command implemented for CellModel objects.
#' @param initialNum how many cells initially
#' @param runTime how long the simulation represents in realtime
#' @param density the density the cells are seeded at
#' @param cycleLengthDist cycle time distribution
#' @param inheritGrowth whether or not daughter cells have the same cycle-length as parents
#' @param outputIncrement time increment to print status at
#' @param randSeed seed for the model
#' @param modelType the name of the cell-based model to use
#' @param ... model specific parameters (depends on modelType)
#' @return A CellModel containing all info from the model run
#' @examples
#' runCancerSim(1,8)
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

    if (modelType != "DrasdoHohme2003") {

      stop("invalid model type")

    } else {

      return (runDrasdoHohme(initialNum, runTime, density, cycleLengthDist,
                inheritGrowth, outputIncrement, randSeed, ...))

    }

}
