#'\code{runCancerSim} Runs the model
#'
#'
#' @param initialNum how many cells initially
#' @param runTime how long the simulation represents in real cellular time (hours)
#' @param density the density the cells are seeded at
#' @param cycleLengthDist cycle time distribution
#' @param drugEffect distribution of drug effects
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
                     drugEffect = getDrugEffect(cycleLengthDist = cycleLengthDist),
                     inheritGrowth = FALSE,
                     outputIncrement = 6,
                     recordIncrement = 0.25,
                     randSeed = 0,
                     modelType = "DrasdoHohme2003",
                     drugTime = 0.0,
                     ...)

{

    if (modelType != "DrasdoHohme2003") {
      stop("invalid model type")
    } else {
      return (runDrasdoHohme(initialNum, runTime, density, cycleLengthDist, inheritGrowth, outputIncrement, recordIncrement, randSeed, drugEffect, drugTime, ...))
    }

}
