#'\code{runModel} Runs the model
#'
#'
#' @param initialNum how many cells initially
#' @param runTime how long the simulation runs
#' @param density the density the cells are seeded at
#' @param meanGrowth mean growth rate of cells
#' @param varGrowth the variance of the growth rate (assumed normal)
#' @param maxMigration farthest a cell can move in a single update
#' @param maxDeform most a cell can deform in a single update
#' @param maxRotate most a cell can rotate in a single update
#' @param epsilon model specific parameter
#' @param delta model specific parameter
#' @param outIncrement time increment to print status at
#' @param randSeed seed for the model
#' @return A CellModel containing all info from the model run
#' @examples
#' runModel(100,10)
#' @export

runModel <- function(initialNum,
                     runTime,
                     density = 0.01,
                     maxMigration = 0.75,
                     maxDeform = 0.1,
                     maxRotate = 0.8,
                     epsilon = 0.05,
                     delta = 5.0,
                     outIncrement = 10,
                     randSeed = 0,
                     growthRates = rep(0.1,initialNum),
                     inheritGrowth = F)

{
  
	if (density > 0.4) {stop()}
	if (length(growthRates) != initialNum) {stop()}

    output = tryCatch({

        CellModel(initialNum, runTime, density, maxMigration,
		maxDeform, maxRotate, epsilon, delta, outIncrement,
		randSeed, growthRates, inheritGrowth)

    }, error = function(cond) {

        message(cond, '\n')
        stop()

    })

    cellMat <- new("CellModel", output)

    return(cellMat)

}
