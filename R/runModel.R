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
					 cycleTimeDist = 12,
					 inheritGrowth = F,
					 timeIncrement = 0.25,
					 outputIncrement = 40,
				     randSeed = 0)
 						
{

	if (density > 0.4) {
		
		message("density too high to seed efficiently\n")
		stop()

	}

	nG <- 10
  
	if (timeIncrement > cycleTimeDist / (2 * nG)) {

		timeIncrement <- cycleTimeDist / (2 * nG)
	
	}

	mcSteps <- ceiling(runTime / timeIncrement)
	grRates <- (sqrt(2) - 1) * timeIncrement * 2 * nG / cycleTimeDist
	maxTranslation <- mean(grRates) / 5
	maxRotate <- 10 * mean(grRates) / 11  
	epsilon <- 0.88
	delta <- 0.5
	maxDeform <- (sqrt(2) - 1) / 5

    output = tryCatch({

        CellModel(initialNum, mcSteps, density, maxTranslation,
		maxDeform, maxRotate, epsilon, delta, outputIncrement,
		randSeed, grRates, inheritGrowth, nG, timeIncrement)

    }, error = function(cond) {

        message(cond, '\n')
        stop()

    })

    cellMat <- new("CellModel", output)

    return(cellMat)

}
