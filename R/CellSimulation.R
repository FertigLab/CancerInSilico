#' Simulates Cell Model
#' @export
#'
#' @details This function provides a centralized R interface to run c++
#' code for cell-based models implemented in this package. Standard
#' parameters, as well as model-specific parameters, are passed in to this
#' function along with a model name. This function then runs the model and
#' returns a CellModel object containing all of the information from the
#' model. This object can then be accessed with various functions designed
#' to interact with the class. To see a list of available functions, there
#' is a show() command implemented for CellModel objects.
#' @param initialNum how many cells initially
#' @param runTime how long the simulation runs in real cellular time (hours)
#' @param cellTypes an array of S4 objects (each representing a different
#' cell type)
#' @param density the density the cells are seeded at
#' @param cycleLengthDist cycle time distribution
#' @param drugEffect distribution of drug effects
#' @param inheritGrowth whether or not daughter cells have the same
#' cycle-length as parents
#' @param outputIncrement time increment to print status at
#' @param randSeed seed for the model
#' @param modelType the name of the cell-based model to use
#' @param cycleSyncProb the probability of cells being seed in interphase
#' (not mitosis)
#' @param ... model specific parameters (depends on modelType)
#' @return A CellModel containing all info from the model run
#' @examples
#' inSilicoCellModel(initialNum=1, runTime=8, density=0.1)
inSilicoCellModel <- function(initialNum, runTime, density, 
modelType = 'DrasdoHohme', ...)
{
    # create cell model object
    model <- new('DrasdoHohmeModel', initialNum = initialNum, 
        runTime = runTime, density = density, ...)
    
    # run model and return result
    return(run(model))
}
