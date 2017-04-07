#'\code{runCancerSim} Runs a cancer simulation given parameters
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
#' runCancerSim(1,8)
#' @export

runCancerSim <- function(initialNum,
                         runTime,
                         cellTypes = c(newCellType('DEFAULT')),
                         cellTypeInitFreq = c(1),
                         drugs = NULL,
                         modelType = "DrasdoHohme",
                         density = 0.01,
                         boundary = TRUE,
                         randSeed = 0,
                         syncCycles = TRUE,
                         outputIncrement = 6,
                         recordIncrement = 0.25,
                         ...)
{
    # store parameters in a list
    params <- list()
    params[['initialNum']] <- initialNum
    params[['runTime']] <- runTime
    params[['cellTypes']] <- cellTypes
    params[['cellTypeInitFreq']] <- cellTypeInitFreq
    params[['modelType']] <- modelType
    params[['density']] <- density
    params[['boundary']] <- boundary
    params[['randSeed']] <- randSeed
    params[['syncCycles']] <- syncCycles
    params[['outputIncrement']] <- outputIncrement
    params[['recordIncrement']] <- recordIncrement

    # make sure all arguments are valid
    checkParameters(params, ...)

    # run model
    return (runModel(params, ...))
}
    
checkParameters <- function(params, ...)
{
    # general parameters check
    if (params[['density']] > 0.7)
    {
       stop("density too high to seed efficiently\n")
    }
}

runModel <- function(params, ...)
{
    # list of valid model types
    validMods <- c('DrasdoHohme2003')

    # check that model is valid type
    if (!(params[['modelType']] %in% validMods))
    {
      stop("invalid model type")
    }
    else if (params[['modelType']] == 'DrasdoHohme2003')
    {
        return (runDrasdoHohme(params, ...))
    }
}

runDrasdoHohme <- function(params, ...)
{  
    ## get model specific parameters
    params[['nG']] <- list(...)$nG
    params[['epsilon']] <- list(...)$epsilon
    params[['delta']] <- list(...)$delta
    params[['drugs']] <- list(...)$drugs

    ## set to defaults if not provided
    if (is.null(params[['nG']])) {params[['nG']] = 24}
    if (is.null(params[['epsilon']])) {params[['epsilon']] = 10}
    if (is.null(params[['delta']])) {params[['delta']] = 0.2}
  
    output <- runCellSimulation(params)
    cellMat <- createCellModel(output[['params']], output[['cells']])

    return(cellMat)
}

