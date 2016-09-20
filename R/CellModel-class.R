#' An S4 class to represent the output of a cell-based model
#'
#' @slot m_cells A list object where each row of the list describes the state
#' of all the cells in the model at a given time. Each cell is described over
#' 6 columns: [1] x-coordinate, [2] y-coordinate, [3] radius, [4] axis length,
#' [5] axis angle, [6] growth rate. For instance, the x-coordinates of the first
#' 3 cells will be in columns 1,7,13.
#' @slot m_initialNumCells the initial number of cells in the model
#' @slot m_runTime the total run time (hours) of the model
#' @slot m_initialDensity the density the cells were seeded at
#' @slot m_inheritGrowth whether or not cells inherit growth rates from their parent
#' @slot m_outputIncrement the frequency of print statements during the run
#' @slot m_randSeed the random seed 
#' @slot m_epsilon model specific parameter 
#' @slot m_nG model specific parameter
#' @slot m_timeIncrement amount of time elapsed in each model step
#' @slot m_cycleLengthDist initial distribution of cell-cycle lengths 
#' @export

setClass("CellModel", representation(
                        mCells = "list",
                        mInitialNumCells = "numeric",
                        mRunTime = "numeric",
                        mInitialDensity = "numeric",
                        mInheritGrowth = "logical",
                        mOutputIncrement = "numeric",
                        mRandSeed = "numeric",
                        mEpsilon = "numeric",
                        mNG = "numeric",
                        mTimeIncrement = "numeric",
                        mCycleLengthDist = "numeric",
                        mRecordIncrement = "numeric" ))

## getters (parameters)

.initialNumCells <- function(model) {

    return (model@mInitialNumCells)

}

.runTime <- function(model) {

    return (model@mRunTime)

}

.initialDensity <- function(model) {

    return (model@mInitialDensity)

}

.inheritGrowth <- function(model) {

    return (model@mInheritGrowth)

}

.outputIncrement <- function(model) {

    return (model@mOutputIncrement)

}

.randSeed <- function(model) {

    return (model@mRandSeed)

}

.epsilon <- function(model) {

    return (model@mEpsilon)

}

.nG <- function(model) {

    return (model@mNG)

}

.timeIncrement <- function(model) {

    return (model@mTimeIncrement)

}

.recordIncrement <- function(model) {

    return (model@mRecordIncrement)

}

.cycleLengthDist <- function(model) {

    return (model@mCycleLengthDist)

}

## getters (cell data)

getCoordinates <- function(model, time) {

    indices <- seq(1,length(model@mCells[[timeToRow(model,time)]]),6)
    ret_mat <- matrix(nrow = length(indices), ncol = 2)
    for (i in 1:length(indices)) {

        ret_mat[i,1] = model@mCells[[timeToRow(model,time)]][indices[i]]
        ret_mat[i,2] = model@mCells[[timeToRow(model,time)]][indices[i]+1]

	}

    return (ret_mat)
         
}

getRadii <-	function(model, time) {

    indices <- seq(3,length(model@mCells[[timeToRow(model,time)]]),6)
    return(model@mCells[[timeToRow(model,time)]][indices])

}
         
getAxisLength <- function(model, time) {

    indices <- seq(4,length(model@mCells[[timeToRow(model,time)]]),6)
    return(model@mCells[[timeToRow(model,time)]][indices])

}

getAxisAngle <- function(model, time) {

    indices <- seq(5,length(model@mCells[[timeToRow(model,time)]]),6)
    return(model@mCells[[timeToRow(model,time)]][indices])

}

getGrowthRates <- function(model, time) {

    indices <- seq(6,length(model@mCells[[timeToRow(model,time)]]),6)
    return(model@mCells[[timeToRow(model,time)]][indices])

}

#' \code{getCycleLengths} return the cycle lengths of each cells at time
#'
#' @param model a CellModel object
#' @param time time of interest
#' @return the cycle lengths of each cell at time
#' @export
#'

getCycleLengths <- function(model, time) {

    gr_rates <- getGrowthRates(model, time)
    return (1 + 2 * (sqrt(2) - 1) * .timeIncrement(model) * .nG(model) / gr_rates)

}

#' \code{show} display summary of CellModel class
#'
#' @param object A CellModel Object
#' @examples show(runCancerSim(1,6))
#' @rdname CellModel-class
#' @export
#'

setMethod("show", "CellModel",

    function (object) {

        cat("model parameters:\n")
        print(getParameters(object))
        cat("available functions:\n")
        cat("interactivePlot\n")
        cat("getCycleLengths\n")

    }

)

#' \code{getParameters} get a named list of parameters in the model
#'
#' @param model A CellModel
#' @param fullDist whether or not to return full distribution of cycle length
#' @return a named list of parameters in the model
#' @examples
#' getParameters(runCancerSim(1,1))
#' @rdname CellModel-class
#' @export
#'

getParameters <- function(model, fullDist=FALSE) {

    retDist = mean(.cycleLengthDist(model))

    if (fullDist) {

     retDist = .cycleLengthDist(model)

    }

    ret_val = list(

        initialNum = .initialNumCells(model),
        runTime = .runTime(model),           
        initialDensity = .initialDensity(model),           
        inheritGrowth = .inheritGrowth(model),           
        outputIncrement = .outputIncrement(model),           
        randSeed = .randSeed(model),           
        epsilon = .epsilon(model),           
        nG = .nG(model),
        timeIncrement = .timeIncrement(model),
        cycleLengthDist = retDist

    )           

    return(ret_val)

}


