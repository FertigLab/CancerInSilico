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
                        m_cells = "list",
                        m_initialNumCells = "numeric",
                        m_runTime = "numeric",
                        m_initialDensity = "numeric",
                        m_inheritGrowth = "logical",
                        m_outputIncrement = "numeric",
                        m_randSeed = "numeric",
                        m_epsilon = "numeric",
                        m_nG = "numeric",
                        m_timeIncrement = "numeric",
                        m_cycleLengthDist = "numeric" ))

## getters

.initialNumCells <- function(model) {

    return (model@m_initialNumCells)

}

.runTime <- function(model) {

    return (model@m_runTime)

}

.initialDensity <- function(model) {

    return (model@m_initialDensity)

}

.inheritGrowth <- function(model) {

    return (model@m_inheritGrowth)

}

.outputIncrement <- function(model) {

    return (model@m_outputIncrement)

}

.randSeed <- function(model) {

    return (model@m_randSeed)

}

.epsilon <- function(model) {

    return (model@m_epsilon)

}

.nG <- function(model) {

    return (model@m_nG)

}

.timeIncrement <- function(model) {

    return (model@m_timeIncrement)

}

.cycleLengthDist <- function(model) {

    return (model@m_cycleLengthDist)

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

        cat("model parameters:")
        print(getParameters(object))
        cat("available functions:\n")
        cat("plotInteractive\n")
        cat("plotCellsAtTime\n")
        cat("getParameters\n")
        cat("getDensity\n")
        cat("getNumberOfCells\n")
        cat("getCycleLengthDistribution\n")

    }

)

