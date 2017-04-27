library(methods)

#' @title CellType
#' @description The properties of a cell type
#'
#' @slot name the name of the cell type
#' @slot size the relative size (volume) of the cell
#' @slot minCycle minimum possible cell cycle length
#' @slot cycleLength function that returns sample from distribution of 
#'  cycle lengths
#' @export
setClass('CellType', slots = c(
    name = 'character',
    size = 'numeric',
    minCycle = 'numeric',
    cycleLength = 'function'
))

setMethod('initialize', 'CellType', 
    function(.Object, ...)
    {
         if (is.null(list(...)$name)) stop('missing name')
        .Object@size <- 1
        .Object@minCycle <- 48
        .Object@cycleLength <- function() {return(48)}
        .Object <- callNextMethod(.Object, ...)
        .Object
    }
)

setValidity('CellType',
    function(object)
    {
        if (length(object@name) == 0)
            "missing 'name'"
        else if (length(object@size) == 0)
            "missing 'size'"
        else if (length(object@minCycle) == 0)
            "missing 'minCycle'"
        else if (is.null(object@cycleLength()))   
            "missing 'cycleLength'"
        else if (object@size < 1)
            "'size' cannot be less than 1"
        else if (object@minCycle < 4)
            "'minCycle' cannot be less than 4"
        else
            TRUE
    }
)

