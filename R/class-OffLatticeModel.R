#' @include class-CellModel.R
NULL

library(methods)

################ Class Definition ################

#' @title OffLatticeModel
#' @description General description of an off-lattice cell-based model.
#'  not quite a full implementation, but contains much of the neccesary
#'  structure for models of this type
#'
#' @slot maxDeformation the largest distance the axis of a cell can increase
#' @slot maxTranslation the largest distance the center of a cell can move
#' @slot maxRotation the largest angle a cell can rotate
#' @export
setClass('OffLatticeModel', contains = c('CellModel', 'VIRTUAL'), slots = c(
    maxDeformation = 'numeric',
    maxTranslation = 'numeric',
    maxRotation = 'numeric'
))

setMethod('initialize', 'OffLatticeModel',
    function(.Object, maxDeformation = 0.1, maxTranslation = 0.1,
    maxRotation = 0.3, ...)
    {
        # store parameters, don't overwrite existing value
        if (!length(.Object@maxDeformation))
            .Object@maxDeformation <- maxDeformation
        if (!length(.Object@maxTranslation))
            .Object@maxTranslation <- maxTranslation
        if (!length(.Object@maxRotation))
            .Object@maxRotation <- maxRotation

        # finish intialization, return object
        .Object <- callNextMethod(.Object, ...)
        return(.Object)
    }
)

setValidity('OffLatticeModel',
    function(object)
    {
        if (length(object@maxDeformation) == 0)
            "missing 'maxDeformation'"
        else if (length(object@maxTranslation) == 0)
            "missing 'maxTranslation'"
        else if (length(object@maxRotation) == 0)
            "missing 'maxRotation'"
        else if (object@maxDeformation <= 0)
            "'maxDeformation' must be greater than zero"
        else if (object@maxTranslation <= 0)
            "'maxTranslation' must be greater than zero"
        else if (object@maxRotation <= 0)
            "'maxRotation' must be greater than zero"
    }
)

##################### Generics ###################

setGeneric('getCoordinates', function(model, time)
    {standardGeneric('getCoordinates')})

setGeneric('getRadii', function(model, time)
    {standardGeneric('getRadii')})

setGeneric('getAxisLength', function(model, time)
    {standardGeneric('getAxisLength')})

setGeneric('getAxisAngle', function(model, time)
    {standardGeneric('getAxisAngle')})

setGeneric('getNumberOfNeighbors', function(model, time, cell, radius)
    {standardGeneric('getNumberOfNeighbors')})

##################### Methods ####################

setMethod('getCoordinates', signature('OffLatticeModel'),
    function(model, time)
    {

    }
)    

setMethod('getRadii', signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getAxisLength', signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getAxisAngle', signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getNumberOfNeighbors', signature('OffLatticeModel'),
    function(model, time, cell, radius)
    {
    
    }
)

setMethod('getCellPhases', signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getCellTypes', signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getCycleLengths', signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getNumberOfCells', signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getDensity', signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('timeToRow', signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getColumn', signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getColumn', signature('OffLatticeModel'),
    function(model, time, col)
    {

    }
)

setMethod('plotCells', signature('OffLatticeModel'),
    function(model, time)
    {
        
    }
)


