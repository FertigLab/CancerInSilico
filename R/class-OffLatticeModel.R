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
setClass('OffLatticeModel', contains = 'CellModel', slots = c(
    maxDeformation = 'numeric',
    maxTranslation = 'numeric',
    maxRotation = 'numeric'
))

##################### Generics ###################

##################### Methods ####################

setMethod('processParameters',
    signature('OffLatticeModel'),
    function(model, ...)
    {
        # for readability
        p <- model@params

        # call function on base class
        base <- new('CellModel', params = p)
        p <- processParameters(base)

        # get extra parameters
        if (is.null(p[['maxDeform']]))
            {p[['maxDeform']] <- list(...)$maxDeform}
        if (is.null(p[['maxTranslation']]))
            {p[['maxTranslation']] <- list(...)$maxTranslation}
        if (is.null(p[['maxRotate']]))
            {p[['maxRotate']] <- list(...)$maxRotate}
        
        # check off lattice parameters
        if (is.null(p[['maxDeform']])) stop('missing maxDeform')
        if (is.null(p[['maxTranslation']])) stop('missing maxTranslation')
        if (is.null(p[['maxRotate']])) stop('missing maxRotate')

        if (p[['maxDeform']] <= 0) stop('invalid maxDeform')
        if (p[['maxTranslation']] <= 0) stop('invalid maxTranslation')
        if (p[['maxRotate']] <= 0) stop('invalid maxRotate')

        # return parameters
        return (p)
    }
)

setMethod('timeToRow',
    signature('OffLatticeModel'),
    function(model, time)
    {
        maxTime <- length(model@cells) * model@params[['runTime']]
        if (time >= maxTime)
        {
            return (length(model@cells))
        }
        else
        {
            return (ceiling(time / model@params[['recordIncrement']]))
        }
    }
)

setMethod('getColumn',
    signature('OffLatticeModel'),
    function(model, time, col)
    {

    }
)

setMethod('getCellPhase',
    signature('OffLatticeModel'),
    function(model, id)
    {

    }
)

setMethod('getCoordinates',
    signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getRadii',
    signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getAxisLength',
    signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getAxisAngle',
    signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getGrowthRates',
    signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getCellTypes',
    signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getCycleLengths',
    signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getNumberOfCells',
    signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('getNumberOfNeighbors',
    signature('OffLatticeModel'),
    function(model, time, index, radius)
    {

    }
)

setMethod('getDensity',
    signature('OffLatticeModel'),
    function(model, time)
    {

    }
)

setMethod('plotCells',
    signature('OffLatticeModel'),
    function(model, time)
    {

    }
)
