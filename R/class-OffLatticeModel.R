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

#' @export
setGeneric('getCoordinates', function(model, time)
    {standardGeneric('getCoordinates')})

#' @export
setGeneric('getRadius', function(model, time)
    {standardGeneric('getRadius')})

#' @export
setGeneric('getAxisLength', function(model, time)
    {standardGeneric('getAxisLength')})

#' @export
setGeneric('getAxisAngle', function(model, time)
    {standardGeneric('getAxisAngle')})

#' @export
setGeneric('getGrowthAcceptRate', function(model, time)
    {standardGeneric('getGrowthAcceptRate')})

##################### Methods ####################

setMethod('timeToRow', signature('OffLatticeModel'),
    function(model, time)
    {
        if (time > model@runTime | time < 0) stop('invalid time')
        else return (floor(time / model@recordIncrement) + 1)
    }
)

setMethod('getColumn', signature('OffLatticeModel'),
    function(model, time, col)
    {
        row <- timeToRow(model, time)
        indices <- seq(col, length(model@cells[[row]]), 9)
        return (model@cells[[row]][indices])
    }
)

setMethod('getCoordinates', signature('OffLatticeModel'),
    function(model, time)
    {
        xCoord <- getColumn(model, time, 1)
        yCoord <- getColumn(model, time, 2)
        return(matrix(c(xCoord, yCoord), nrow = length(xCoord)))
    }
)    

setMethod('getRadius', signature('OffLatticeModel'),
    function(model, time)
    {
        return(getColumn(model, time, 3))
    }
)

setMethod('getAxisLength', signature('OffLatticeModel'),
    function(model, time)
    {
        return(getColumn(model, time, 4))
    }
)

setMethod('getAxisAngle', signature('OffLatticeModel'),
    function(model, time)
    {
        return(getColumn(model, time, 5))
    }
)

setMethod('getCycleLengths', signature('OffLatticeModel'),
    function(model, time)
    {
        return(getColumn(model, time, 6))
    }
)

setMethod('getCellPhases', signature('OffLatticeModel'),
    function(model, time)
    {
        indices <- getColumn(model, time, 7)
        phases <- c('I', 'M', 'G0', 'G1', 'S', 'G2')
        get_phase <- function(i) phases[i+1]
        return(sapply(indices, get_phase))
    }
)

setMethod('getCellTypes', signature('OffLatticeModel'),
    function(model, time)
    {
        return(getColumn(model, time, 8) + 1)
    }
)

setMethod('getGrowthAcceptRate', signature('OffLatticeModel'),
    function(model, time)
    {
        return(getColumn(model, time, 9))
    }
)

setMethod('getNumberOfCells', signature('OffLatticeModel'),
    function(model, time)
    {
        return(sum(getRadii(model, time) > 0))
    }
)

setMethod('getDensity', signature('OffLatticeModel'),
    function(model, time)
    {
        radii <- getRadii(model, time)
        if (model@boundary > 0)
        {
            return(sum(radii ** 2) / model@boundary ^ 2)
        }
        else
        {
            coords <- getCoordinates(model, time)
            d <- max(sqrt(coords[,1] ** 2 + coords[,2] ** 2) + radii)
            return(sum(radii ** 2) / (d ^ 2))
        }
    }
)

setMethod('getLocalDensity', signature('OffLatticeModel'),
    function(model, time, cell, radius)
    {
#        coords <- getCoordinates(model, time)
#        axisLen <- getAxisLength(model, time)
#        axisAng <- getAxisAngle(model, time)
#        rad <- getRadius(model, time)
#        types <- getCellTypes(model, time)
#        cellSize <- sapply(types, function(i) model@cellTypes[[i]]@size)
#
#        ndx <- setdiff(1:getNumberOfCells(model, time), cell)
#        x_1 <- coords[ndx,1] + (0.5 * axisLen - rad) * cos(axisAng)
#        x_2 <- coords[ndx,1] - (0.5 * axisLen - rad) * cos(axisAng)
#        y_1 <- coords[ndx,2] + (0.5 * axisLen - rad) * sin(axisAng)
#        y_2 <- coords[ndx,2] - (0.5 * axisLen - rad) * sin(axisAng)
#
#        dist <- function(a, b) (a[1] - b[1])^2 + (a[2] - b[2])^2
#
#        d1 <- apply(cbind(x_1, y_1), 1, dist, b=coords[cell,]) - rad[ndx]
#        d2 <- apply(cbind(x_2, y_2), 1, dist, b=coords[cell,]) - rad[ndx] - 
#        d <- c(d1, d2)
#        w <- c(rep(1, length(d1)), 1 - axisLen / sqrt(16 * cellSize))
    }
)

setMethod('plotCells', signature('OffLatticeModel'),
    function(model, time)
    {
        # get all the cell information
        coords <- getCoordinates(model, time)
        radii <- getRadii(model, time)
        axisLen <- getAxisLength(model, time)
        axisAng <- getAxisAngle(model, time)
        mitNdx <- rep(getCellPhases(model, time), 2) == 'M'

        # calculate plot bounds
        mn <- ifelse(model@boundary > 0, -model@boundary-2, min(coords)-2)
        mx <- ifelse(model@boundary > 0,  model@boundary+2, max(coords)+2)

        # create the plot template
        plot(c(mn, mx), c(mn, mx), main=paste("Plot of CellModel At Time",
            time), xlab="", ylab="", type="n", asp=1)
          
        # get all (x,y) pairs for each of the cell centers
        x_1 <- coords[,1] + (0.5 * axisLen - radii) * cos(axisAng)
        x_2 <- coords[,1] - (0.5 * axisLen - radii) * cos(axisAng)
        y_1 <- coords[,2] + (0.5 * axisLen - radii) * sin(axisAng)
        y_2 <- coords[,2] - (0.5 * axisLen - radii) * sin(axisAng)

        # combine all coordinate pairs along with the radii
        x <- c(x_1,x_2)
        y <- c(y_1,y_2)
        rad <- c(radii, radii)
    
        # plot the cells
        if (sum(mitNdx))
            symbols(x[mitNdx], y[mitNdx], circles=rad[mitNdx],
                inches=FALSE, add=TRUE, bg="black", fg="black")
        if (sum(!mitNdx))
            symbols(x[!mitNdx], y[!mitNdx], circles=rad[!mitNdx],
                inches=FALSE, add=TRUE, bg="bisque4", fg="bisque4")

        # draw boundary
        symbols(0, 0, circles = model@boundary, inches = FALSE, add = TRUE,
            lwd = 2)
    }
)


