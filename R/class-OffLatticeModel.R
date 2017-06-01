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
setGeneric('getRadii', function(model, time)
    {standardGeneric('getRadii')})

#' @export
setGeneric('getAxisLength', function(model, time)
    {standardGeneric('getAxisLength')})

#' @export
setGeneric('getAxisAngle', function(model, time)
    {standardGeneric('getAxisAngle')})

#' @export
setGeneric('getContactInhibition', function(model, time)
    {standardGeneric('getContactInhibition')})

#' @export
setGeneric('getNumberOfNeighbors', function(model, time, cell, radius)
    {standardGeneric('getNumberOfNeighbors')})

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

setMethod('getRadii', signature('OffLatticeModel'),
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

setMethod('getContactInhibition', signature('OffLatticeModel'),
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

setMethod('getNumberOfNeighbors', signature('OffLatticeModel'),
    function(model, time, cell, radius)
    {
        num <- 0
        cds <- getCoordinates(model, time)
        for (i in setdiff(1:nrow(cds), cell))
        {
            dist2 <- (cds[i,1]-cds[cell,1])^2 + (cds[i,2]-cds[cell,2])^2
            if (dist2 < radius^2) num <- num + 1
        }
        return(num)
    }
)

setMethod('plotCells', signature('OffLatticeModel'),
    function(model, time)
    {
        # get all the cell information
        coords <- getCoordinates(model, time)
        radii <- getRadii(model, time)
        axis_len <- getAxisLength(model, time)
        axis_ang <- getAxisAngle(model, time)
        mitNdx <- rep(getCellPhases(model, time), 2) == 'M'

        # calculate plot bounds
        if (model@boundary > 0)
        {
            mn <- -model@boundary - 2
            mx <- model@boundary + 2
        }
        else
        {
            mn <- min(coords) - 2
            mx <- max(coords) + 2
        }

        # create the plot template
        plot(c(mn, mx), c(mn, mx), main = paste("Plot of CellModel At Time",
            time), xlab = "", ylab = "", type = "n", asp = 1)
          
        # get all (x,y) pairs for each of the cell centers
        x_1 <- coords[,1] + (0.5 * axis_len - radii) * cos(axis_ang)
        x_2 <- coords[,1] - (0.5 * axis_len - radii) * cos(axis_ang)
        y_1 <- coords[,2] + (0.5 * axis_len - radii) * sin(axis_ang)
        y_2 <- coords[,2] - (0.5 * axis_len - radii) * sin(axis_ang)

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


