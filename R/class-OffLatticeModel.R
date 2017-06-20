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
setGeneric('getCoordinates', function(model, time, cell)
    {standardGeneric('getCoordinates')})

#' @export
setGeneric('getRadius', function(model, time, cell)
    {standardGeneric('getRadius')})

#' @export
setGeneric('getAxisLength', function(model, time, cell)
    {standardGeneric('getAxisLength')})

#' @export
setGeneric('getAxisAngle', function(model, time, cell)
    {standardGeneric('getAxisAngle')})

#' @export
setGeneric('getGrowthAcceptRate', function(model, time, cell)
    {standardGeneric('getGrowthAcceptRate')})

##################### Methods ####################

setMethod('timeToRow', signature(model='OffLatticeModel'),
    function(model, time)
    {
        if (time > model@runTime | time < 0) stop('invalid time')
        else return (floor(time / model@recordIncrement) + 1)
    }
)

setMethod('getEntry', signature(model='OffLatticeModel'),
    function(model, time, cell, col)
    {
        row <- timeToRow(model, time)
        col <- col + 9 * (cell - 1)
        return(model@cells[[row]][col])
    }
)

setMethod('getCoordinates', signature(model='OffLatticeModel'),
    function(model, time, cell)
    {
        return(c(getEntry(model,time,cell,1), getEntry(model,time,cell,2)))
    }
)    

setMethod('getRadius', signature(model='OffLatticeModel'),
    function(model, time, cell)
    {
        return(getEntry(model, time, cell, 3))
    }
)

setMethod('getAxisLength', signature(model='OffLatticeModel'),
    function(model, time, cell)
    {
        return(getEntry(model, time, cell, 4))
    }
)

setMethod('getAxisAngle', signature(model='OffLatticeModel'),
    function(model, time, cell)
    {
        return(getEntry(model, time, cell, 5))
    }
)

setMethod('getCycleLength', signature(model='OffLatticeModel'),
    function(model, time, cell)
    {
        return(getEntry(model, time, cell, 6))
    }
)

setMethod('getCellPhase', signature(model='OffLatticeModel'),
    function(model, time, cell)
    {
        phases <- c('I', 'M', 'G0', 'G1', 'S', 'G2')
        return(phases[getEntry(model, time, cell, 7)+1])
    }
)

setMethod('getCellType', signature(model='OffLatticeModel'),
    function(model, time, cell)
    {
        return(getEntry(model, time, cell, 8) + 1)
    }
)

setMethod('getGrowthAcceptRate', signature('OffLatticeModel'),
    function(model, time, cell)
    {
        return(getEntry(model, time, cell, 9))
    }
)

setMethod('getNumberOfCells', signature('OffLatticeModel'),
    function(model, time)
    {
        return(length(model@cells[[timeToRow(model, time)]]) / 9)
    }
)

setMethod('getDensity', signature('OffLatticeModel'),
    function(model, time)
    {
        nCells <- getNumberOfCells(model, time)
        radii <- sapply(1:nCells, getRadius, model=model, time=time)
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
        # subset to nearby cells
        cells <- setdiff(1:getNumberOfCells(model, time), cell)
        cells <- cells[sapply(cells, function(c)
        {
            sizeA <- model@cellTypes[[getCellType(model, time, cell)]]@size
            sizeB <- model@cellTypes[[getCellType(model, time, c)]]@size
        
            posA <- getCoordinates(model, time, cell)
            posB <- getCoordinates(model, time, c)

            dist <- sqrt((posA[1] - posB[1])^2 + (posA[2] - posB[2])^2)
            return(dist - 2 * (sqrt(sizeA) + sqrt(sizeB)) < radius)
        })]

        # generate uniform grid of points within radius
        dist <- seq(0,1,0.01)
        ang <- dist * 2 * pi
        allCombo <- expand.grid(dist, ang)
        localGrid <- matrix(nrow=nrow(allCombo), ncol=2) 
        localGrid[,1] <- radius * sqrt(allCombo[,1]) * cos(allCombo[,2])
        localGrid[,2] <- radius * sqrt(allCombo[,1]) * sin(allCombo[,2])

        # return true if point inside cell
        inside <- function(c, p)
        {
            coords <- getCoordinates(model, time, c)
            rad <- getRadius(model, time, c)
            axisLen <- getAxisLength(model, time, c)
            axisAng <- getAxisAngle(model, time, c)
    
            x1 <- coords[1] + (0.5 * axisLen - rad) * cos(axisAng)  
            y1 <- coords[2] + (0.5 * axisLen - rad) * sin(axisAng)
            x2 <- coords[1] - (0.5 * axisLen - rad) * cos(axisAng)
            y2 <- coords[2] - (0.5 * axisLen - rad) * sin(axisAng)

            dist2 <- min((x1-p[1])^2+(y1-p[2])^2, (x2-p[1])^2+(y2-p[2])^2)
            return(dist2 < rad ^ 2)            
        }

        # calculate proportion of points inside of a cell
        # TODO: include points outside boundary in density
        numPoints <- sum(apply(localGrid, 1, function(p)
            sum(sapply(cells, inside, p=p)) > 0))
        numPointsBase <- sum(apply(localGrid, 1, function(p) inside(cell, p)))
        return(numPoints / (nrow(localGrid) - numPointsBase))
    }
)

setMethod('plotCells', signature('OffLatticeModel'),
    function(model, time)
    {
        # get all the cell information
        nCells <- getNumberOfCells(model, time)
        coords <- sapply(1:nCells, getCoordinates, model=model, time=time)
        radii <- sapply(1:nCells, getRadius, model=model, time=time)
        axisLen <- sapply(1:nCells, getAxisLength, model=model, time=time)
        axisAng <- sapply(1:nCells, getAxisAngle, model=model, time=time)
        phases <- sapply(1:nCells, getCellPhase, model=model, time=time)
        mitNdx <- rep(phases, 2) == 'M'

        # calculate plot bounds
        mn <- ifelse(model@boundary > 0, -model@boundary-2, min(coords)-2)
        mx <- ifelse(model@boundary > 0,  model@boundary+2, max(coords)+2)

        # create the plot template
        plot(c(mn, mx), c(mn, mx), main=paste("Plot of CellModel At Time",
            time), xlab="", ylab="", type="n", asp=1)
          
        # get all (x,y) pairs for each of the cell centers
        x_1 <- coords[1,] + (0.5 * axisLen - radii) * cos(axisAng)
        x_2 <- coords[1,] - (0.5 * axisLen - radii) * cos(axisAng)
        y_1 <- coords[2,] + (0.5 * axisLen - radii) * sin(axisAng)
        y_2 <- coords[2,] - (0.5 * axisLen - radii) * sin(axisAng)

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


