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

setMethod('getCellDistance', signature(model='OffLatticeModel'),
    function(model, time, cellA, cellB)
    {
        centers <- function(model, time, cell)
        {
            crds <- getCoordinates(model, time, cell)
            rad <- getRadius(model, time, cell)
            axisLen <- getAxisLength(model, time, cell)
            axisAng <- getAxisAngle(model, time, cell)
            
            x1 <- crds[1] + (0.5 * axisLen - rad) * cos(axisAng)
            y1 <- crds[2] + (0.5 * axisLen - rad) * sin(axisAng)
            x2 <- crds[1] - (0.5 * axisLen - rad) * cos(axisAng)
            y2 <- crds[2] - (0.5 * axisLen - rad) * sin(axisAng)
            return(matrix(c(x1,x2,y1,y2), ncol=2))
        }

        cA <- centers(model, time, cellA)
        cB <- centers(model, time, cellB)

        minDist <- (cA[1,1]-cB[1,1])^2 + (cA[1,2]-cB[1,2])^2
        minDist <- min(minDist, (cA[1,1]-cB[2,1])^2 + (cA[1,2]-cB[2,2])^2)
        minDist <- min(minDist, (cA[2,1]-cB[1,1])^2 + (cA[2,2]-cB[1,2])^2)
        minDist <- min(minDist, (cA[2,1]-cB[1,1])^2 + (cA[2,2]-cB[1,2])^2)
        return(sqrt(minDist) - getRadius(model, time, cellA) - 
            getRadius(model, time, cellB))
    }
)

setMethod('getLocalDensity', signature('OffLatticeModel'),
    function(model, time, cell, radius)
    {
        dis <- function(a,b) sqrt((a[1]-b[1])^2 + (a[2]-b[2])^2)

        # generate grid around point
        genGrid <- function(p1, rad, p2=NULL)
        {
            width <- seq(-rad, rad, length.out=10)
            grid <- as.matrix(unname(expand.grid(width, width)))
            grid <- grid[apply(grid, 1, dis, b=c(0,0)) < rad,]
            grid <- t(t(grid) + p1)

            if (!is.null(p2))
                grid <- grid[apply(grid, 1, dis, b=p1) < apply(grid, 1, dis, b=p2),]
            return(grid)
        }

        # find nearby cells
        cells <- setdiff(1:getNumberOfCells(model, time), cell)
        cells <- cells[sapply(cells, function(c) getRadius(model, time, cell) +
            getCellDistance(model, time, cell, c) < radius)]
        if (!length(cells)) return(0)

        # get cell info
        coords <- sapply(cells, getCoordinates, model=model, time=time)
        rad <- sapply(cells, getRadius, model=model, time=time)
        axisLen <- sapply(cells, getAxisLength, model=model, time=time)
        axisAng <- sapply(cells, getAxisAngle, model=model, time=time)
        type <- sapply(cells, getCellType, model=model, time=time)
        sz <- sapply(type, function(t) model@cellTypes[[t]]@size)

        # find cell center coordinates
        term <- 0.5 * axisLen - rad
        p1 <- cbind(coords[1]+term*cos(axisAng), coords[2]+term*sin(axisAng))
        p2 <- cbind(coords[1]-term*cos(axisAng), coords[2]-term*sin(axisAng))

        grid <- matrix(nrow=0, ncol=2)
        for (c in 1:length(cells))
        {
            if (all.equal(2 * rad[c], axisLen[c], tol=1e-3) == TRUE)
                grid <- rbind(grid, genGrid(coords[,c], rad[c]))
            else
                grid <- rbind(grid, rbind(genGrid(p1[c,], rad[c], p2[c,]),
                    genGrid(p2[c,], rad[c], p1[c,])))
        }

        # check points for being in radius, return proportion of area
        cellCoords <- getCoordinates(model, time, cell)
        numPoints <- apply(grid, 1, function(p) dis(p, cellCoords) < radius)
        prop <- sum(numPoints) / nrow(grid)
        area <- sapply(1:length(cells), function(c) ifelse(all.equal(2*rad[c],
            axisLen[c], tol=1e-3)==TRUE, rad[c]^2, 2*sz[c]))
        return(prop * sum(area) / (radius^2 - getRadius(model, time, cell)^2))
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
        x1 <- coords[1,] + (0.5 * axisLen - radii) * cos(axisAng)
        x2 <- coords[1,] - (0.5 * axisLen - radii) * cos(axisAng)
        y1 <- coords[2,] + (0.5 * axisLen - radii) * sin(axisAng)
        y2 <- coords[2,] - (0.5 * axisLen - radii) * sin(axisAng)

        # combine all coordinate pairs along with the radii
        x <- c(x1,x2)
        y <- c(y1,y2)
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



