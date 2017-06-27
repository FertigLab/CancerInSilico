library(methods)

#' @include class-CellModel.R
NULL

################ Class Definition ################

#' Pathway Class
#' @export
#' @name Pathway-class
#' @rdname Pathway-class
#'
#' @description Describes the basic properties of a gene pathway
#' @details expressionScale is a function that accepts three arguments:
#'  model, cell, and time. It should return a number in [0,1] that describes
#'  how active the genes are in this pathway for a given cell in the model
#'  at a given time. In bulk data, the pathway activity is averaged
#'  and transformed by 1 / (1 + exp(-k * (x - M))) where 
#'  k = transformSlope and M = transformMidpoint.
#'  The scale determines how expressed genes in this pathway are. i.e.
#'  and scale of 0 means all genes will have minExpression value and a scale
#'  of 1 means all genes will have maxExpression value. In between these
#'  values the gene expression scales linearly.
#' @slot genes names of genes in the pathway
#' @slot expressionScale function descibing how this pathway is affected
#'  by the state of the model
#' @slot minExpression minimum expression value for each gene (vector)
#' @slot maxExpression maximum expression value for each gene (vector)
#' @slot transformSlope parameter for transforming bulk data
#' @slot transformMidpoint parameter for transforming bulk data
setClass('Pathway', slots = c(
    genes = 'character',
    expressionScale = 'function',
    minExpression = 'numeric',
    maxExpression = 'numeric',
    transformSlope = 'numeric',
    transformMidpoint = 'numeric'
))

#' Pathway Class Constructor
#' @rdname Pathway-class 
setMethod('initialize', 'Pathway',
    function(.Object, ...)
    {
        if (is.null(list(...)$genes))
            stop('missing genes')
        if (is.null(list(...)$expressionScale))
            stop('missing expressionScale')

        .Object <- callNextMethod(.Object, ...)
        return(.Object)
    }
)

# valid object check for pathway class
setValidity('Pathway',
    function(object)
    {
        if (!length(object@genes))
            'no genes'
        if (length(object@genes) != length(object@minExpression))
            'minExpression length doesn\'t match number of genes'
        if (length(object@genes) != length(object@maxExpression))
            'maxExpression length doesn\'t match number of genes'
        if (length(object@transformSlope) & !length(object@transformMidpoint)
        |  !length(object@transformSlope) & length(object@transformMidpoint))
            'need both/neither of transformSlope and transformMidpoint'
    }
)

##################### Methods ####################

#' simulate pathway activity based in model behavior
#' @keywords internal
#'
#' @description simulates pathway activity over the course of a cell model
#' @param pathway a 'Pathway' object
#' @param model a 'CellModel' object
#' @param sampFreq time in between samples
#' @param singleCell simulate single cell activity
#' @param sampleSize number of cells to sample in the case of single cell
#' @return vector with values in [0,1] indicatin pathway activity
simulatePathwayActivity <- function(pathway, model, sampFreq,
singleCell=FALSE, sampSize=0)
{
    # check parameters
    if (sampFreq <= 0) stop('invalid sampling frequency')
    if (!length(pathway@minExpression)) stop('pathway not calibrated')
    if (!length(pathway@maxExpression)) stop('pathway not calibrated')
    if (sampSize < 1) stop('invalid sample size')

    # find closest, valid, sampling frequency
    freq <- model@recordIncrement * ceiling(sampFreq / model@recordIncrement)

    # loop through each time
    scaleVector <- c()
    times <- seq(0, model@runTime, freq)
    for (t in times)
    {
        # determine which cells to calculate expression for
        cells <- sort(sample(1:getNumberOfCells(model, t), sampSize))

        # calculate scale and add to matrix
        scale <- sapply(cells, pathway@expressionScale, model=model, time=t)
        if (singleCell)
        {
            scaleVector <- c(scaleVector, scale)
        }
        else if (length(pathway@transformSlope))
        {
            logistic <- function(k, x0, x) 1 / (1 + exp(-k * (x - x0)))
            scaleVector <- c(scaleVector, logistic(pathway@transformSlope,
                pathway@transformMidpoint, mean(scale)))
        }
        else
        {
            scaleVector <- c(scaleVector, mean(scale))
        }
    }

    # add names of each cell/time combination
    if (!singleCell)
        names(scaleVector) <- paste('t', times, sep='_')
    else
        names(scaleVector) <- apply(expand.grid(1:sampSize,times), 1,
            function(r) paste('c', r[1], '_t', r[2], sep=''))

    # return vector of expression scale
    return(scaleVector)
}

#' simulate gene expression data for a single pathway
#' @keywords internal
#'
#' @description applies given pathway activity to the gene expression
#'  ranges to produce simulated gene expression data
#' @param pathway object of 'Pathway' class
#' @param activity named vector of pathway activity in [0,1]
#' @return scale applied to min/max expression values of the pathway
simulatePathwayExpression <- function(pathway, activity)
{
    gsMatrix <- matrix(0, length(pathway@genes), length(activity))
    rownames(gsMatrix) <- pathway@genes
    colnames(gsMatrix) <- names(activity)
    for (col in 1:length(activity))
    {
        gsMatrix[,col] <- (pathway@maxExpression - pathway@minExpression) *
            activity[col] + pathway@minExpression
    }
    return(gsMatrix)
}

#' calibrate pathway with data
#' @export
#'
#' @description sets the min and max values for each gene in a pathway
#'  based on a data set
#' @param pathway a 'Pathway' object
#' @param dataSet reference data set
#' @return pathway with min/max values for expression based on data set
calibratePathway <- function(pathway, dataSet)
{
    if (missing(dataSet)) stop('need a data set')
    
    checkDataSet(dataSet, pathway@genes)
    data <- dataSet[pathway@genes,]
    pathway@minExpression <- unname(apply(data, 1, min))
    pathway@maxExpression <- unname(apply(data, 1, max))
    return(pathway)
}
