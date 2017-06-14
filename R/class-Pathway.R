#' @include class-CellModel.R
NULL

library(methods)

################ Class Definition ################

#' @title Pathway
#' @description Describes the basic properties of a gene pathway
#'
#' @slot genes names of genes in the pathway
#' @slot expressionScale function descibing how this pathway is affected
#'  by the state of the model
#' @slot minExpression minimum expression value for each gene (vector)
#' @slot maxExpression maximum expression value for each gene (vector)
#' @export
setClass('Pathway', slots = c(
    genes = 'character',
    expressionScale = 'function',
    minExpression = 'numeric',
    maxExpression = 'numeric'   
))

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

setValidity('Pathway',
    function(object)
    {
        if (!length(object@genes)) 'no genes'
    }
)

##################### Methods ####################

simulatePathwayExpression <- function(pathway, model, sampFreq,
singleCell=FALSE, sampSize=0)
{
    # check parameters
    if (sampFreq <= 0) stop('invalid sampling frequency')
    if (!singleCell) {sampSize <- 1}
    else if (sampSize < 1) stop('invalid sample size for single cell')

    # find closest, valid, sampling frequency
    sampFreq <- model@recordIncrement *
        ceiling(sampFreq / model@recordIncrement)

    # get number of time points and create the return matrix
    times <- seq(0, model@runTime, sampFreq)
    gsMatrix <- matrix(0,length(pathway@genes),length(times)*sampSize)
    rownames(gsMatrix) <- pathway@genes
    colnames(gsMatrix) <- ifelse(singleCell, 
        apply(expand.grid(1:sampSize,times), 1, function(r)
        paste('c',r[1],'_t',r[2],sep='')), paste('t',times,sep='_'))

    # loop through each time
    for (t in 1:length(times))
    {
        # determine which cells to calculate expression for
        cells <- 1:getNumberOfCells(model, times[t])
        if (singleCell) cells <- sort(sample(cells, sampSize))

        # calculate gene expression
        scale <- sapply(cells, function(c)
            pathway@expressionScale(model, c, times[t]))
        exp <- (pathway@maxExpression - pathway@minExpression) %*%
            t(scale) + pathway@minExpression

        # add expression to matrix
        cols <- (sampSize * (t-1) + 1):(sampSize * t)            
        if (singleCell) gsMatrix[,cols] <- exp
        else gsMatrix[,cols] <- rowMeans(exp)
    }
    return(gsMatrix)
}

calibratePathway <- function(pathway, dataSet=NULL, dist=NULL)
{
    if (is.null(dataSet)) stop('need a data set')
    
    checkDataSet(dataSet, pathway@genes)
    data <- dataSet[pathway@genes,]
    pathway@minExpression <- unname(apply(data, 1, min))
    pathway@maxExpression <- unname(apply(data, 1, max))
    return(pathway)
}
