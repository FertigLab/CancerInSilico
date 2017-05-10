#' \code{inSilicoGeneExpression}
#' @param model a CellModel object
#' @param pathways list of genes pathways
#' @param sampFreq how often to generate data
#' @param nGenes total number of genes used (matrix padded with dummys)
#' @param combineFUN function used to combine gene expression data
#' @param singleCell logical: simulate single cell data
#' @param nCells number of cells to use for single cell data
#' @param perError percent error?
#' @param microArray true if micro array data
#' @export
#'
inSilicoGeneExpression <- function(model, pathways, sampFreq=1, nGenes=NULL,
combineFUN=max, singleCell=FALSE, nCells=96, perError=0.1, microArray=FALSE)
{
    # run simulation for each pathway
    pathwayOutput <- lapply(pathways, function(p) simulateExpression(p,
        model, sampFreq, singleCell, nCells))

    # combine expression matrices, add dummy genes, and shuffle order
    meanExp <- combineGeneExpression(pathwayOutput, combineFUN))
    if (!is.null(nGenes)) meanExp <- padExpMatrix(meanExp, nGenes)
    meanExp <- meanExp[sample(nrow(meanExp)),]

    # add simulated error
    return (simulateError(meanExp, pathways, fasta, dataSet, perError,
        attrsep, microArray))
}

checkDataSet <- function(dataSet, genes)
{
    if (length(setdiff(genes, row.names(dataSet))))
        stop('dataSet does not contain all neccesary genes')
    else if (sum(is.na(unname(apply(dataSet[genes,], 2, median)))) > 0)
        stop('dataSet has NA values for pathway genes')
    else if (max(dataSet[!is.na(dataSet)]) > 50)
        warning('check that the data is log transformed (large max value)')
    else if (min(dataSet[!is.na(dataSet)]) < 0)
        stop('dataSet should be strictly non-negative')
}

calibratePathways <- function(pathways, dataSet=NULL, distFunc=NULL)
{
    if (is.null(dataSet) & is.null(distFunc))
        stop('need either a data set or a distribution to calibrate')

    for (p in names(pathways))
    {
        if (!is.null(dataSet))
        {
            checkDataSet(dataSet, pathways[[p]]@genes)
            data <- dataSet[pathways[[p]]@genes,]
            pathways[[p]]@minExpression <- unname(apply(data, 1, min))
            pathways[[p]]@maxExpression <- unname(apply(data, 1, max))
        }
        else if (!is.null(dist))
        {
            distribution <- distFunc(pathways[[p]])
            pathways[[p]]@minExpression <- distribution[['min']]
            pathways[[p]]@maxExpression <- distribution[['max']]
        }
    }
    return(pathways)
}

# combine gene expression matrices according to some function
combineGeneExpression <- function(expList, combineFUN=max)
{
    # verify same number of columns
    nCols <- sapply(expList, ncol)
    if (!all(nCols == nCols[1]))
        stop("unequal number of columns in gene expression matrices")

    # initalize matrix
    allGenes <- unique(unlist(lapply(expList, row.names)))
    output <- matrix(nrow = length(allGenes), ncol = nCols[1]) 
    row.names(output) <- allGenes
    colnames(output) <- colnames(expList[[1]])
 
    # combine info from all matrices containing the gene
    for (g in allGenes)
    {
        valid <- sapply(expLst, function(m) g %in% row.names(m))
        exp <- sapply(expList[valid], function(m) m[g,])
        output[g,] <- apply(exp, 1, combineFUN)
    }   
    return(output)
}

padExpMatrix <- function(meanExp, nGenes)
{
    dummyExp <- matrix(nrow=max(nGenes-nrow(meanExp),0), ncol=ncol(meanExp))
    dummyExp[] <- sample(meanExp, length(dummyExp), replace=TRUE)
    rownames(dummyExp) <- paste('dummy', 1:nrow(dummyExp), sep='_')
    colnames(dummyExp) <- colnames(meanExp)
    return(combineGeneExpression(list(meanExp, dummyExp)))
}

simulateError <- function(meanExp, dataSet, perError, microArray)
{
    if (microArray)
    {
        error <- pmax(perError * meanExp, perError) * 
            matrix(rnorm(length(meanExp)), nrow(meanExp), ncol(meanExp))
        output <- meanExp + error
        output[output < 0] <- 0
    }
    else
    {
        meanExp <- round(2 ^ meanExp - 1)
        output <- apply(meanExp, 2, function(exp) NBsim(exp, dataSet, TRUE))
    }
    dimnames(output) <- dimnames(meanExp)
    return(output)
}

#TODO: clarify this error model (limma-voom)
NBsim <- function(pwyMean, dataSet = NULL, invChisq = TRUE)
{
    # get number of genes and mean expression for each one
    if (!is.null(dataSet) & checkDataset(dataSet, row.names(pwyMean)))
        {referenceMean <- pmax(round(rowMeans(2 ^ dataSet - 1)),1)}
    else
        {referenceMean <- pmax(pwyMean,1)}
    nGenes <- length(referenceMean)

    # Biological variation
    BCV0 <- 0.2 + 1 / sqrt(referenceMean)
    if (invChisq)
        BCV <- BCV0 * sqrt(40 / rchisq(nGenes, df=40))
    else
        BCV <- BCV0 * exp(rnorm(nGenes, mean=0, sd=0.25) / 2)

    shape <- 1 / BCV^2
    scale <- pwyMean / shape
    lambda <- matrix(rgamma(nGenes, shape=shape, scale=scale), nGenes)

    # Technical variation
    return(matrix(rpois(ngenes, lambda=lambda), nGenes))
}

