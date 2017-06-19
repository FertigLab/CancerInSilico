#' Simulate Gene Expression Data
#' @export
#'
#' @description simulate gene expression data for a set of pathways, using
#'  the behavior of a CellModel as the basis for the simulation
#' @param model a CellModel object
#' @param pathways list of genes pathways
#' @param sampFreq how often to generate data
#' @param nGenes total number of genes used (matrix padded with dummys)
#' @param combineFUN function used to combine gene expression data
#' @param singleCell logical: simulate single cell data
#' @param nCells number of cells to use for single cell data
#' @param perError percent error?
#' @param microArray true if micro array data
#' @param randSeed random seed for simulation
#' @return matrix of gene expression data
inSilicoGeneExpression <- function(model, pathways, sampFreq=1,
nGenes=NULL, combineFUN=max, singleCell=FALSE, nCells=96, perError=0.1,
microArray=FALSE, randSeed=0, dataSet=NULL)
{
    # run simulation for each pathway
    pathwayOutput <- lapply(pathways, function(p)
        simulatePathwayExpression(p, model, sampFreq, singleCell, nCells))

    # combine expression matrices, add dummy genes, and shuffle order
    meanExp <- combineGeneExpression(pathwayOutput, combineFUN)
    if (!is.null(nGenes)) meanExp <- padExpMatrix(meanExp, nGenes)
    meanExp <- simulateError(meanExp, dataSet, perError, microArray)

    # shuffle order of genes
    return (meanExp[sample(nrow(meanExp)),])
}

#' Verify Gene Expression Data Set
#' @export
#'
#' @description Checks a data set before it is used to calibrate the 
#'  pathway values for min/max expression
#' @param dataSet matrix of gene expression data where row names are genes
#' @param genes names of all genes being simulated
#' @return no value is return, but errors/warnings are thrown related to
#'  potential problems in the data set
checkDataSet <- function(dataSet, genes)
{
    if (length(setdiff(genes, row.names(dataSet))))
        stop('dataSet does not contain all neccesary genes')
    else if (sum(is.na(unname(apply(dataSet[genes,], 2, median)))) > 0)
        stop('dataSet has NA values for pathway genes') #TODO: warning?
    else if (max(dataSet[!is.na(dataSet)]) > 50)
        warning('check that the data is log transformed (large max value)')
    else if (min(dataSet[!is.na(dataSet)]) < 0)
        stop('dataSet should be strictly non-negative')
}

#' Combine Gene Expression Matrices
#' @keywords internal
#'
#' @description combines mutliple matrices by in a similiar way to rbind,
#'  but combines muttiple values for a single gene by a given function
#' @param expList list of expression matrices
#' @param combineFUN function to use when combining multiple values
#' @return matrix containing combined expression values
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
        valid <- sapply(expList, function(m) g %in% row.names(m))
        exp <- sapply(expList[valid], function(m) m[g,])
        output[g,] <- apply(exp, 1, combineFUN)
    }   
    return(output)
}

#' Add Noise to Gene Expression Matrices
#' @export
#'
#' @description add genes with random expression values to the matrix
#' @param mat matrix of gene expression values
#' @param nGenes final number of genes in the expression matrix
#' @param dist distribution to sample random values from
#' @return gene expression matrix with dummy genes added in
padExpressionMatrix <- function(mat, nGenes, distr)
{
    if (nGenes <= nrow(mat)) return(mat)
    dummyExp <- matrix(nrow=nGenes - nrow(mat), ncol=ncol(mat))
    dummyExp[] <- sapply(1:length(dummyExp), function(x) distr())
    rownames(dummyExp) <- paste('dummy', 1:nrow(dummyExp), sep='_')
    colnames(dummyExp) <- colnames(mat)
    return(combineGeneExpression(list(mat, dummyExp)))
}

#' Add Simulated Error to Expression Data
#' @export
#'
#' @description add noise to all values in expression matrix
#' @param meanExp matrix of gene expression data
#' @param dataSet matrix of gene expression data where row names are genes
#' @param perError TODO
#' @param microArray whether this data is RNA-seq or microarray
#' @return gene expression matrix with error
simulateError <- function(meanExp, dataSet=NULL, perError, microArray)
{
    if (microArray)
    {
        normalError <- matrix(rnorm(length(meanExp)), nrow=nrow(meanExp))
        meanExp <- meanExp + pmax(perError*meanExp, perError) * normalError
        return(pmax(meanExp, 0))
    }
    else
    {
        meanExp <- round(2 ^ meanExp - 1) # convert to counts
        return(negBinError(meanExp, dataSet))
    }
}

#' Negative Binomial Error model
#' @keywords internal
#'
#' @description #TODO
#' @param meanExp matrix of mean expression data
#' @param dataSet
#' @param invChisq
negBinError <- function(meanExp, dataSet=NULL, invChisq=TRUE)
{
    # get number of genes and mean expression for each one
    if (is.null(dataSet))
    {
        referenceMean <- pmax(rowMeans(meanExp),1)
    }
    else
    {
        checkDataSet(dataSet, row.names(meanExp))
        dataSet <- dataSet[row.names(meanExp),]
        referenceMean <- pmax(rowMeans(round(2 ^ dataSet - 1)),1)
    }
    nGenes <- length(referenceMean)

    # Biological variation
    BCV0 <- 0.2 + 1 / sqrt(referenceMean)
    if (invChisq)
        BCV <- BCV0 * sqrt(40 / rchisq(nGenes, df=40))
    else
        BCV <- BCV0 * exp(rnorm(nGenes, mean=0, sd=0.25) / 2)

    shape <- 1 / BCV^2
    scale <- meanExp / shape
    lambda <- matrix(rgamma(nGenes, shape=shape, scale=scale), nGenes)

    # Technical variation
    return(matrix(rpois(nGenes, lambda=lambda), nGenes))
}
