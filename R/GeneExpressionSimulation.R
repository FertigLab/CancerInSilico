#' simulate gene expression data
#' @export
#'
#' @description simulate gene expression data for a set of pathways, using
#' the behavior of a CellModel as the basis for the simulation
#' @param model a CellModel object
#' @param pathways list of genes pathways
#' @param params GeneExpressionParams object
#' @return list of pathway activity and gene expression
inSilicoGeneExpression <- function(model, pathways,
params=new('GeneExpressionParams'))
{
    # calculate activity in each pathway at every sample point
    pwyActivity <- lapply(pathways, function(p)
        simulatePathwayActivity(p, model, params@sampleFreq,
        params@nCells, params@singleCell))

    # get mean expression matrix
    meanExp <- getMeanExp(pathways, pwyActivity)

    # call correct error model
    if (params@RNAseq & params@singleCell)
        exp <- simulateWithSplatter(meanExp, params)
    else if (params@RNAseq & !params@singleCell)
        exp <- simulateWithLimmaVoom(meanExp, params)
    else if (!params@RNAseq & !params@singleCell)
        exp <- simulateMicroArray(meanExp, params)

    # pad expression matrix with dummy genes
    if (params@nDummyGenes > 0)
    {
        dummyExp <- t(sapply(1:params@nDummyGenes, function(x)
            params@dummyDist(ncol(exp))))
        rownames(dummyExp) <- paste('dummy', 1:nrow(dummyExp), sep='_')
        colnames(dummyExp) <- colnames(exp)
        exp <- rbind(exp, dummyExp)
    }

    # return raw pathway activity as well as expression matrix
    return(list(pathways=pwyActivity, expression=exp))
}

#' verify gene expression data set is valid for this package
#' @export
#'
#' @description checks a data set before it is used to calibrate the 
#' pathway values for min/max expression
#' @param dataSet matrix of gene expression data where row names are genes
#' @param genes names of all genes being simulated
#' @return no value is return, but errors/warnings are thrown related to
#' potential problems in the data set
#' @importFrom stats median
#' @examples
#' data(referenceGeneExpression)
checkDataSet <- function(dataSet, genes)
{
    if (length(setdiff(genes, row.names(dataSet))))
        stop('dataSet does not contain all neccesary genes')
    else if (sum(is.na(unname(apply(dataSet[genes,], 2, median)))) > 0)
        stop('dataSet has NA values for pathway genes') #TODO: warning?
    else if (min(dataSet[!is.na(dataSet)]) < 0)
        stop('dataSet should be strictly non-negative')
}

#' Combine Gene Expression Matrices
#' @keywords internal
#'
#' @description combines mutliple matrices by in a similiar way to rbind,
#' but combines muttiple values for a single gene by a given function
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
    allGenes <- unique(unlist(lapply(expList, rownames)))
    output <- matrix(nrow = length(allGenes), ncol = nCols[1]) 
    rownames(output) <- allGenes
    colnames(output) <- colnames(expList[[1]])

    # combine info from all matrices containing the gene
    for (g in allGenes)
    {
        valid <- sapply(expList, function(m) g %in% rownames(m))
        exp <- sapply(expList[valid], function(m) m[g,])
        output[g,] <- apply(exp, 1, combineFUN)
    }   
    return(output)
}

#' Calculate Mean Expression Value
#' @keywords internal
#'
#' @description calculates the mean expression values for each gene and sample
#' by combining the pathway activity with the gene expression range
#' @param pathways list of pathway objects
#' @param activity list of pathway activity
#' @return matrix of mean expression values
getMeanExp <- function(pathways, activity)
{
    exp <- list()
    for (p in 1:length(pathways))
    {
        pwy <- pathways[[p]]
        mat <- (pwy@maxExpression - pwy@minExpression) %*% t(activity[[p]])
            + pwy@minExpression
        rownames(mat) <- pwy@genes
        colnames(mat) <- names(activity[[p]])
        exp[[p]] <- mat
    }
    return(combineGeneExpression(exp))
}

#' Bulk MicroArray Error Model
#' @keywords internal
#'
#' @description adds a normally distributed error across the entire expression
#' matrix
#' @param meanExp matrix of mean expression values
#' @param params GenExpressionParams object
#' @return matrix with error model incorporated
#' @importFrom stats rnorm
simulateMicroArray <- function(meanExp, params)
{
    error <- matrix(rnorm(length(meanExp)), ncol=ncol(meanExp))
    exp <- meanExp + pmax(params@perError * meanExp, params@perError) * error
    return(pmax(exp, 0))
}

#' Error Model found in Limma-Voom
#' @keywords internal
#'
#' @description Core error model presented in original voom paper, used for 
#' multiple RNA-seq simulations
#' @param mu mean expression matrix
#' @param bcvCommon biological variation parameter
#' @param bcvDF degrees of freedom for the chisq distribution
#' @return list of the cell means and the simulated counts
#' @importFrom stats rpois rgamma rchisq
voomErrorModel <- function(mu, bcvCommon, bcvDF)
{
    mu <- floor(mu) + 1
    nGenes <- nrow(mu)
    bcv <- (bcvCommon + (1 / sqrt(mu))) * sqrt(bcvDF / rchisq(nGenes, df=bcvDF))
    shape <- 1 / bcv^2
    scale <- mu / shape
    cellMeans <- matrix(rgamma(length(mu), shape=shape, scale=scale),
        nrow=nGenes)
    trueCounts <- matrix(rpois(length(mu), lambda=cellMeans), nrow=nGenes)
    rownames(trueCounts) <- rownames(mu)
    colnames(trueCounts) <- colnames(mu)
    return(list("cellMeans"=cellMeans, "trueCounts"=floor(pmax(trueCounts, 0))))
}

#' Bulk RNA-seq Error Model
#' @keywords internal
#'
#' @description adds the limma-voom error model to the data
#' @param meanExp matrix of mean expression values
#' @param params GenExpressionParams object
#' @return matrix with error model incorporated
simulateWithLimmaVoom <- function(meanExp, params)
{
    voomErrorModel(meanExp, params@bcvCommon, params@bcvDF)$trueCounts
}

#' Single Cell RNA-seq Error Mpdel
#' @keywords internal
#'
#' @description adds both the limma-voom error model and calculates dropout
#' @param meanExp matrix of mean expression values
#' @param params GenExpressionParams object
#' @return matrix with error model incorporated
#' @importFrom stats rbinom
simulateWithSplatter <- function(meanExp, params)
{
    exp <- voomErrorModel(meanExp, params@bcvCommon, params@bcvDF)
    counts <- exp$trueCounts
    if (params@dropoutPresent)
    {
        # Generate probabilites based on expression
        logistic <- function(x, x0, k) 1 / (1 + exp(-k * (x - x0)))
        prob <- sapply(1:ncol(exp$cellMeans), function(n)
        {
            eta <- log(exp$cellMeans[,n])
            return(logistic(eta, x0=params@dropoutMid, k=params@dropoutShape))
        })
        
        # Decide which counts to keep
        keep <- matrix(rbinom(length(counts), 1, 1 - prob), nrow=nrow(counts))
        counts <- counts * keep
    }
    return(counts)
}

