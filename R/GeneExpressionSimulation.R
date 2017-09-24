#' simulate gene expression data
#' @export
#'
#' @description simulate gene expression data for a set of pathways, using
#'  the behavior of a CellModel as the basis for the simulation
#' @param model a CellModel object
#' @param pathways list of genes pathways
#' @return list of pathway activity and gene expression
inSilicoGeneExpression <- function(model, pathways,
params=new('GeneExpressionParams'))
{
    # calculate activity in each pathway at every sample point
    pwyActivity <- lapply(pathways, function(p)
        simulatePathwayActivity(p, model, params@sampleFreq,
        params@nCells, params@singleCell))

    # simulate gene expression
    if (params@RNAseq & params@singleCell)
        exp <- simulateWithSplatter(pathways, pwyActivity, params)
    else if (params@RNAseq & !params@singleCell & !length(params@fasta))
        exp <- simulateWithLimmaVoom(pathways, pwyActivity, params)
    else if (params@RNAseq & !params@singleCell & length(params@fasta))
        exp <- simulateWithPolyester(pathways, pwyActivity, params)
    else
        exp <- simulateMicroArray(pathways, pwyActivity, params)

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
    rownames(output) <- allGenes
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

simulateWithSplatter <- function(pathways, activity, params)
{
    # generate mean expression
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
    exp <- floor(combineGeneExpression(exp)) + 1

    # splatter parameters
    splatParams <- params@splatParams
    splatParams <- splatter:::expandParams(splatParams)
    validObject(splatParams)
    seed <- splatter:::getParam(splatParams, 'seed')
    bcv.common <- splatter:::getParam(splatParams, "bcv.common")
    bcv.df <- splatter:::getParam(splatParams, "bcv.df")
    dropout.present <- splatter:::getParam(splatParams, "dropout.present")
    nGenes <- nrow(exp)
    nCells <- ncol(exp)
    cell.names <- colnames(exp)
    gene.names <- rownames(exp)

    set.seed(seed)

    # create SCEset
    phenos <- new('AnnotatedDataFrame', data = data.frame(Cell = colnames(exp)))
    rownames(phenos) <- cell.names
    features <- new('AnnotatedDataFrame', data = data.frame(Gene = rownames(exp)))
    rownames(features) <- gene.names
    sim <- scater:::newSCESet(countData = exp, phenoData = phenos, featureData = features)

    scater:::set_exprs(sim, 'BaseCellMeans') <- exp

    # simulate BCV
    base.means.cell <- scater:::get_exprs(sim, "BaseCellMeans")
    if (is.finite(bcv.df))
    {
        bcv <- (bcv.common + (1 / sqrt(base.means.cell))) *
            sqrt(bcv.df / rchisq(nGenes, df = bcv.df))
    }
    else
    {
        warning("'bcv.df' is infinite. This parameter will be ignored.")
        bcv <- (bcv.common + (1 / sqrt(base.means.cell)))
    }

    means.cell <- matrix(rgamma(nGenes * nCells, shape = 1 / (bcv ^ 2),
        scale = base.means.cell * (bcv ^ 2)), nrow = nGenes, ncol = nCells)
    colnames(means.cell) <- cell.names
    rownames(means.cell) <- gene.names

    scater:::set_exprs(sim, "BCV") <- bcv
    scater:::set_exprs(sim, "CellMeans") <- means.cell

    # simulate true counts
    cell.means <- scater:::get_exprs(sim, "CellMeans")
    true.counts <- matrix(rpois(nGenes * nCells, lambda = cell.means),
      nrow = nGenes, ncol = nCells)

    colnames(true.counts) <- cell.names
    rownames(true.counts) <- gene.names

    scater:::set_exprs(sim, "TrueCounts") <- true.counts

    # simulate dropout
    if (dropout.present) 
    {
        dropout.mid <- splatter:::getParam(splatParams, "dropout.mid")
        dropout.shape <- splatter:::getParam(splatParams, "dropout.shape")

        # Generate probabilites based on expression
        logistic <- function(x, x0, k) 1 / (1 + exp(-k * (x - x0)))
        drop.prob <- sapply(seq_len(nCells), function(idx)
        {
            eta <- log(cell.means[, idx])
            return(logistic(eta, x0 = dropout.mid, k = dropout.shape))
        })

        # Decide which counts to keep
        keep <- matrix(rbinom(nCells * nGenes, 1, 1 - drop.prob),
           nrow = nGenes, ncol = nCells)

        counts <- true.counts * keep

        colnames(drop.prob) <- cell.names
        rownames(drop.prob) <- gene.names
        colnames(keep) <- cell.names
        rownames(keep) <- gene.names

        scater:::set_exprs(sim, "DropProb") <- drop.prob
        scater:::set_exprs(sim, "Dropout") <- !keep
    }
    else
    {
        counts <- true.counts
    }
    scater::counts(sim) <- counts

    # Create new SCESet to make sure values are calculated correctly
    sce <- scater::newSCESet(countData = scater::counts(sim),
                     phenoData = new("AnnotatedDataFrame", data = Biobase:::pData(sim)),
                     featureData = new("AnnotatedDataFrame", data = Biobase:::fData(sim)))

    for (assay.name in names(Biobase:::assayData(sim))) {
        if (!(assay.name %in% names(Biobase:::assayData(sce)))) {
            scater:::set_exprs(sce, assay.name) <- scater:::get_exprs(sim, assay.name)
        }
    }
    return(scater::counts(sce))
}

oldSimulateWithSplatter <- function(pathways, activity, params)
{
    # common splatter parameters
    splatParams <- params@splatParams
    splatParams <- splatter::setParam(splatParams, "groupCells", 500)
    splatParams <- splatter::setParam(splatParams, "path.nonlinearProb", 0)
    splatParams <- splatter::setParam(splatParams, "path.skew", 0.5)
    splatParams <- splatter::setParam(splatParams, "de.prob", 1)
    splatParams <- splatter::setParam(splatParams, "de.downProb", 0)
    
    # simulate expression for each pathway
    exp <- list()
    for (p in 1:length(pathways))
    {
        # create dummy matrix of gene expression ranges
        totalGenes <- 10000
        rep <- max(floor(totalGenes / length(pathways[[p]]@genes)), 1)
        dummy <- matrix(nrow=(rep * length(pathways[[p]]@genes)),
            ncol=length(activity[[p]]))
        for (g in 1:length(pathways[[p]]@genes))
        {
            dummy[(rep*(g-1)+1):(rep*g),] <- floor(runif(ncol(dummy) * rep,
                pathways[[p]]@minExpression[g], pathways[[p]]@maxExpression[g]))
        }

        # normalize dummy matrix to library size
        dummyLibSizes <- colSums(dummy)
        dummyLibMed <- median(dummyLibSizes)
        dummyNormCounts <- t(t(dummy) / dummyLibSizes * dummyLibMed)
        dummyNormCounts <- dummyNormCounts[rowSums(dummyNormCounts > 0) > 1, ]

        # combine given parameters with mean and lib size estimates
        splatParams <- splatter::setParam(splatParams, "nGenes", nrow(dummy))
        splatParams <- splatter:::splatEstMean(dummyNormCounts, splatParams)
        splatParams <- splatter:::splatEstLib(dummy, splatParams)

        # simulate reference data
        refData <- splatter::splatSimulate(splatParams, method='paths',
            verbose=FALSE)

        # match cells to reference data        
        matchedCols <- c()
        for (i in 1:length(activity[[p]]))
        {
            diff <- abs(100 * activity[[p]][i] - refData@phenoData@data$Step)
            candidates <- which(diff < 5)
            if (length(candidates) > 0)
                matchedCols[i] <- sample(candidates, 1)
            else
                matchedCols[i] <- which.min(diff)
        } 

        # match genes to reference data
        matchedRows <- c()
        for (i in 1:length(pathways[[p]]@genes))
        {
            meanExp <- (pathways[[p]]@maxExpression[i] - pathways[[p]]@minExpression[i]) / 2
            diff <- abs(refData@featureData@data$BaseGeneMean - meanExp)
            matchedRows[i] <- which.min(diff)
        }

        # return final expression matrix
        exp[[p]] <- scater::counts(refData)[matchedRows, matchedCols]
        rownames(exp[[p]]) <- pathways[[p]]@genes
        colnames(exp[[p]]) <- names(activity[[p]])
    }

    # combine expression from all pathways
    return(combineGeneExpression(exp, params@combineFUN))
}

simulateWithPolyester <- function(pathways, activity, params)
{
    stop('not implemented')
}

simulateWithLimmaVoom <- function(pathways, activity, params)
{
    simError <- function(mu)
    {
        nGenes <- length(mu)
        mu.ref <- pmax(mu,1)
        BCV0 <- 0.2 + 1 / sqrt(mu.ref)
        df.BCV <- 40
        BCV <- BCV0 * sqrt(df.BCV / rchisq(nGenes, df=df.BCV))
        shape <- 1 / BCV^2
        scale <- mu / shape
        rpois(nGenes, lambda=rgamma(nGenes, shape=shape, scale=scale))
    }

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
    meanExp <- combineGeneExpression(exp)
    meanExp[] <- floor(pmax(apply(meanExp, 2, simError), 0))
    return(meanExp)
}

simulateMicroArray <- function(pathways, activity, params)
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

    meanExp <- combineGeneExpression(exp)
    error <- matrix(rnorm(length(meanExp)), ncol=ncol(meanExp))
    exp <- meanExp + pmax(params@perError * meanExp, params@perError) * error
    return(pmax(exp, 0))
}


