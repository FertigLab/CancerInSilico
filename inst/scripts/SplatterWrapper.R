library(splatter)
library(MASS)

# pathways - list of lists where there inner list is a named list (genes are
#   names) and the values of the list are mean expression,
# patterns - list of matrices where rows are cells and columns are timepoints,
#   values in the matrix are in [0,1] and represent the pathway activity of the
#   cell at each time point
# each pathway must have the same column names (cells) and may have overlapping
#   genes (row names)
# each pattern should have the same length, the name of each entry is the name
#   of the timepoint - these must also be identical across patterns
splatSimulateExtended <- function(pathways, patterns, combineFun=max)
{
    # combine genes appearing in mutliple pathways
    pathways <- removeGeneDuplicates(pathways, combineFun)

    # process pathway information
    dist <- suppressWarnings(fitdistr(unname(sapply(pathways[[1]], function(p) p[1])), 'gamma'))

    # generate reference gene expression
    params <- newSplatParams()
    params <- setParam(params, 'mean.shape', unname(dist$estimate[1]))
    params <- setParam(params, 'mean.rate', unname(dist$estimate[2]))
    params <- setParam(params, 'nGenes', length(pathways[[1]]))
    params <- setParam(params, 'groupCells', c(10))
    params <- setParam(params, 'out.prob', 0)
    params <- setParam(params, 'de.prob', 0)
    sim <- splatter::splatSimulate(params)
    referenceMatrix <- counts(sim)
   
    # sort reference matrix, first by mean expression, then by relative    
    baseMean <- sim@featureData$BaseGeneMean
    referenceMatrix <- referenceMatrix[order(baseMean),]
    baseMean <- order(baseMean)
    referenceMatrix <- t(apply(referenceMatrix, 1, sort))
    return(list(sim=sim, mat=referenceMatrix))

    # first generate matrix with values in [0,1]
    # then combine repeated genes
    # anchor to reference matrix


    # generate reference matrix
    # match row with closest mean expression
    # match column with closest relative expression in [0,1]

    # simulate dropout directly with splatter:::splatSimDropout
}

# eliminate duplicate genes by appying FUN across pathways that contain
#   the gene and store the result in the first pathway gene appeays in
removeGeneDuplicates <- function(pathways, FUN=max)
{
    pathways
}

testSplatter <- function()
{
    pwy1_genes <- paste('pwy1', letters, sep='_')
    pwy1 <- sapply(pwy1_genes, function(x) rgamma(1, shape=0.3, rate=0.6),
        USE.NAMES=T, simplify=F)

    pwy2_genes <- paste('pwy2', letters, sep='_')
    pwy2 <- sapply(pwy2_genes, function(x) rgamma(1, shape=0.3, rate=0.6))

    splatSimulateExtended(list(pwy1, pwy2), NULL)
}


