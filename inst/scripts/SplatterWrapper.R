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

library(splatter)
data('sc_example_counts')
params <- splatEstimate(sc_example_counts)
sim <- splatSimulate(params, dropout.present = FALSE)
sim@featureData@data


params <- newSplatParams()
paramValues <- list(
  #nGroups = 1,
  groupCells = 10,
  mean.shape = 0.6,
  mean.rate = 0.3,
  lib.loc = 11,
  lib.scale = 0.2,
  out.prob = 0,
  out.facLoc = 4,
  out.facScale = 0.5,
  de.prob = 0,
  de.downProb = 0,
  de.facLoc = 0.1,
  de.facScale = 0.4,
  bcv.common = 0,
  bcv.df = 1,
  dropout.present = FALSE,
  dropout.mid = 0,
  dropout.shape = -1,
  path.from = 0,
  path.length = 1,
  path.skew = 0.5,
  path.nonlinearProb = 0.1,
  path.sigmaFac = 0.8,
  nGenes = 10,
  #nCells = 10,
  seed = 827593
)
params <- setParams(params, update = paramValues)
sim <- splatSimulate(params, method = "single")

mat <- counts(sim)
mat <- sim@assayData$BaseCellMeans
geneMin <- apply(mat, 1, min)
geneMax <- apply(mat, 1, max)
mat2 <- matrix(nrow=10, ncol=10)

for (r in 1:10)
{
  mat2[r,] <- (mat[r,] - geneMin[r]) / (geneMax[r] - geneMin[r])
}
mat
geneMin
geneMax
mat2
apply(mat2, 2, mean)
apply(mat2, 2, sd)

getPathwayActivity <- function(sim, id=1)
{
  gene1 <- sim@assayData$BaseCellMeans[id,]
  #gene1 <- sim@assayData$TrueCounts[id,]
  return((gene1 - min(gene1)) / (max(gene1) - min(gene1)))
}
getPathwayActivity(sim)


library(splatter)
sim <- splatSimulate(nGenes=5, groupCells=5, dropout.present=TRUE, out.prob=0.5)
mat <- counts(sim)
mat
sim@assayData$BaseCellMeans
sim@assayData$TrueCounts
sim@featureData@data$BaseGeneMean
sim@phenoData@data$ExpLibSize

round(getPathwayActivity(sim,5), 4)

params <- newSplatParams()
params <- setParam(params, "nGenes", 10)
params <- setParam(params, "groupCells", 10)
means <- rgamma(params@nGenes, shape = params@mean.shape,
   rate = params@mean.rate)
means
sim@featureData@data$BaseGeneMean
counts <- matrix(rnbinom(params@nGenes * params@nCells, mu = means,
    size = 1 / params@bcv.common), nrow = params@nGenes)
counts

sim <- simpleSimulate(nCells=10, nGenes=10)
counts(sim)

sim <- splatSimulate(nGenes=1000, groupCells=100, mean.shape = 1.0,
                     mean.rate = 0.3)
mat <- counts(sim)
mat <- sim@featureData@data$BaseGeneMean
gene_means <- mat #rowMeans(mat)
plot(density(rgamma(1e6, shape = 1.0, rate = 0.3)), col='red')
lines(density(gene_means))

sim@featureData@data$BaseGeneMean[1:10]
sim@featureData@data$GeneMean[1:10]


# notes

# pathway sim function takes single time
# loop over times in GeneExpression code

# method
# for each pathway, at each sample (time)
#   - fit gamma distribution to expression range parameter
#   - simulate base expression with minimal possible variance
#   - @assayData$BaseCellMeans contains "pathway activity"
#   - line up genes based on mean expression
#   - line up cells based on pathway activity
#   - simulate all other effects directly
# alternatively, do full simulation and match genes & cells by
#   looking at intermediate data

# params to estimate from data
#   - dropout 



