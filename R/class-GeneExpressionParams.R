#' @title GeneExpressionParams
#' @export
#'
#' @description Parameters for simulating gene expression
#' @slot sampleFreq how often to generate data
#' @slot RNAseq generate RNA-seq data
#' @slot singleCell generate single cell data
#' @slot nCells number of cells to sample at each time point
#' @slot nDummyGenes number of dummy genes
#' @slot dummyDist function to determine expression of dummy genes
#' @slot combineFUN function used to combine gene expression data
#' @slot randSeed random seed
#' @slot perError error for normal error model
#' @slot bcvCommon error for voom error model
#' @slot bcvDF degrees of freedom for voom error model
#' @slot dropoutPresent whether to simulate dropout in single cell data
#' @slot dropoutMid parameter for dropout distribution
#' @slot dropoutShape parameter for dropout distribution
setClass("GeneExpressionParams", slots = c(
    sampleFreq = "numeric",
    RNAseq = "logical",
    singleCell = "logical",
    nCells = "numeric",
    nDummyGenes = "numeric",
    dummyDist = "function",
    combineFUN = "function",
    randSeed = "numeric",
    perError = "numeric",
    bcvCommon = "numeric",
    bcvDF = "numeric",
    dropoutPresent = "logical",
    dropoutMid = "numeric",
    dropoutShape = "numeric"
))

setMethod("initialize", "GeneExpressionParams", 
    function(.Object, ...)
    {
        .Object@sampleFreq <- 1
        .Object@RNAseq <- FALSE
        .Object@singleCell <- FALSE
        .Object@nCells <- 100
        .Object@nDummyGenes <- 0
        .Object@dummyDist <- function(N) rep(0, N)
        .Object@combineFUN <- max
        .Object@randSeed <- 0
        .Object@perError <- 0.1
        .Object@bcvCommon <- 0.2
        .Object@bcvDF <- 40
        .Object@dropoutPresent <- FALSE
        .Object@dropoutMid <- 0
        .Object@dropoutShape <- -1

        .Object <- callNextMethod(.Object, ...)
        .Object
    }
)

setValidity('GeneExpressionParams',
    function(object)
    {
        if (object@singleCell & !is.null(object@fasta))
            'polyester is not for single cell'
        else if (!object@RNAseq & object@singleCell)
            'can\'t generate microarray data for single cell'
        else if (!object@RNAseq & !is.null(object@fasta))
            'polyester is only for RNA-seq data'
        else
            TRUE
    }
)

##################### Generics ###################

##################### Methods ####################