#' @title GeneExpressionParams
#' @description Parameters for simulating gene expression
#'
#' @slot sampFreq how often to generate data
#' @slot nGenes total number of genes used (matrix padded with dummys)
#' @slot combineFUN function used to combine gene expression data
#' @slot singleCell logical: simulate single cell data
#' @slot nCells number of cells to use for single cell data
#' @slot perError percent error?
#' @slot microArray true if micro array data
#' @slot randSeed random seed for simulation
#' @export
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