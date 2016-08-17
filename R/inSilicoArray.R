#' \code{inSilicoArray} Simulate gene expression data that would be generated from a microarray
#'
#' @param model A \code{\link{CellModel}}
#' @param pathway Optional; A list of genes associated with pathways, defaulting to gene symbols in the package data object \code{\link{inSilicoPathways}}. If specified by the user, names of pathways should be GtoM for genes associated with the G to M cell cycle transition, GtoS for genes associated with G to S, Prox for genes associated with contact, and Growth for genes associated with growth factor receptor signaling.
#' @param perError Optional; Percentage of the mean value to be used as the standard deviation in the error model to simulate normally distributed array data for a CellModels simulation. Defaults to 10% of the average expression value.
#' @param ReferenceDataSet Optional; Reference gene expression dataset to use to calculate the gene expression values for genes in the pathway to use as mean values in the simulation. Defaults values randomly selected from an exponential distribution with parameter \code{lambda}. If specified by the user, \code{row.names} of the ReferenceDataset must match gene names in the \code{pathway} argument.
#' @param lambda Optional; Parameter of the exponential distribution used to determine the maximum expression value of each simulated gene in the pathway. Defaults to 1/3. Not used if values are determined from a dataset in \code{ReferenceDataSet}.
#' @param samplingFrequency Optional; Distance between time points in hours at which gene expression values are simulated. Defaults to every hour.
#' @param combineFUN Optional; One of max, sum, mean, or median to specify the function used to combine expression values from genes shared in multiple pathways. Defaults to maximum value in any condition (max).
#' @param ... other input arguments to combineFUN
#' @return Matrix simulating the gene expression values for pathway genes from a \code{\link{CellModel}} simulation at the specified frequency in \code{samplingFrequency}.
#' @export
#' 

inSilicoArray <- function(CellModels, pathways=NULL, lambda=1/3, perError=0.1, ReferenceDataSet=NULL, 
                          samplingFrequency=1, combineFUN="max", ...){
  
  # call standard function to simulate mean expression values for all pathways
  simMeanExprs <- simulateMeanExpression(CellModels=CellModels, pathways=pathways, 
                                         lambda=lambda, ReferenceDataSet=ReferenceDataSet, 
                                         samplingFrequency=samplingFrequency, combineFUN=combineFUN,
                                         ...)
  
  # simulate data with a normal error model
  sdMatrix <- pmax(perError*simMeanExprs, perError)
  outData <- simMeanExprs + sdMatrix*rnorm(length(simMeanExprs), nrow=nrow(simMeanExprs))
  
  # return simulated data matrix
  return(outData)
}