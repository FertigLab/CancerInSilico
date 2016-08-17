#' \code{simulateMeanExpression} Standardized function to simulate mean gene expression values for all datatypes
#'
#' @param model A \code{\link{CellModel}}
#' @param pathway Optional; A list of genes associated with pathways, defaulting to gene symbols in the package data object \code{\link{inSilicoPathways}}. If specified by the user, names of pathways should be GtoM for genes associated with the G to M cell cycle transition, GtoS for genes associated with G to S, Prox for genes associated with contact, and Growth for genes associated with growth factor receptor signaling.
#' @param ReferenceDataSet Optional; Reference gene expression dataset to use to calculate the gene expression values for genes in the pathway to use as mean values in the simulation. Defaults values randomly selected from an exponential distribution with parameter \code{lambda}. If specified by the user, \code{row.names} of the ReferenceDataset must match gene names in the \code{pathway} argument.
#' @param lambda Optional; Parameter of the exponential distribution used to determine the maximum expression value of each simulated gene in the pathway. Defaults to 1/3. Not used if values are determined from a dataset in \code{ReferenceDataSet}.
#' @param samplingFrequency Optional; Distance between time points in hours at which gene expression values are simulated. Defaults to every hour.
#' @param combineFUN Optional; Function used to combine expression values from genes shared in multiple pathways. Defaults to maximum value in any condition.
#' @param ... other input arguments to combineFUN
#' @return Matrix simulating the mean gene expression values for pathway genes from a \code{\link{CellModel}} simulation at the specified frequency in \code{samplingFrequency}.
#' 
simulateMeanExpression <- function(CellModels, pathways=NULL, lambda=1/3, ReferenceDataSet=NULL,
                                   samplingFrequency=1, combineFUN=max, ...){
  

  # determine which pathways to simulate from the input values, subsetting only to simulated pathways
  pathways <- getPathwaysSim(pathways)

  # get gene expression values to use for each pathway
  pathwayExprsValues <- getPathwayExpressionValues(ReferenceDataset=ReferenceDataset,
                                                   lambda=lambda, pathways=pathways)

  # run simulation for each pathway
  pathSimOutput <- list()
  for (path in names(pathwayExprsValues)) {
    pathSimOutput[[path]] <- eval(parse(text=sprintf("simulate%sPathGroup(model, pathwayExprsValues[[path]], samplingFrequency)",
                                                     path)))
  }

  # merge values from all simulated pathways
  pathwaysCombine <- combinePathways(pathSimOutput, combineFUN=combineFUN, ...)
  
  return(pathwaysCombine)
}

#' \code{combinePathways} Combine gene expression values for all simulated pathways
#'
#' @param pathSimOutput List of matrices containing gene by time expression values for each pathway generated from \code{simulateMeanExpression}.
#' @param combineFUN Optional; Function used to combine expression values from genes shared in multiple pathways. Defaults to maximum value in any condition.
#' @param ... other input arguments to combineFUN
#' @return Matrix simulating the mean gene expression values for pathway genes from a \code{\link{CellModel}} simulation at the specified frequency in \code{samplingFrequency}.
#' 
combinePathways <- function(pathSimOutput, combineFUN=max, ...){
  
  # initialize matrix
  outMatrix <- matrix(0., nrow=length(unique(unlist(sapply(pathSimOutput,row.names)))),
                      ncol=ncol(pathSimOutput[[1]]),
                      dimnames=list(unique(unlist(sapply(pathSimOutput,row.names))),
                                    colnames(pathSimOutput[[1]])))
  
  # apply the combination function to merge the data simulated from different pathways
  pathSimCombine <- pathSimOutput
  for (path in names(pathSimOutput)) {
    pathSimCombine[[path]] <- outMatrix
    pathSimCombine[[path]][row.names(pathSimOutput[[path]]),] <- pathSimOutput[[path]]
  }
  
  # combine all the pathways according to the specified function
  for (t in 1:ncol(outMatrix)) {
    for (g in row.names(outMatrix)){
      outMatrix[g,t] <- combineFUN(sapply(pathSimCombine, function(x){x[g,t]}), ...)
    }
  }
  
  return(outMatrix)
  
}
