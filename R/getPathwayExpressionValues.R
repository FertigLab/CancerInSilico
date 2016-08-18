#' \code{getPathwayExpressionValues} Determine the gene expression values to use for each gene in each pathway.
#' @param pathway A list of genes associated with pathways.
#' @param ReferenceDataSet Optional; Reference gene expression dataset to use to calculate the gene expression values for genes in the pathway to use as mean values in the simulation. Defaults values randomly selected from an exponential distribution with parameter \code{lambda}. If specified by the user, \code{row.names} of the ReferenceDataset must match gene names in the \code{pathway} argument. In this case, gene expression values will be set as the values for genes in the pathway of the sample with max median value of raw RNA expression and expression Z-score for pathway genes. 
#' @param lambda Optional; Parameter of the exponential distribution used to determine the maximum expression value of each simulated gene in the pathway. Defaults to 1/3. Not used if values are determined from a dataset in \code{ReferenceDataSet}.
#' @return list of numeric objects for the maximum expression value of each gene in each pathway. 
#' 

getPathwayExpressionValues <- function(ReferenceDataset=NULL, lambda=1/3, pathways) {
  
  if (is.null(ReferenceDataset)) {
    
    pathValues <- lapply(pathways,function(x){
      val <- rexp(length(x),rate=lambda)
      names(val) <- x
      return(val)
    })
    
    return(pathValues)

  } else {
    D <- sweep(ReferenceDataset,1,apply(ReferenceDataset,1,max),FUN="/")
    
    # limit pathways to genes contained in the reference dataset
    pathways <- lapply(pathways, intersect, row.names(ReferenceDataset))
    
    # check that the reference dataset is a valid, log transformed dataset
    checkReferenceDataset(ReferenceDataset,pathways)
    
    if (any(sapply(pathways,length)==0)) {
      stop(paste('The following pathways do not have any genes in the reference dataset: ',
                 paste(names(which(sapply(pathways,length)==0)), collapse=",")))
    }
    
    # for each pathway return the gene expression values in the reference dataset for the sample
    # that has the maximum median gene expression for all pathway genes
    return(lapply(pathways,function(x){ReferenceDataset[x,names(which.max(apply(D[x,],2,median)))]}))

  }
}

