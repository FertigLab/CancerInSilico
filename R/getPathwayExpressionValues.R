###CancerInSilico 8.10.2016
###Emily Flam / Elana J. Fertig
###Find sample with max median value of raw RNA expression and expression Z-score for genes 
### for all pathways to use in gene expression simulation
getPathwayExpressionValues <- function(ReferenceDataset=NULL, exponentialModel=T, pathways) {
  
  if (is.null(ReferenceDataset)) {
    if (exponentialModel) {
        ### Raymond fill this in or kill it as needed
        ### This would be generating the values for each gene from
        ### The exponenential model 
    } else {
       ### Use specified values from an internal package dataset
    }
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

