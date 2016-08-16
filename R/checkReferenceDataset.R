checkReferenceDataset <- function(ReferenceDataset, pathways=NULL) {

  ### Reference dataset is a gene expression dataset used to find values
  ### should be assumed to be log transformed -- may want to have a check that does this in case!!!
  if (max(ReferenceDataset)>50) {
    warning(paste('ReferenceDataset has a maximum value of',
                  max(ReferenceDataset),
                  'check that the data is log transformed'))
  }
  
  if (min(ReferenceDataset) < 0) {
    warning('ReferenceDataset should be strictly non-negative.')
  }
  
  if (!is.null(pathways)) {
    if (any(is.na(ReferenceDataset[intersect(unlist(pathways),row.names(ReferenceDataset)),]))) {
      stop('ReferenceDataset contains NA values for pathway genes.')
    }
  }
}