#' \code{getPathwaysSim} Determine the pathways from which to simulate gene expression data and check its accuracy
#'
#' @param pathway Optional; A list of genes associated with pathways, defaulting to gene symbols in the package data object \code{\link{inSilicoPathways}}. If specified by the user, names of pathways should be GtoM for genes associated with the G to M cell cycle transition, GtoS for genes associated with G to S, Prox for genes associated with contact, and Growth for genes associated with growth factor receptor signaling.
#' @return A list of genes associated with pathways that will be used for the simulation.
#' 
getPathwaysSim <- function(pathways) {

  pathSim <- c('GtoM', 'GtoS', 'Prox', 'Growth')
  
  if (is.null(pathways)) {
    # load default pathways 
    data('inSilicoPathways', package = 'CancerInSilico')
  } else {
    # find pathways that should be simulated with CancerInSilico
    pathways <- pathways[intersect(names(pathways), pathSim)]
    
    # check that at least 1 pathway matches the simulated pathway in the model
    if (length(pathways) == 0) {
      stop(paste('Names of user defined pathways do not match',
                 paste(pathSim, collapse=', '), 'simulated in CancerInSilico.'))
    }
    
    # warn users if any pathways will be excluded from the simulation
    if (length(pathways) < length(pathSim)) {
      warning(paste(paste(setdiff(pathSim,names(pathways)), collapse=", "),
                    "are not specified in the user defined pathways. CancerInSilico will simulate array data from only",
                    paste(names(pathways), collapse=", ")))
    }
    
  }
  
  return(pathways)
  
}
