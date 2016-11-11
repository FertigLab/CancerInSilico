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

combineGeneExpression <- function(geneExpression, combineFUN=max) {

    # verify same number of columns
    col_num <- sapply(geneExpression, ncol)
    if (min(col_num) != max(col_num)) {

        stop("unequal number of columns in gene expression matrices")

    }

    # find number of unique genes
    gene_list <- sapply(geneExpression, row.names)
    total_genes <- unique(unlist(gene_list))
    
    # initalize matrix
    output <- matrix(nrow = length(total_genes), ncol = col_num[1]) 
    row.names(output) <- genes
    colnames(output) <- colnames(geneExpression[[1]])
 
    # combine matrices, iterate through every gene
    for (g in total_genes) {

        # find which pathways have this gene
        pwys <- sapply(sapply(gene_list, is.element, g), sum) > 0
    
        # get expression values from each pathway
        exp <- sapply(geneExpression[pwys], function(x) {x[g,]})

        # combine expression values
        output[g, ] <- apply(exp, 1, combineFUN)

    }
 
}
