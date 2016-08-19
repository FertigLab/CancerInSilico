#' \code{inSilicoRNASeq} Simulate gene expression data that would be generated from a RNA-seq
#'
#' @param model A \code{\link{CellModel}}
#' @param pathways Optional; A list of genes associated with pathways, defaulting to gene symbols in the package data object \code{\link{inSilicoPathways}}. If specified by the user, names of pathways should be GtoM for genes associated with the G to M cell cycle transition, GtoS for genes associated with G to S, Prox for genes associated with contact, and Growth for genes associated with growth factor receptor signaling.
#' @param ReferenceDataSet Optional; Reference gene expression dataset to use to calculate the gene expression values for genes in the pathway to use as mean values in the simulation. Defaults values randomly selected from an exponential distribution with parameter \code{lambda}. If specified by the user, \code{row.names} of the ReferenceDataset must match gene names in the \code{pathway} argument.
#' @param fasta Optional; NULL value will create a matrix of counts with a negative binomial error model. Otherwise, path to FASTA file containing transcripts from which to simulate reads. See details in \code{\link{simulate_experiment_countmat}}.
#' @param lambda Optional; Parameter of the exponential distribution used to determine the maximum expression value of each simulated gene in the pathway. Defaults to 1/3. Not used if values are determined from a dataset in \code{ReferenceDataSet}.
#' @param samplingFrequency Optional; Distance between time points in hours at which gene expression values are simulated. Defaults to every hour.
#' @param combineFUN Optional; One of max, sum, mean, or median to specify the function used to combine expression values from genes shared in multiple pathways. Defaults to maximum value in any condition (max).
#' @param attrsep Optional; string separating gene attributes in the fasta file for \code{\link{polyester}}. Defaults to " ".
#' @param ... other input arguments to combineFUN and \code{\link{polyester}}
#' @return Matrix simulating the gene expression values for pathway genes from a \code{\link{CellModel}} simulation at the specified frequency in \code{samplingFrequency} if fasta is the default value of NULL. Otherwise, fastq files will be generated to simulate RNA-seq data as described in \code{\link{simulate_experiment_countmat}}.
#' @export
#' 

inSilicoRNASeq <- function(CellModels, pathways=NULL, ReferenceDataSet=NULL, fasta=NULL,
                          samplingFrequency=1, lambda=1/3, combineFUN=max, attrsep = " ", ...){
  
  # call standard function to simulate mean expression values for all pathways
  simMeanExprs <- simulateMeanExpression(CellModels=CellModels, pathways=pathways, 
                                         lambda=lambda, ReferenceDataSet=ReferenceDataSet, 
                                         samplingFrequency=samplingFrequency, combineFUN=combineFUN,
                                         ...)
  
  # convert to counts from estimated log-transformed data
  simMeanExprs <- round(2^simMeanExprs-1)
  
  
  if (is.null(fasta)) {
    # simulate data with a negative binomial error model
    message('Simulating time-course RNA-seq data with a negative binomial error model')
    
    outData <- apply(simMeanExprs,2,function(mu){NBsim(mu=mu,
                                                       ReferenceDataset=ReferenceDataSet,
                                                       ...)})
    
    dimnames(outData) <- dimnames(simMeanExprs)
    
    return(outData)
    
  } else { 
    
    # run polyester
    message('Simulating time-course RNA-seq data with Polyester')
    
    simulatePolyester(simMeanExprs=simMeanExprs,fasta=fasta, attrsep = attrsep, ...)
  }
  
}