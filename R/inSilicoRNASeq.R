#' \code{inSilicoRNASeq} Simulate gene expression data that would be generated from a RNA-seq
#'
#' @param model A \code{\link{CellModel}}
#' @param pathway Optional; A list of genes associated with pathways, defaulting to gene symbols in the package data object \code{\link{inSilicoPathways}}. If specified by the user, names of pathways should be GtoM for genes associated with the G to M cell cycle transition, GtoS for genes associated with G to S, Prox for genes associated with contact, and Growth for genes associated with growth factor receptor signaling.
#' @param ReferenceDataSet Optional; Reference gene expression dataset to use to calculate the gene expression values for genes in the pathway to use as mean values in the simulation. Defaults values randomly selected from an exponential distribution with parameter \code{lambda}. If specified by the user, \code{row.names} of the ReferenceDataset must match gene names in the \code{pathway} argument.
#' @param fasta Optional; NULL value will create a matrix of counts with a negative binomial error model. Otherwise, path to FASTA file containing transcripts from which to simulate reads. See details in \code{\link{simulate_experiment_countmat}}.
#' @param samplingFrequency Optional; Distance between time points in hours at which gene expression values are simulated. Defaults to every hour.
#' @param combineFUN Optional; One of max, sum, mean, or median to specify the function used to combine expression values from genes shared in multiple pathways. Defaults to maximum value in any condition (max).
#' @param attrsep Optional; string separating gene attributes in the fasta file for \code{\link{polyester}}. Defaults to " ".
#' @param ... other input arguments to combineFUN and \code{\link{polyester}}
#' @return Matrix simulating the gene expression values for pathway genes from a \code{\link{CellModel}} simulation at the specified frequency in \code{samplingFrequency} if fasta is the default value of NULL. Otherwise, fastq files will be generated to simulate RNA-seq data as described in \code{\link{simulate_experiment_countmat}}.
#' @export
#' 

inSilicoRNASeq <- function(CellModels, pathway=NULL, ReferenceDataSet=NULL, fasta=NULL,
                          samplingFrequency=1, combineFUN="max", attrsep = " ", ...){
  
  # call standard function to simulate mean expression values for all pathways
  simMeanExprs <- simulateMeanExpression(CellModels=CellModels, pathways=pathways, 
                                         lambda=lambda, ReferenceDataSet=ReferenceDataSet, 
                                         samplingFrequency=samplingFrequency, combineFUN=combineFUN,
                                         ...)
  
  
  if (is.null(fasta)) {
    # simulate data with a negative binomial error model
    message('Simulating time-course RNA-seq data with a negative binomial error model')
    
    outData <- apply(simMeanExprs,2,function(mu){NBsim(mu=mu,
                                                       ReferenceDataset=ReferenceDataSet,
                                                       ...)})
    
    return(outData)
    
  } else { 
    
    # run polyester
    message('Simulating time-course RNA-seq data with Polyester')
    
    # process fasta to create mean values for transcripts
    transcripts = readDNAStringSet(myfasta) # read in the fasta
    tmpsplit=strsplit(names(transcripts),split=attrsep) # split the fasta file into elements
    
    tmpl=unlist(lapply(tmpsplit, function(x) x[grep(fastaGeneColumn,x)])) #grab the gene name element
    
    geti=function(x,i) return(x[i]) #supportive function
    
    #gnames may be created with something too specific to this particular one?
    gnames=sapply(strsplit(as.character(tmpl),split=":",fixed=T),geti,i=2) #grab the second element after the colon, the gene name
    
    # check that using valid columns of the fasta file
    if (!all(unique(unlist(pathways)) %in% gnames)) {
      stop(paste(fastaGeneColumn, 'does not correspond to the annotation of genes in pathways'))
    }
    
    # set the newmu to the length of the fasta, set all values to 0, 
    # the mus are then put in, divided by the number of transcripts for each gene
    newmu <- matrix(0, nrow=length(gnames), ncol=ncol(simMeanExprs),
                    dimnames = list(gnames, colnames(simMeanExprs)))
    
    for (j in 1:ncol(simMeanExprs)) {
 
      
      for(i in intersect(names(newmu),names(mu))){ #for each gene that is in common between the input and the fasta,
        newmu[which(gnames==i),j]=simMeanExprs[i,j]/length(which(row.names(newmu)==i)) #take the provided gene data to each element with that gene of the fasta, divide by how many elements are that gene in the fasta
        
      }
      
    }
    

    # run Polyester to simulate gene expression data
    simulate_experiment_countmat(fasta=myfasta,
                                 num_reps=1, readmat= newmu,
                                 outdir=outdirs,idfield=fastaGeneColumn, ...)
  }
  
}