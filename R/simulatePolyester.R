simulatePolyester <- function(simMeanExprs=simMeanExprs,fasta=fasta, attrsep = attrsep,
                              idfield='gene_symbol', ...) {
  
  # process fasta to create mean values for transcripts
  transcripts = readDNAStringSet(fasta) # read in the fasta
  tmpsplit=strsplit(names(transcripts),split=attrsep) # split the fasta file into elements
  
  tmpl=unlist(lapply(tmpsplit, function(x) x[grep(idfield,x)])) #grab the gene name element
  
  geti=function(x,i) return(x[i]) #supportive function
  
  #gnames may be created with something too specific to this particular one?
  gnames=sapply(strsplit(as.character(tmpl),split=":",fixed=T),geti,i=2) #grab the second element after the colon, the gene name
  names(gnames) <- sapply(tmpsplit,geti,1) # name as transcript ids
  
  # check that using valid columns of the fasta file
  if (!any(unique(unlist(pathways)) %in% gnames)) {
    stop(paste(idfield, 'does not correspond to the annotation of genes in pathways'))
  }
  
  # set the newmu to the length of the fasta, set all values to 0, 
  # the mus are then put in, divided by the number of transcripts for each gene
  newmu <- matrix(0, nrow=length(gnames), ncol=ncol(simMeanExprs),
                  dimnames = list(names(gnames), colnames(simMeanExprs)))
  
  for (j in 1:ncol(simMeanExprs)) {
    
    for (i in intersect(unique(gnames),row.names(simMeanExprs))){ #for each gene that is in common between the input and the fasta,
      newmu[which(gnames==i),j]=simMeanExprs[i,j]/length(which(gnames==i)) #take the provided gene data to each element with that gene of the fasta, divide by how many elements are that gene in the fasta
      
    }
    
  }
  
  # run Polyester to simulate gene expression data
  simulate_experiment_countmat(fasta=fasta,
                               num_reps=1, readmat= newmu, 
                               attrsep = attrsep, ...)
  
}