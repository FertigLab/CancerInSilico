simulatePolyester <- function(simMeanExprs=simMeanExprs,fasta=fasta, attrsep = attrsep,
                              idfield='gene_symbol', outdir=".", ...) {
  
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
    
    for (i in intersect(unique(gnames),
                        row.names(simMeanExprs)[simMeanExprs[,j]>0])) { #for each gene that is in common between the input and the fasta,
      # distribute reads randomly across transcripts for the same gene
      # distributed according to the relative size of each transcript
      newmu[which(gnames==i),j] <- round(simMeanExprs[i,j]*width(transcripts)[gnames==i] / 
        (sum(width(transcripts)[gnames==i])))
    }
    
  }
  
  # run Polyester to simulate gene expression data
  for (i in 1:ncol(newmu)) {
    simulate_experiment(fasta=fasta,reads_per_transcript=1,
                        num_reps=1, fold_changes = newmu[,i], 
                        attrsep = attrsep, outdir=outdir) #, ...)
    
    fnnew <- paste0('s_',sprintf('%02d',i))
    if (!is.null(colnames(newmu))) {
      fnnew <- sprintf('s_%s_%02d',colnames(newmu)[i],i)
    }
  
    fns <- list.files(outdir,pattern=paste0('sample_01'))
    file.rename(file.path(outdir, fns),
                file.path(outdir, sub('sample_01',fnnew,fns)))
  }
  
  
}