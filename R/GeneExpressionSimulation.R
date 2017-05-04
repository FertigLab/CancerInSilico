#' \code{inSilicoGeneExpression}
#' @param model a CellModel object
#' @param singleCell logical: simulate single cell data
#' @param pathways list of genes pathways
#' @param nCells number of cells to use for single cell data
#' @param fasta fasta file
#' @param ReferenceDataSet reference gene expression data
#' @param lambda parameter used for exponential
#' @param sampFreq how often to generate data
#' @param perError unknown
#' @param combineFUN function used to combine gene expression data
#' @param attrsep separation character
#' @export
#'
inSilicoGeneExpression <- function(model, singleCell = FALSE,
pathways = NULL, nCells = 96, fasta = NULL, ReferenceDataSet = NULL,
lambda = 1/3, sampFreq = 1, perError = 0.1, combineFUN = max,
attrsep = " ", microArray = FALSE, nGenes = NULL) {

    # simulate mean expression values for all pathways
    simMeanExprs <- simulateMeanExpression(model=model,
                        singleCell = singleCell, pathways=pathways, 
                        nCells = nCells, lambda = lambda,
                        ReferenceDataSet = ReferenceDataSet,
                        sampFreq = sampFreq, perError = perError,
                        combineFUN = combineFUN)
    
    # add dummy genes
    simMeanExprs <- padExpMatrix(simMeanExprs, nGenes)

    # add simulated error
    return (simulateError(simMeanExprs, pathways, fasta,
                ReferenceDataSet, perError, attrsep, microArray))


}

padExpMatrix <- function(meanExp, nGenes) {

    if (!is.null(nGenes) && nGenes > nrow(meanExp)) {

        newData <- matrix(nrow = nGenes, ncol = ncol(meanExp))
        newGenes <- c()
        for (i in 1:(nGenes - nrow(meanExp))) {
          
            newGenes[i] <- paste('dummy_', i, sep = "")
          
        }

        rownames(newData) <- sample(c(newGenes, rownames(meanExp)))

        for (g in rownames(newData)) {
      
            if (grepl('dummy', g)) {
         
               newData[g,] <- rexp(ncol(meanExp), 1) + 3.5

            } else {

                newData[g,] <- meanExp[g,]

            }

        }

    } else {

        newData <- meanExp

    }

    return (newData)

}

simulateError <- function(meanExp, pathways, fasta, 
ReferenceDataSet, perError, attrsep, microArray) {

    if (microArray) {

        sdMatrix <- pmax(perError * meanExp, perError)  
        outData <- meanExp + sdMatrix * matrix(rnorm(
                        length(meanExp)), nrow=nrow(meanExp))
    
    } else {

        # convert to counts from estimated log-transformed data
        meanExp <- round(2 ^ meanExp - 1)

        if (is.null(fasta)) {

            # simulate data with a negative binomial error model
            message(paste('Simulating time-course RNA-seq data with a ',
                            'negative binomial error model'))

            outData <- apply(meanExp, 2,
                            function(mu) {
                               NBsim(mu=mu,
                               ReferenceDataset=ReferenceDataSet)
                            })
        } else {

            # run polyester
            message('Simulating time-course RNA-seq data with Polyester')
    
            simulatePolyester(simMeanExprs=meanExp,fasta=fasta,
                            attrsep = attrsep, pathways = pathways)

        }

    }

    dimnames(outData) <- dimnames(meanExp)

    return(outData)

}

simulateMeanExpression <- function (model, singleCell = FALSE, nCells = 96,
pathways = NULL, ReferenceDataSet = NULL, lambda = 1/3, perError = 0.1,
sampFreq = 1, combineFUN = max) {

    # determine which pathways to simulate from the input values,
    # subsetting only to simulated pathways
    pathways <- getPathways(pathways)

    # set gene expression values to use for each pathway
    pathways <- setPathwayExpressionRange(
                    ReferenceDataSet = ReferenceDataSet, lambda = lambda,
                    pathways = pathways)

    # run simulation for each pathway
    pathwayOutput <- list()
    for (path in names(pathways)) {

        pathwayOutput[[path]] <- simulatePathway(model = model,
                                    pathway = pathways[[path]],
                                    type = path, sampFreq = sampFreq,
                                    sampSize = nCells,
                                    singleCell = singleCell)

    }

    # return combined expression
    return (combineGeneExpression(pathwayOutput, combineFUN))

}

NBsim <- function(mu, ReferenceDataset = NULL, invChisq = TRUE) {

    if (is.null(ReferenceDataset)) {

        ngenes<-length(mu)
        mu.Reference <- pmax(mu,1)

    } else {
    
        # check that the reference dataset is a valid log transformed dataset
        # with values for every gene in mu
        checkReferenceDataset(ReferenceDataset, list(row.names(mu)))

        #number of genes  
        ngenes<-dim(ReferenceDataset)[1]

        # assume that the reference dataset is log transformed so 
        # convert to counts
        mu.Reference<-pmax(round(rowMeans(2^ReferenceDataset - 1)),1)

    }

    # Biological variation
    BCV0 <- 0.2+1/sqrt(mu.Reference)
  
    if (invChisq) {
    
        df.BCV <- 40
        BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV))

    } else {

        BCV <- BCV0*exp(rnorm(ngenes,mean=0,sd=0.25)/2 )

    }


    #if(NCOL(BCV)==1) BCV <- matrix(BCV,ngenes)
    shape <- 1/BCV^2
    scale <- mu/shape
    mu <- matrix(rgamma(ngenes,shape=shape,scale=scale),ngenes)

    # Technical variation
    counts <- matrix(rpois(ngenes,lambda=mu),ngenes)
    return(counts)

}

# combine gene expression matrices according to some function
combineGeneExpression <- function(geneExpression, combineFUN=max) {

    # verify same number of columns
    col_num <- sapply(geneExpression, ncol)
    if (min(col_num) != max(col_num)) {

        stop("unequal number of columns in gene expression matrices")

    }

    # find number of unique genes
    gene_list <- lapply(geneExpression, row.names)
    total_genes <- unique(unlist(gene_list))
    
    # initalize matrix
    output <- matrix(nrow = length(total_genes), ncol = col_num[1]) 
    row.names(output) <- total_genes
    colnames(output) <- colnames(geneExpression[[1]])
 
    # combine matrices, iterate through every gene
    for (g in total_genes) {

        # find which pathways have this gene
        pwys <- sapply(lapply(gene_list, is.element, g), sum) > 0
    
        # get expression values from each pathway
        exp <- sapply(geneExpression[pwys], function(x) {x[g,]})

        # combine expression values
        output[g, ] <- apply(exp, 1, combineFUN)

    }   

    return (output)
 
}

simulatePolyester <- function(simMeanExprs=simMeanExprs,fasta=fasta, attrsep = attrsep,
                              idfield='gene_symbol', outdir=".", pathways=NULL) {
  
  # process fasta to create mean values for transcripts
  transcripts = readDNAStringSet(fasta) # read in the fasta
  tmpsplit=strsplit(names(transcripts),split=attrsep) # split the fasta file into elements
  
  tmpl=unlist(lapply(tmpsplit, function(x) x[grep(idfield,x)])) #grab the gene name element
  
  geti=function(x,i) return(x[i]) #supportive function
  
  #gnames may be created with something too specific to this particular one?
  gnames=sapply(strsplit(as.character(tmpl),split=":",fixed=T),geti,i=2) #grab the second element after the colon, the gene name
  names(gnames) <- sapply(tmpsplit,geti,1) # name as transcript ids
  
  # check that using valid columns of the fasta file
  if (is.null(pathways)) {
    pathways <- inSilicoPathways
  }
  
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
