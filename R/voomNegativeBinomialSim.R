### Function to simulate negative binomial distribution from count data
NBsim<- function(mu, # mean expression values for each gene
                 ReferenceDataset=NULL, # data matrix of gene expression by sample
                 invChisq=TRUE, # Use inverse chi-square or log-normal dispersion
                 ...){
  
  
  
  if(is.null(ReferenceDataset)){
    ngenes<-length(mu)
    mu.Reference <- pmax(mu,1)
  }
  else {
    
    # check that the reference dataset is a valid log transformed dataset
    # with values for every gene in mu
    checkReferenceDataset(ReferenceDataset, pathways = list(row.names(mu)))
    
    ngenes<-dim(ReferenceDataset)[1]
    
    # assume that the reference dataset is log transformed so 
    # convert to counts
    mu.Reference<-pmax(round(rowMeans(2^ReferenceDataset - 1)),1)
  }
  
  # Biological variation
  BCV0 <- 0.2+1/sqrt(mu.Reference)
  if(invChisq){
    df.BCV <- 40
    BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV))
  } else {
    BCV <- BCV0*exp(rnorm(ngenes,mean=0,sd=0.25)/2 )
  }
  if(NCOL(BCV)==1) BCV <- matrix(BCV,ngenes)
  shape <- 1/BCV^2
  
  
  scale <- mu/shape
  mu <- matrix(rgamma(ngenes,shape=shape,scale=scale),ngenes)
  
  # Technical variation
  counts <- matrix(rpois(ngenes,lambda=mu),ngenes)
  
  return(counts)
}