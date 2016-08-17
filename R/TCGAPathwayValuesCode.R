###CancerInSilico 8.10.2016
###Emily Flam
###Find sample with max median value of raw RNA expression and expression Z-score for genes in given pathway
###Generate heatmap

library(gplots)
library('GSA')

#Load TCGA HNC RNA expression data
load("hnRNAseqNormV2Jul2013.Rda")
RNAexp <- hnRNAseqNormV2

#load HPV +/- pheno data
pheno <- read.table("hpv_status_tumor_site_1_14_14.txt", header = T, sep = "\t")

#Create matrix of RNA expression data for only Tumor samples
#remove "-01" from sample name to match sample name in pheno dataframe 
TumorExp <- RNAexp[,which(substr(colnames(RNAexp),14,15)=="01")]
colnames(TumorExp) <- c(substr(colnames(TumorExp),1,12))

#Create vector of HPV status phenotype of tumor samples
# +/- as value, sample name as name
#change sample name to have - instead of . separation in order to match the rest of the data
HPVpheno <- as.character(substr(pheno$New.HPV.Status,4,4))
names(HPVpheno) <- as.character(gsub("\\.","-",pheno$Barcode))

#Create RNA expression matrix for only HPV- tumor samples
TumorNegExp <- ReferenceDataset <- TumorExp[,which(colnames(TumorExp) %in% names(which(HPVpheno == "-")))]

#Create vector of desired gene names
#Random 7 genes chosen here for example 
GeneSetMatrix <-read.table('E2FpathwayGenes_hallmark.txt', header=T, sep="\t",
                           stringsAsFactors = F)
GeneSet <- as.character(GeneSetMatrix[2:nrow(GeneSetMatrix),1])

pathways <- list('E2F'=GeneSet)

## call the function to get pathway expression values
pathwayValues <- getPathwayExpressionValues(ReferenceDataset=TumorNegExp,
                                            pathways=pathways)

## check the negative binomial function for one value
source('voomNegativeBinomialSim.R')
sapply(1:1000,NBsim,mu=2^pathwayValues$E2F-1,ReferenceDataset=NULL,invChisq=T) -> x
par(mfrow=c(1,2))
plot(rowMeans(x),2^pathwayValues$E2F-1)
lines(rowMeans(x),rowMeans(x),type='l',col='red')
plot(apply(x,1,sd),2^pathwayValues$E2F-1,xlim=c(0,max(2^pathwayValues$E2F-1)))
lines(rowMeans(x),rowMeans(x),type='l',col='red')

## checking negative binomial function with reference
sapply(1:1000,NBsim,mu=2^pathwayValues$E2F-1,
       ReferenceDataset=TumorNegExp[names(pathwayValues$E2F),],invChisq=T) -> x
plot(rowMeans(x),2^pathwayValues$E2F-1)
lines(rowMeans(x),rowMeans(x),type='l',col='red')
plot(apply(x,1,sd),2^pathwayValues$E2F-1,xlim=c(0,max(2^pathwayValues$E2F-1)))
lines(rowMeans(x),rowMeans(x),type='l',col='red')

## comparing variability in negative binomial with and without reference
## for different values of uncertainty

sapply(1:1000,NBsim,mu=(2^pathwayValues$E2F-1)/(2^pathwayValues$E2F-1),
       ReferenceDataset=TumorNegExp[names(pathwayValues$E2F),],invChisq=T) -> x
sapply(1:1000,NBsim,mu=(2^pathwayValues$E2F-1)/(2^pathwayValues$E2F-1),ReferenceDataset=NULL,
       invChisq=T) -> y

plot(apply(x,1,sd),apply(y,1,sd),xlim=c(0,max(apply(rbind(x,y),1,sd))),
     ylim=c(0,max(apply(rbind(x,y),1,sd))),
     xlab='sd reference', ylab='sd no reference')
lines(c(0,max(apply(rbind(x,y),1,sd))),c(0,max(apply(rbind(x,y),1,sd))),type='l',col='red')

sapply(1:1000,NBsim,mu=(2^pathwayValues$E2F-1),
       ReferenceDataset=TumorNegExp[names(pathwayValues$E2F),],invChisq=T) -> x
sapply(1:1000,NBsim,mu=(2^pathwayValues$E2F-1),ReferenceDataset=NULL,
       invChisq=T) -> y

plot(apply(x,1,sd),apply(y,1,sd),xlim=c(0,max(apply(rbind(x,y),1,sd))),
     ylim=c(0,max(apply(rbind(x,y),1,sd))),
     xlab='sd reference', ylab='sd no reference')
lines(c(0,max(apply(rbind(x,y),1,sd))),c(0,max(apply(rbind(x,y),1,sd))),type='l',col='red')

## testing out the function on a mean data matrix with apply
mu.mat <- cbind(2^pathwayValues$E2F-1,0.5*(2^pathwayValues$E2F-1),0.5*(2^pathwayValues$E2F-1))
mu.mat[1:10,3] <- 1
mu.mat.sim.ref <- apply(mu.mat,2,NBsim,
                        ReferenceDataset=TumorNegExp[names(pathwayValues$E2F),])
mu.mat.sim <- apply(mu.mat,2,NBsim)
