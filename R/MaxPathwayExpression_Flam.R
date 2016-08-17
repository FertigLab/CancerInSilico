###CancerInSilico 8.10.2016
###Emily Flam
###Find sample with max median value of raw RNA expression and expression Z-score for genes in given pathway
###Generate heatmap

library(gplots)

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
TumorNegExp <- TumorExp[,which(colnames(TumorExp) %in% names(which(HPVpheno == "-")))]

#Create vector of desired gene names
#Random 7 genes chosen here for example 
GeneSetMatrix <- read.table("geneset-2.txt", header=T, sep = "\t")
GeneSet <- as.character(GeneSetMatrix[2:nrow(GeneSetMatrix),1])

#Create matrix of RNA expression values of just pathway genes for all HPV - tumor samples 
GeneSetExp <- TumorNegExp[which(rownames(TumorNegExp) %in% GeneSet),]

#Get Z-score for GeneSetExp data
GeneSetZ <- sweep(GeneSetExp,1,apply(GeneSetExp,1,max),FUN="/")

#Create Vector of median values of RNA expression of pathway genes for all samples
MedianVec <- numeric()
MedianZvec <- numeric()
for (i in 1:ncol(GeneSetExp)){
  #based on raw median
  MedianVec[i] <- median(GeneSetExp[,i])
  names(MedianVec)[i] <- colnames(GeneSetExp)[i]
  
  #based on Z score
  MedianZvec[i] <- median(GeneSetZ[,i])
  names(MedianZvec)[i] <- colnames(GeneSetZ)[i]
}

#Get sample with max Median Value for pathway genes (both raw and Z-score)
#Generate list of expression values of pathway genes for max sample
MaxSamp <- names(which(MedianVec == max(MedianVec)))
MaxSampZ <- names(which(MedianZvec == max(MedianZvec)))

#raw median
MaxMedian <- MedianVec[MaxSamp]
GeneList <- GeneSetExp[,MaxSamp]

#Z-score
MaxMedianZ <- MedianZvec[MaxSampZ]
GeneListZ <- GeneSetZ[,MaxSampZ]

#generate heatmap
ColorVec <- c(rep("aquamarine", ncol(GeneSetExp)))
names(ColorVec) <- colnames(GeneSetExp)
ColorVec[names(MaxMedianZ)] <- "darkorchid"
heatmap.2(GeneSetExp, Colv=T, Rowv=T, trace="none", scale='row',
          col=redblue, ColSideColors = ColorVec,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2))
