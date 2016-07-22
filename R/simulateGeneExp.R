#' \code{getTotalCells} Simulate Gene Expression Data (Beginning)
#'
#' @param model A CellModel
#' @param genes A vector of gene names
#' @return the size of the cell population over time
#' @examples
#' getTotalCells(runModel(10,100))
#' @export

setGeneric("simulateGeneExp", function(model,genes)
    standardGeneric("simulateGeneExp"))

setMethod("simulateGeneExp", "CellModel",
          
          function(model,genes) {
              #Randomized Values for # of Each Genes
              numogenes = rexp((genes),1/3) #rep(1,length(genes))
              #Total Number of Cells
              numcells = getTotalCells(model)
              
              #Matrix Creation
              geneMatrix = matrix(NA,length(model@cells),length(numogenes))
              for(t in 1:length(model@cells)){
                  #Genes Per Cell
                  gpcell = length(splitCells(model,t)) * numogenes
                  #Check if there are cells dividing at the time
                  if(length(splitCells(model,t)) == 0){
                      #Case: None are dividing
                      gpcell = rep(0,length(numogenes))
                      geneMatrix[t,] = gpcell
                  }
                  else{
                      #Case: Some cells are dividing
                      #Calculate the average of each gene at the time
                      geneMatrix[t,] = gpcell/numcells[t]
                  }
              }
              colnames(geneMatrix)<-genes
              rownames(geneMatrix)<-c(1:401)
              return(geneMatrix)
          }
)