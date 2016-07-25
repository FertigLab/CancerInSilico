#' \code{simulateGeneExp} Simulate Gene Expression Data (Beginning)
#'
#' @param model A CellModel
#' @param pathway A list of pathways, Format:(GtoM, GtoS, Prox)
#' @return the size of the cell population over time
#' @examples
#' getTotalCells(runModel(10,100))
#' @export

setGeneric("simulateGeneExp", function(model,pathway)
    standardGeneric("simulateGeneExp"))

setMethod("simulateGeneExp", "CellModel",
          
          function(model,pathway) {
              
              output = list()
              gtm = pathway[[1]]
              gts = pathway[[2]]
              prox = pathway[[3]]
              #Total Number of Cells
              numcells = getTotalCells(model)
              
              

              #Randomized Values for # of Each Genes
              nummgenes = rexp((gtm),1/3) #rep(1,length(gtm))
              numsgenes = rexp((gts),1/3)
              #Matrix Creation
              gmMatrix = matrix(NA,length(model@cells),length(nummgenes))
              gsMatrix = matrix(NA,length(model@cells),length(numsgenes))
              for(t in 1:length(model@cells)){
                  #Genes Per Cell
                  gmcell = length(getCellPhasePos(model,t)) * nummgenes
                  gscell = length(getCellPhasePos(model,t,2)) * numsgenes
                  #Check if there are cells dividing at the time
                  if(length(getCellPhasePos(model,t)) == 0){
                      #Case: None are dividing
                      x = rep(0,length(nummgenes))
                      gmMatrix[t,] = x
                      gsMatrix[t,] = x
                  }
                  else{
                      #Case: Some cells are dividing
                      #Calculate the average of each gene at the time
                      gmMatrix[t,] = gmcell/numcells[t]
                      gsMatrix[t,] = gscell/numcells[t]
                  }
              }
              colnames(gmMatrix)<-gtm
              rownames(gmMatrix)<-c(1:401)
              
              colnames(gsMatrix)<-gts
              rownames(gsMatrix)<-c(1:401)
              
              
              
              output[["GtoMPath"]] = gmMatrix
              output[["GtoSPath"]] = gsMatrix
          }
)