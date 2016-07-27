#' \code{simulateGeneExpExt} Simulate Gene Expression Data (Beginning)
#'
#' @param model A CellModel
#' @param pathway A list of pathways, Format:(GtoM, GtoS, Prox)
#' @return the size of the cell population over time
#' @examples
#' getTotalCells(runModel(10,100))
#' @export

setGeneric("simulateGeneExpExt", function(model,pathway)
    standardGeneric("simulateGeneExpExt"))

setMethod("simulateGeneExpExt", "CellModel",
          
          function(model,pathway) {
              
                output = list()
                gtm = pathway[[1]]
                gts = pathway[[2]]
                prox = pathway[[3]]
                
                temp = 0
                for(t in 1:length(model@cells)){
                    radii <- seq(3,length(model@cells[[t]]),6)
                    if(temp < max(model@cells[[t]][radii])){
                        temp = max(model@cells[[t]][radii])
                    }
                }
                #Randomized Values for # of Each Genes
                nummgenes = rexp((gtm),1/3) #rep(1,length(gtm))
                numsgenes = rexp((gts),1/3)
                #Matrix Creation
                gmMatrix = matrix(NA,length(model@cells),length(nummgenes))
                gsMatrix = matrix(NA,length(model@cells),length(numsgenes))
                for(t in 1:length(model@cells)){
                    #Radii of Model
                    radii <- seq(3,length(model@cells[[t]]),6)
                    radius <- data@cells[[t]][radii]
                    #Total Number of Cells
                    numcells = sum(model@cells[[t]][radii] > 0)
                    #List of ranges for a cell to be in synthesis phase
                    toprange <- max(radius)/sqrt(2) + 0.1
                    botrange <- max(radius)/sqrt(2) - 0.1
                    test = subset(radius,radius<toprange & radius>botrange)
                    #Genes Per Cell
                    gmcell = length(getCellPhasePos(model,t)) * nummgenes
                    gscell = length(test) * numsgenes
                    #GTOM PATHWAY
                    #Check if there are cells dividing at the time
                    if(length(getCellPhasePos(model,t)) == 0){
                        #Case: None are dividing
                        x = rep(0,length(nummgenes))
                        gmMatrix[t,] = x
                    }
                    else{
                        #Case: Some cells are dividing
                        #Calculate the average of each gene at the time
                        gmMatrix[t,] = gmcell/numcells
                    }
                    #GTOS PATHWAY
                    if(length(test) == 0){
                        x = rep(0,length(numsgenes))
                        gsMatrix[t,] = x
                    }
                    else{
                        #Case: Some cells are in target range
                        #Calculate the average of each gene at the time
                        gsMatrix[t,] = gscell/numcells
                    }
                }
                # colnames(gmMatrix)<-gtm
                # rownames(gmMatrix)<-c(1:401)
                # 
                # colnames(gsMatrix)<-gts
                # rownames(gsMatrix)<-c(1:401)
                
                
                output[["GtoMPath"]] = gmMatrix
                output[["GtoSPath"]] = gsMatrix
                return(output)
          }
)