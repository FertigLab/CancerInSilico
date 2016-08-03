#' \code{simulateGToSPathGroup} Simulate G1 to Synthesis Phase Gene Expression (Average)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGToSPathGroup", function(model,pathway)
    standardGeneric("simulateGToSPath"))

setMethod("simulateGToSPathGroup", "CellModel",
            function(model,pathway) {
                numsgenes = rexp(length(pathway),1/3)
                gsMatrix = matrix(NA,model@parameters[2],length(numsgenes))
                for(t in 1:model@parameters[2]){
                    radii <- seq(3,length(model@cells[[timeToRow(model,t)]]),6)
                    radius <- data@cells[[timeToRow(model,t)]][radii]
                    #Total Number of Cells
                    numcells = sum(model@cells[[timeToRow(model,t)]][radii] > 0)
                    #Range for determining a cell
                    toprange <- max(radius)/sqrt(2) + 0.1
                    botrange <- max(radius)/sqrt(2) - 0.1
                    test = subset(radius,radius<toprange & radius>botrange)
                    #Genes Per Cell
                    gscell = length(test) * numsgenes
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
                colnames(gsMatrix)<-pathway
                rownames(gsMatrix)<-c(1:model@parameters[2])
                return(gsMatrix)
            }
)