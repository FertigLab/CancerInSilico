#' \code{simulateGToSPathGroup} Simulate G1 to Synthesis Phase Gene Expression (Average)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGToSPathGroup", function(model,pathway)
    standardGeneric("simulateGToSPathGroup"))

setMethod("simulateGToSPathGroup", "CellModel",
            function(model,pathway) {
                numsgenes = rexp(length(pathway),1/3)
                gsMatrix = matrix(0,model@parameters[2],length(numsgenes))
                for(t in 1:model@parameters[2]){
                    radii <- seq(3,length(model@cells[[timeToRow(model,t-1)]]),6)
                    prevradius <- model@cells[[timeToRow(model,t-1)]][radii]
                    curradius <- model@cells[[timeToRow(model,t)]][radii]
                    #Total Number of Cells
                    numcells = sum(model@cells[[timeToRow(model,t)]][radii] > 0)
                    
                    test = which(curradius > sqrt(3/2) & prevradius < sqrt(3/2))
                    #Genes Per Cell
                    gscell = length(test) * numsgenes
                    
                    if(length(test) != 0){
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

radii <- seq(3,length(model@cells[[timeToRow(model,t-1)]]),6)
prev_radius <- model@cells[[timeToRow(model,t-1)]][radii]
cur_radius <- model@cells[[timeToRow(model,t)]][radii]
