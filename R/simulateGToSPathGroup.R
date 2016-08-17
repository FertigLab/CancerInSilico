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
                numsgenes = pathway
                gsMatrix = matrix(0,model@parameters[2],length(numsgenes))
                for(t in 1:model@parameters[2]){
                    radii <- seq(3,length(model@cells[[timeToRow(model,t)]]),6)
                    currradius <- model@cells[[timeToRow(model,t)]][radii]
                    test = vector();
                    if(is(try(model@cells[[timeToRow(model,t+1)]][radii],TRUE),'try-error')==FALSE){
                        nextradius <- model@cells[[timeToRow(model,t+1)]][radii]
                        test = which(nextradius > sqrt(3/2) & currradius < sqrt(3/2))
                    }
                    #Number of Cells
                    numcells = sum(model@cells[[timeToRow(model,t)]][radii] > 0)
                    
                    #Genes Per Cell
                    gscell = length(test) * numsgenes
                    
                    if(length(test) != 0){
                        #Case: Some cells are in target range
                        #Calculate the average of each gene at the time
                        gsMatrix[t,] = gscell/numcells
                    }
                    
                }
                colnames(gsMatrix)<-names(pathway)
                rownames(gsMatrix)<-c(1:model@parameters[2])
                return(gsMatrix)
            }
)
