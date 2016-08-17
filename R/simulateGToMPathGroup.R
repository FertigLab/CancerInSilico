#' \code{simulateGToMPathGroup} Simulate G2 to Mitosis Phase Gene Expression (Average)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGToMPathGroup", function(model,pathway)
    standardGeneric("simulateGToMPathGroup"))

setMethod("simulateGToMPathGroup", "CellModel",
        function(model,pathway) {
            nummgenes = pathway
            gmMatrix = matrix(0,model@parameters[2],length(nummgenes))
            for(t in 1:model@parameters[2]){
                radii <- seq(3,length(model@cells[[timeToRow(model,t)]]),6)
                #Total Number of Cells
                numcells = sum(model@cells[[timeToRow(model,t)]][radii] > 0)
                #Genes Per Cell
                gmcell = length(getCellPhasePos(model,t)) * nummgenes
                #GTOM PATHWAY
                if(length(getCellPhasePos(model,t)) != 0){
                    #Case: Some cells are dividing
                    #Calculate the average of each gene at the time
                    gmMatrix[t,] = gmcell/numcells
                }
            }
            colnames(gmMatrix)<-names(pathway)
            rownames(gmMatrix)<-c(1:model@parameters[2])
            #output = gmMatrix[rowSums(gmMatrix[,-1]) != 0,]
            return(gmMatrix)
        }
)