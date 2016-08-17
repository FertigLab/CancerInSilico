#' \code{simulateGToMPathSing} Simulate G2 to Mitosis Phase Gene Expression (Per Cell)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGToMPathSing", function(model,pathway)
    standardGeneric("simulateGToMPathSing"))

setMethod("simulateGToMPathSing", "CellModel",
        function(model,pathway) {
            nummgenes = pathway
            output = list()
            for(t in 1:model@parameters[2]){
                radii <- seq(3,length(model@cells[[timeToRow(model,t)]]),6)
                #vector for the specific time
                cells = matrix(0,length(radii),length(pathway))
                rownames(cells,TRUE,prefix = "cell ")
                colnames(cells)<-names(pathway)
                if(length(getCellPhasePos(model,t)) == 0){
                    #Case: None are dividing
                    output[[t]] = t(cells)
                }
                else{
                    #Case: Some cells are dividing
                    #Calculate the average of each gene at the time
                    cells[getCellPhasePos(model,t),] = nummgenes
                    output[[t]] = t(cells)
                }
            }
            return(output)
        }
)