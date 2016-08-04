#' \code{simulateGToSPathSing} Simulate G1 to Synthesis Phase Gene Expression (Per Cell)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGToSPathSing", function(model,pathway)
    standardGeneric("simulateGToSPathSing"))

setMethod("simulateGToSPathSing", "CellModel",
            function(model,pathway) {
                numsgenes = rexp(length(pathway),1/3)
                output = list()
                for(t in 1:model@parameters[2]){
                    radii <- seq(3,length(model@cells[[timeToRow(model,t)]]),6)
                    radius <- data@cells[[timeToRow(model,t)]][radii]
                    #Total Number of Cells
                    numcells = sum(model@cells[[timeToRow(model,t)]][radii] > 0)
                    #Range for determining a cell
                    toprange <- max(radius)/sqrt(2) + 0.1
                    botrange <- max(radius)/sqrt(2) - 0.1
                    test = which(radius<toprange & radius>botrange)
                    #Matrix Calculation
                    cells = matrix(0,length(radii),length(pathway))
                    rownames(cells,TRUE,prefix = "cell ")
                    colnames(cells)<-pathway
                    if(length(test) == 0){
                        output[[t]] = t(cells)
                    }
                    else{
                        #Case: Some cells are in target range
                        #Calculate the average of each gene at the time
                        cells[test,] = numsgenes
                        output[[t]] = t(cells)
                    }
                }
                return(output)
            }
)