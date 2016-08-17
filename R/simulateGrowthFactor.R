#' \code{simulateGrowthFactor} Simulate Gene Expression Data (Average)
#'
#' @param model A CellModel
#' @param pathway A vector of gene names
#' @return the size of the cell population over time
#' @export


setGeneric("simulateGrowthFactor", function(model,pathway)
    standardGeneric("simulateGrowthFactor"))

setMethod("simulateGrowthFactor", "CellModel",
          
        function(model,pathway) {
            
            numfgenes = rexp(length(pathway),1/3)
            gfmatrix = matrix(0,length(model@cells),length(numfgenes))
            
            for(t in 1:length(model@cells)){
                #Get Model Cell Data
                radii <- seq(3,length(model@cells[[timeToRow(model,t)]]),6)
                radius <- model@cells[[timeToRow(model,t)]][radii]
                growRate <- seq(6,length(model@cells[[timeToRow(model,t)]]),6)
                rates <- model@cells[[timeToRow(model,t)]][growRate]
                #Total Number of Cells
                numcells = sum(model@cells[[timeToRow(model,t)]][radii] > 0)
                
                avgrates = sum(rates)/(length(rates)*max(rates))
                gfcell = avgrates * numfcells / numcells
                
                gfmatrix[t,] = gfcell
                
            }
            return(gfmatrix)
        }
)