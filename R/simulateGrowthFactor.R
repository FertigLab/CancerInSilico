#' \code{simulateGrowthFactor} Simulate Gene Expression Data (Average)
#'
#' @param model A CellModel
#' @param pathway A vector of gene names
#' @return the size of the cell population over time
#' @export


setGeneric("simulateGrowthFactor", function(model,pathway)
    standardGeneric("simulateGrowthFactor"))

setMethod("simulateGeneExpGroup", "CellModel",
          
        function(model,pathway) {
            
            numfgenes = rexp(length(pathway),1/3)
            gfmatrix = matrix(0,length(model@cells),length(numfgenes))
            
            for(t in 1:length(model@cells)){
                
                growFact <- seq(6,length(model@cells[[timeToRow(model,t)]]),6)
                factors <- model@cells[[timeToRow(model,t)]][growFact]
                
            }
            
            
            
        }
)