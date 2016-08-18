#' \code{simulateGrowthPathGroup} Simulate Gene Expression Data (Average)
#'
#' @param model A CellModel
#' @param pathway A vector of gene names
#' @param samplingFrequency Time (in hours) at which to simulate gene expression data
#' @return the size of the cell population over time
#' @export


setGeneric("simulateGrowthPathGroup", function(model,pathway,samplingFrequency)
    standardGeneric("simulateGrowthPathGroup"))

setMethod("simulateGrowthPathGroup", "CellModel",
          
        function(model,pathway,samplingFrequency) {
            
            numfgenes = pathway 
            gfmatrix = matrix(0,model@parameters[2],length(numfgenes))
            
            for(t in 1:model@parameters[2]){
                #Get Model Cell Data
                radii <- seq(3,length(model@cells[[timeToRow(model,t)]]),6)
                radius <- model@cells[[timeToRow(model,t)]][radii]
                growRate <- seq(6,length(model@cells[[timeToRow(model,t)]]),6)
                rates <- model@cells[[timeToRow(model,t)]][growRate]
                #Total Number of Cells
                numcells = sum(model@cells[[timeToRow(model,t)]][radii] > 0)
                
                avgrates = sum(rates)/(length(rates)*max(rates))
                gfcell = avgrates * numfgenes / numcells
                
                gfmatrix[t,] = gfcell
                
            }
   
            colnames(gfmatrix) <- names(pathway)

            return(gfmatrix)
        }
)
