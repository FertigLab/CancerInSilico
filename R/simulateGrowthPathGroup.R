#' \code{simulateGrowthPathGroup} Simulate Gene Expression Data (Average)
#'
#' @param model A CellModel
#' @param pathway A vector of gene names
#' @param sampFreq Time (in hours) at which to simulate gene expression data
#' @return the size of the cell population over time
#' @export


setGeneric("simulateGrowthPathGroup", function(model,pathway,sampFreq = 1)
    standardGeneric("simulateGrowthPathGroup"))

setMethod("simulateGrowthPathGroup", "CellModel",
          
        function(model,pathway,sampFreq = 1) {
            numfgenes = pathway
            times = seq(sampFreq,model@parameters[2],sampFreq)
            t = sampFreq
            gfmatrix = matrix(0,length(times),length(numfgenes))
            
            while(t < model@parameters[2]){
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
                t = t + sampFreq
            }
   
            colnames(gfmatrix) <- names(pathway)
            rownames(gfmatrix) <- times
            return(gfmatrix)
        }
)
