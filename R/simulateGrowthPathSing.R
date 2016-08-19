#' \code{simulateGrowthPathSing} Simulate Gene Expression Data (Average)
#'
#' @param model A CellModel
#' @param pathway A vector of gene names
#' @param sampFreq Time (in hours) at which to simulate gene expression data
#' @param nCell Number of cells selected for a random sample at given time
#' @return the size of the cell population over time
#' @export


setGeneric("simulateGrowthPathSing", function(model,pathway,sampFreq = 1,ncell = 0)
    standardGeneric("simulateGrowthPathGroup"))

setMethod("simulateGrowthPathSing", "CellModel",
          
          function(model,pathway,sampFreq = 1,ncell = 0) {
              numfgenes = pathway
              times = seq(sampFreq,model@parameters[2],sampFreq)
              t = sampFreq
              gfmatrix = matrix(0,length(times),length(numfgenes))
              altout = matrix(NA, length(pathway), length(times) * model@cells[[times[length(times)]]])
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
