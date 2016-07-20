#' \code{getTotalCells} Simulate Gene Expression Data (Beginning)
#'
#' @param model A CellModel
#' @return the size of the cell population over time
#' @examples
#' getTotalCells(runModel(10,100))
#' @export

setGeneric("simulateGeneExp", function(model,time,genes)
    standardGeneric("simulateGeneExp"))

setMethod("simulateGeneExp", "CellModel",
          
          function(model,time,genes) {
              numogenes = rexp(length(genes),1/3)
              placeholder = seq(3,ncol(model),6)
              numcells = sum(model[time,placeholder]>0)
              
              #Check if there are cells dividing at the time
              if(length(splitCells(model,time)) == 0){
                  return(0)
              }
              else{
                  #Calculate the average of each gene at the time
                  gpcell = length(splitCells(model,time)) * numogenes
                  distocells = rnorm(numcells)
                  for(t in 1:ncol(gpcells)){
                      x = distrocells * gpcell[t]
                      gpcell[t] = sum(distrocells)/numcells
                  }
              }

          }
          
)