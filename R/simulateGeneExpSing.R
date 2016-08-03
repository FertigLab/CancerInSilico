#' \code{simulateGeneExpSing} Simulate Gene Expression Data (Each Cell)
#'
#' @param model A CellModel
#' @param pathway A list of pathways, Format:(GtoM, GtoS, Prox)
#' @param time Time Specified by user
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGeneExpSing", function(model,pathway)
    standardGeneric("simulateGeneExpSing"))

setMethod("simulateGeneExpSing", "CellModel",
          
          function(model,pathway) {
              
              nummgenes = rexp(length(pathway),1/3)
              output = list()
              for(t in 1:model@parameters[2]){
                  radii <- seq(3,length(model@cells[[timeToRow(model,t)]]),6)
                  #vector for the specific time
                  cells = matrix(0,length(radii),length(pathway))
                  rownames(cells,TRUE,prefix = "cell ")
                  colnames(cells)<-pathway
                  if(length(getCellPhasePos(model,t)) == 0){
                      #Case: None are dividing
                      output[[t]] = t(cells)
                  }
                  else{
                      #Case: Some cells are dividing
                      #Calculate the average of each gene at the time
                      x = rep(0,length(nummgenes))
                      cells[getCellPhasePos(model,t),] = nummgenes
                      output[[t]] = t(cells)
                  }
              }
              return(output)
          }
)