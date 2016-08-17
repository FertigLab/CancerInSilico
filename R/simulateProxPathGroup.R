#' \code{simulateProxPathGroup} Simulate Proximity Based Gene Expression (Average)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @return the size of the cell population over time
#' @export

setGeneric("simulateProxPathGroup", function(model,pathway)
    standardGeneric("simulateProxPathGroup"))

setMethod("simulateProxPathGroup", "CellModel",
          function(model,pathway) {
              proxgenes = pathway
              proxMatrix = matrix(0,model@parameters[2],length(proxgenes))
              for(t in 1:model@parameters[2]){
                  xcoords = seq(1,length(model@cells[[timeToRow(model,t)]]),6)
                  ycoords = xcoords + 1
                  radii = xcoords + 2
                  #Total Number of Cells
                  numcells = sum(model@cells[[timeToRow(model,t)]][radii] > 0)
                  #Getting distances for cells
                  x = model@cells[[timeToRow(model,t)]][xcoords]
                  y = model@cells[[timeToRow(model,t)]][ycoords]
                  rad = model@cells[[timeToRow(model,t)]][radii]
                  pairs = cbind(x,y)
                  #Rows are Cells, Columns are the pairing to a different cell
                  test = data.matrix(dist(pairs))
                  test = t(t(test) - rad) - rad
                  tester = apply(test<0.2 & test>0,1,sum)
                  #Genes Per Cell
                  proxcheck = which(tester<6)
                  proxcell = sum(tester[proxcheck]/6) * proxgenes
                  #PROXIMITY PATHWAY
                  if(sum(tester/6) != 0){
                      proxMatrix[t,] = proxcell/numcells
                  }
              }
              colnames(proxMatrix) <- names(pathway)
              rownames(proxMatrix) <- c(1:model@parameters[2])
              return(proxMatrix)
          }
)