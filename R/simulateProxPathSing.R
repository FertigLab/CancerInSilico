#' \code{simulateProxPathSing} Simulate Proximity Based Gene Expression (Per Cell)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @return the size of the cell population over time
#' @export

setGeneric("simulateProxPathSing", function(model,pathway)
    standardGeneric("simulateProxPathSing"))

setMethod("simulateProxPathSing", "CellModel",
        function(model,pathway) {
            proxgenes = pathway
            output = list()
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
                #Matrix Calculation
                cells = matrix(0,length(xcoords),length(pathway))
                rownames(cells,TRUE,prefix = "cell ")
                colnames(cells)<-names(pathway)
                if(sum(tester/6) == 0){
                    x = rep(0,length(proxgenes))
                    output[[t]] = t(cells)
                }
                else{
                    cells[proxcheck,] = proxgenes
                    output[[t]] = t(cells)
                }
            }
            return(output)
        }
)