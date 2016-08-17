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
                    currradius <- model@cells[[timeToRow(model,t)]][radii]
                    test = vector();
                    if(is(try(model@cells[[timeToRow(model,t+1)]][radii],TRUE),'try-error')==FALSE){
                        nextradius <- model@cells[[timeToRow(model,t+1)]][radii]
                        test = which(nextradius > sqrt(3/2) & currradius < sqrt(3/2))
                    }
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