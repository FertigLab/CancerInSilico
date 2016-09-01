#' \code{simulateGtoSPathGroup} Simulate G1 to Synthesis Phase Gene Expression (Average)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @param sampFreq Time (in hours) at which to simulate gene expression data
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGtoSPathGroup", function(model,pathway,sampFreq = 1)
    standardGeneric("simulateGtoSPathGroup"))

setMethod("simulateGtoSPathGroup", "CellModel",

    function(model,pathway,sampFreq = 1) {

        times <- seq(0,model@parameters[2] - sampFreq,sampFreq)
        gsMatrix <- matrix(0,length(times),length(pathway))
        colnames(gsMatrix) <- names(pathway)
        rownames(gsMatrix) <- times

        for (t in times) {

            cur_rad <- getRadii(model, t)
            next_rad <- getRadii(model, t + sampFreq)
            indices <- which(next_rad > sqrt(3/2) & cur_rad < sqrt(3/2))               

            gsMatrix[t,] = pathway * length(indices) / getNumberOfCells(model, t) 

        }

        return (gsMatrix)

    }

)
