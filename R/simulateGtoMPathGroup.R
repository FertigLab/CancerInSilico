#' \code{simulateGtoMPathGroup} Simulate G2 to Mitosis Phase Gene Expression (Average)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @param sampFreq Time (in hours) at which to simulate gene expression data
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGtoMPathGroup", function(model,pathway,sampFreq = 1)
    standardGeneric("simulateGtoMPathGroup"))

setMethod("simulateGtoMPathGroup", "CellModel",
        function(model,pathway,sampFreq = 1) {
            nummgenes = pathway
            times = seq(sampFreq,model@parameters[2],sampFreq)
            t = sampFreq
            gmMatrix = matrix(0,length(times),length(nummgenes))
            while(t < model@parameters[2]){
                radii <- seq(3,length(model@cells[[timeToRow(model,t)]]),6)
                #Total Number of Cells
                numcells = sum(model@cells[[timeToRow(model,t)]][radii] > 0)
                #Genes Per Cell
                gmcell = length(getCellPhasePos(model,t)) * nummgenes
                #GTOM PATHWAY
                if(length(getCellPhasePos(model,t)) != 0){
                    #Case: Some cells are dividing
                    #Calculate the average of each gene at the time
                    gmMatrix[t,] = gmcell/numcells
                }
                t = t + sampFreq
            }
            colnames(gmMatrix)<-names(pathway)
            rownames(gmMatrix)<-times
            return(gmMatrix)
        }
)
