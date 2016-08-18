#' \code{simulateGtoMPathSing} Simulate G2 to Mitosis Phase Gene Expression (Per Cell)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @param sampFreq The frequency at which data is collected
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGtoMPathSing", function(model,pathway,sampFreq = 1)
    standardGeneric("simulateGtoMPathSing"))

setMethod("simulateGtoMPathSing", "CellModel",
        function(model,pathway,sampFreq = 1) {
            nummgenes = pathway
            times = seq(sampFreq,model@parameters[2],sampFreq)
            altout = matrix(NA,length(pathway),length(times) * model@cells[[times[length(times)]]])
            cnames = vector()
            count = 1
            t = sampFreq
            while(t < model@parameters[2]){
                radii <- seq(3,length(model@cells[[timeToRow(model,t)]]),6)
                #vector for the specific time
                cells = matrix(0,length(radii),length(pathway))
                #Generate Column Names
                test = paste("t",t,rownames(cells,FALSE,"c"))
                cnames = append(cnames,test)
                #Add to Matrix
                if(length(getCellPhasePos(model,t,sampFreq)) == 0){
                    #Case: None are dividing
                    altout[,count:(count+length(radii)-1)] = t(cells)
                    count = count + length(radii)
                }
                else{
                    #Case: Some cells are dividing
                    #Calculate the average of each gene at the time
                    cells[getCellPhasePos(model,t,sampFreq),] = nummgenes
                    altout[,count:(count+length(radii)-1)] = t(cells)
                    count = count + length(radii)
                }
                t = t + sampFreq
            }
            altout <- altout[,colSums(is.na(altout))<nrow(altout)]
            rownames(altout) <- names(pathway)
            colnames(altout) <- cnames
            return(altout)
        }
)
