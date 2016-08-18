#' \code{simulateGToSPathSing} Simulate G1 to Synthesis Phase Gene Expression (Per Cell)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGToSPathSing", function(model,pathway,sampFreq = 1)
    standardGeneric("simulateGToSPathSing"))

setMethod("simulateGToSPathSing", "CellModel",
            function(model,pathway,sampFreq = 1) {
                numsgenes = pathway
                times = seq(sampFreq,model@parameters[2],sampFreq)
                altout = matrix(NA,length(pathway),length(times) * model@cells[[times[length(times)]]])
                cnames = vector()
                count = 1
                t = sampFreq
                while(t < model@parameters[2]){
                    radii <- seq(3,length(model@cells[[timeToRow(model,t)]]),6)
                    currradius <- model@cells[[timeToRow(model,t)]][radii]
                    test = vector();
                    if(is(try(model@cells[[timeToRow(model,t+1)]][radii],TRUE),'try-error')==FALSE){
                        nextradius <- model@cells[[timeToRow(model,t+1)]][radii]
                        test = which(nextradius > sqrt(3/2) & currradius < sqrt(3/2))
                    }
                    #Matrix Calculation
                    cells = matrix(0,length(radii),length(pathway))
                    #Generate Column Names
                    test = paste("t",t,rownames(cells,FALSE,"c"))
                    cnames = append(cnames,test)
                    #Add to Matrix
                    if(length(test) == 0){
                        altout[,count:(count+length(radii)-1)] = t(cells)
                        count = count + length(radii)
                    }
                    else{
                        #Case: Some cells are in target range
                        #Calculate the average of each gene at the time
                        cells[test,] = numsgenes
                        altout[,count:(count+length(radii)-1)] = t(cells)
                        count = count + length(radii)
                    }
                    t = t + 1
                }
                altout <- altout[,colSums(is.na(altout))<nrow(altout)]
                rownames(altout) <- names(pathway)
                colnames(altout) <- cnames
                return(altout)
            }
)