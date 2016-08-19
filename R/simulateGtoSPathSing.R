#' \code{simulateGtoSPathSing} Simulate G1 to Synthesis Phase Gene Expression (Per Cell)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @param sampFreq Time (in hours) at which to simulate gene expression data
#' @param nCell Number of cells selected for a random sample at given time
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGtoSPathSing", function(model,pathway,sampFreq = 1,ncell = 0)
    standardGeneric("simulateGtoSPathSing"))

setMethod("simulateGtoSPathSing", "CellModel",
            function(model,pathway,sampFreq = 1,ncell = 0) {
                numsgenes = pathway
                times = seq(sampFreq,model@parameters[2],sampFreq)
                altout = matrix(NA,length(pathway),length(times) * model@cells[[times[length(times)]]])
                cnames = vector()
                count = 1
                t = sampFreq
                while(t < model@parameters[2]){
                    radii <- seq(3,length(model@cells[[timetoRow(model,t)]]),6)
                    currradius <- model@cells[[timetoRow(model,t)]][radii]
                    test = vector();
                    if(is(try(model@cells[[timetoRow(model,t+1)]][radii],TRUE),'try-error')==FALSE){
                        nextradius <- model@cells[[timetoRow(model,t+1)]][radii]
                        test = which(nextradius > sqrt(3/2) & currradius < sqrt(3/2))
                    }
                    #Matrix Calculation
                    cells = matrix(0,length(radii),length(pathway))
                    #Generate Column Names
                    tests = paste("t",t,rownames(cells,FALSE,"c"))
                    if(ncell != 0){
                        if(length(radii)< ncell){
                            ncell2 = length(radii)
                        }
                        else{
                            ncell2 = ncell
                        }
                        samp = sample(1:length(radii),ncell2)
                        tests = tests[samp]
                        cells = matrix(0,length(samp),length(pathway))
                        x = samp[which(samp %in% test)]
                    }
                    else{
                        x = test
                    }
                    cnames = append(cnames,tests)
                    #Add to Matrix
                    if(length(x) == 0){
                        altout[,count:(count+nrow(cells)-1)] = t(cells)
                        count = count + nrow(cells)
                    }
                    else{
                        #Case: Some cells are in target range
                        #Calculate the average of each gene at the time
                        cells[x,] = numsgenes
                        altout[,count:(count+nrow(cells)-1)] = t(cells)
                        count = count + nrow(cells)
                    }
                    t = t + 1
                }
                altout <- altout[,colSums(is.na(altout))<nrow(altout)]
                rownames(altout) <- names(pathway)
                colnames(altout) <- cnames
                return(altout)
            }
)
