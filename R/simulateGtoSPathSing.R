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

                times = seq(0,model@parameters[2],sampFreq)
                altout = matrix(NA,length(pathway),length(times) * model@cells[[times[length(times)]]])
                cnames = vector()
                count = 1

                for (t in times) {

                    cur_rad <- getRadii(model, t)
                    next_rad <- getRadii(model, t + sampFreq)
                    indices <- which(next_rad > sqrt(3/2) & cur_rad < sqrt(3/2))                    

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
