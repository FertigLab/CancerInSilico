#' \code{simulateGtoMPathSing} Simulate G2 to Mitosis Phase Gene Expression (Per Cell)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @param sampFreq Time (in hours) at which to simulate gene expression data
#' @param nCell Number of cells selected for a random sample at given time
#' @return the size of the cell population over time
#' @export

setGeneric("simulateGtoMPathSing", function(model,pathway,sampFreq = 1,ncell = 0)
    standardGeneric("simulateGtoMPathSing"))

setMethod("simulateGtoMPathSing", "CellModel",
        function(model,pathway,sampFreq = 1,ncell = 0) {
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
                tests = gsub(" ","",(paste("t",t,rownames(cells,FALSE,"c"))))
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
                    x = samp[which(samp %in% getCellPhasePos(model,t,sampFreq))]
                }
                else{
                    x = getCellPhasePos(model,t,sampFreq)
                }
                cnames = append(cnames,tests)
                
                #Add to Matrix
                if(length(getCellPhasePos(model,t,sampFreq)) == 0){
                    #Case: None are dividing
                    altout[,count:(count+nrow(cells)-1)] = t(cells)
                    count = count + nrow(cells)
                }
                else{
                    #Case: Some cells are dividing
                    #Calculate the average of each gene at the time
                    cells[x,] = nummgenes
                    altout[,count:(count+nrow(cells)-1)] = t(cells)
                    count = count + nrow(cells)
                }
                t = t + sampFreq
            }
            altout <- altout[,colSums(is.na(altout))<nrow(altout)]
            rownames(altout) <- names(pathway)
            colnames(altout) <- cnames
            return(altout)
        }
)
