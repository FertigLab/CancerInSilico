#' \code{simulateProxPathSing} Simulate Proximity Based Gene Expression (Per Cell)
#'
#' @param model A CellModel
#' @param pathway A gene pathway
#' @param sampFreq Time (in hours) at which to simulate gene expression data
#' @return the size of the cell population over time
#' @export

setGeneric("simulateProxPathSing", function(model,pathway)
    standardGeneric("simulateProxPathSing"))

setMethod("simulateProxPathSing", "CellModel",
        function(model,pathway) {
            proxgenes = pathway
            times = seq(sampFreq,model@parameters[2],sampFreq)
            altout = matrix(NA,length(pathway),length(times) * model@cells[[times[length(times)]]])
            cnames = vector()
            count = 1
            t = sampFreq
            while(t < model@parameters[2]){
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
                #Generate Column Names
                test = paste("t",t,rownames(cells,FALSE,"c"))
                cnames = append(cnames,test)
                #Add to Matrix
                if(sum(tester/6) == 0){
                    altout[,count:(count+length(radii)-1)] = t(cells)
                    count = count + length(radii)
                }
                else{
                    cells[proxcheck,] = proxgenes
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