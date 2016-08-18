#' \code{getCellPhasePos} Determines which cells are in the process of dividing
#'
#' @param model A CellModel
#' @param time Time
#' @param sampFreq Time (in hours) at which to simulate gene expression data
#' @return the size of the cell population over time
#' @export

setGeneric("getCellPhasePos", function(model,time,sampFreq)
    standardGeneric("getCellPhasePos"))

setMethod("getCellPhasePos", "CellModel",
            #returns a vector of positions of the cells that have divided
            function(model,time,sampFreq) {
                positions = vector()
                axislen <- seq(4,length(model@cells[[timeToRow(model,time)]]),6)
                #Checks if can take a step, and then fills in positions of divisions
                #axlen <- model@cells[[timeToRow(model,time)]][axislen]
                #Checks the hour between current and next to see if mitosis occurs and at which cells
                if(is(try(model@cells[[timeToRow(model,time+sampFreq)]][axislen],TRUE),'try-error')==FALSE){
                    for(t in timeToRow(model,time):timeToRow(model,time+sampFreq)){
                        axtemp <- model@cells[[t]][axislen]
                        axtemp2 <- model@cells[[t-sampFreq]][axislen]
                        positions = c(positions,which(axtemp2 > axtemp))
                    }  
                }
                positions = positions[!duplicated(positions)]
                return(positions)
                
          }
          
)