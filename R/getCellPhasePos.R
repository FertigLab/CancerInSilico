#' \code{getCellPhasePos} Determines which cells are in the process of dividing
#'
#' @param model A CellModel
#' @param time Time
#' @return the size of the cell population over time
#' @export

setGeneric("getCellPhasePos", function(model,time)
    standardGeneric("getCellPhasePos"))

setMethod("getCellPhasePos", "CellModel",
            #returns a vector of positions of the cells that have divided
            function(model,time) {
                positions = vector()
                axislen <- seq(4,length(model@cells[[timeToRow(model,time)]]),6)
                #Checks if can take a step, and then fills in positions of divisions
                #axlen <- model@cells[[timeToRow(model,time)]][axislen]
                #Checks the hour between current and next to see if mitosis occurs and at which cells
                if(is(try(model@cells[[timeToRow(model,time+1)]][axislen],TRUE),'try-error')==FALSE){
                    for(t in timeToRow(model,time):timeToRow(model,time+1)){
                        axtemp <- model@cells[[t]][axislen]
                        axtemp2 <- model@cells[[t-1]][axislen]
                        positions = c(positions,which(axtemp2 > axtemp))
                    }  
                }
                positions = positions[!duplicated(positions)]
                return(positions)
                
          }
          
)