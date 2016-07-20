#' \code{splitCells} Determines which cells are in the process of dividing
#'
#' @param model A CellModel
#' @param time Time
#' @return the size of the cell population over time
#' @export

setGeneric("splitCells", function(model,time)
    standardGeneric("splitCells"))

setMethod("splitCells", "CellModel",
          #returns a vector of positions of the cells that have divided
          function(model,time) {
              positions = vector()
              axislen <- seq(4,ncol(model),6)
              axlen <- model[time,axislen]
              axtemp1 <- model[time+1,axislen]
              axtemp2 <- model[time+2,axislen]
              axtemp3 <- model[time+3,axislen]
              axtemp4 <- model[time+4,axislen]
              #Checks next 4 timesteps if a division will occur
              if(sum(axlen > axtemp1) == 1){
                  positions = which(axlen > axtemp1)
                  positions = positions[!duplicated(positions)]
              }
              if(sum(axlen > axtemp2) == 1){
                  positions = which(axlen > axtemp2)
                  positions = positions[!duplicated(positions)]
              }
              if(sum(axlen > axtemp3) == 1){
                  positions = which(axlen > axtemp3)
                  positions = positions[!duplicated(positions)]
              }
              if(sum(axlen > axtemp4) == 1){
                  positions = which(axlen > axtemp4)
                  positions = positions[!duplicated(positions)]
              }
              return(positions)
          }
          
)