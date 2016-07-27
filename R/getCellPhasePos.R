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
                axislen <- seq(4,length(model@cells[[time]]),6)
                
                #Checks if can take a step, and then fills in positions of divisions
                if(is(try(model@cells[[time]][axislen],TRUE),'try-error')==FALSE){
                    axlen <- model@cells[[time]][axislen]
                }
                #STEP 1
                if(is(try(model@cells[[time+1]][axislen],TRUE),'try-error')==FALSE){
                    axtemp1 <- model@cells[[time+1]][axislen]
                    positions = c(positions,which(axlen > axtemp1))
                }
                #STEP 2
                if(is(try(model@cells[[time+2]][axislen],TRUE),'try-error')==FALSE){
                    axtemp2 <- model@cells[[time+2]][axislen]
                    positions = c(positions,which(axlen > axtemp2))
                }
                #STEP 3
                if(is(try(model@cells[[time+3]][axislen],TRUE),'try-error')==FALSE){
                    axtemp3 <- model@cells[[time+3]][axislen]
                    positions = c(positions,which(axlen > axtemp3))
                }
                #STEP 4
                if(is(try(model@cells[[time+4]][axislen],TRUE),'try-error')==FALSE){
                    axtemp4 <- model@cells[[time+4]][axislen]
                    positions = c(positions,which(axlen > axtemp4))
                }
                #Removes Duplicates
                positions = positions[!duplicated(positions)]
                    return(positions)
                
                
                # #GTOS PATHWAY
                # else if(opt == 2){
                #     localmax = rep(0,length(radii))
                #     intmax = 0
                #     t = time
                #     #Check if localmax is full or not
                #     while(sum(which(localmax>0)) != length(localmax) & t <= length(model@cells)){
                #         #Check in order to update localmax
                #         if(sum(length(which(localmax > 0))) < sum(length(model@cells[[t]][radii]))){
                #             #Location to update localmax
                #             temp = which(localmax < model@cells[[t]][radii])
                #             if((sum(localmax[temp])>0) == TRUE){
                #                 rem = which(localmax>0)
                #                 temp = which(!(temp==rem))
                #             }
                #             localmax[temp] = model@cells[[t]][radii[temp]]
                #         }
                #         t = t + 1
                #         print(localmax)
                #     }
                #     return(localmax)
                # }
              
          }
          
)