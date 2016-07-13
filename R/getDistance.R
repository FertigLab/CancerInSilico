setGeneric("plotInteractive", function(mat,time)
  standardGeneric("plotInteractive"))

setMethod("plotInteractive", "CellMatrix",
          
          function(mat){
            #d^2 = (x2 - xtempplace)^2 + (y2 - ytempplace)^2
            radii = seq(3,ncol(mat),6)
            numCells = sum(mat[time,radii])
            
            xcoords = seq(tempplace,ncol(mat),6)
            ycoords = xcoords + tempplace
            
            place = tempplace
            while(place<=nrow(mat)){
              tempCell = c(mat[place,xcoords[place]],mat[place,ycoords[place]],mat[place,radii[place]])
              tempplace = tempplace
              while(tempplace <= nrow(mat)){
                if(tempplace == place){
                  tempplace = tempplace + 1
                }
                placeCell = c(mat[tempplace,xcoords[tempplace]],mat[tempplace,ycoords[tempplace]],mat[tempplace,radii[tempplace]])
              }
            }
            
          }
          
)