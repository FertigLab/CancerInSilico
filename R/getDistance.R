setGeneric("getDistance", function(mat)
  standardGeneric("getDistance"))

setMethod("getDistance", "CellMatrix",
          
          function(mat){
            #d^2 = (x2 - xtempplace)^2 + (y2 - ytempplace)^2
            radii = seq(3,ncol(mat),6)
            xcoords = seq(1,ncol(mat),6)
            ycoords = xcoords + 1
            
            place = 1
            while(place<=nrow(mat)){
              tempCell = c(mat[place,xcoords[place]],mat[place,ycoords[place]],mat[place,radii[place]])
              tempplace = 1
              while(tempplace <= nrow(radii)){
                if(tempplace == place & tempplace + 1 <= 100){
                  tempplace = tempplace + 1
                }
                placeCell = c(mat[tempplace,xcoords[tempplace]],mat[tempplace,ycoords[tempplace]],mat[tempplace,radii[tempplace]])
                distance = sqrt(((placeCell[1]-tempCell[1])**2) + ((placeCell[2]-tempCell[2])**2))
                radiisum = placeCell[3] + tempCell[3]
                print(paste(distance,radiisum))
                if(distance < radiisum){
                  return(FALSE)
                  break
                }
                else{
                  print(paste(place,tempplace))
                  tempplace = tempplace + 1
                }
              }
              place = place + 1
            }
            return(TRUE)
          }
          
)