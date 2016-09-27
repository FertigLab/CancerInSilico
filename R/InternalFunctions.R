timeToRow <- function(model, time) {

    return (floor(time / .recordIncrement(model)) + 1)

}

getDensity <- function(model,time) {

    coords <- getCoordinates(model, time)
    ycoords <- xcoords + 1
    radii <- getRadii(model, time)
    
    # farthest distance from (0,0) of cell
    d <- max(sqrt(coords[,1]**2 + coords[,2]**2))
       
    return(sum(radii ** 2) / (d ^ 2))

}

getNumberOfCells <- function(model, time) {

    return(sum(getRadii(model, time) > 0))
    
}

getNumNeighbors <- function(model, time, index, radius) {

    num <- 0 
    coords <- getCoordinates(model, time)
    
    for (i in setdiff(1:nrow(coords), index)) {

        if ((coords[i,1] - coords[index,1])^2 + (coords[i,2] - coords[index,2])^2 < radius^2) {

            num <- num + 1

        }

    }

    return (num)

}

plotCellsAtTime <- function(model,time)  {

    coords <- getCoordinates(model, time)
    radii <- getRadii(model, time)
    axis_len <- getAxisLength(model, time)
    axis_ang <- getAxisAngle(model, time)

    mn <- min(coords) - 2
    mx <- max(coords) + 2
    
    plot(c(mn,mx),c(mn,mx),main=paste("Plot of CellModel At Time",time),xlab = "",ylab="",type="n",asp=1)
          
    x_1 <- coords[,1] + (0.5 * axis_len - radii) * cos(axis_ang)
    x_2 <- coords[,1] - (0.5 * axis_len - radii) * cos(axis_ang)
    y_1 <- coords[,2] + (0.5 * axis_len - radii) * sin(axis_ang)
    y_2 <- coords[,2] - (0.5 * axis_len - radii) * sin(axis_ang)
    
    x <- c(x_1,x_2)
    y <- c(y_1,y_2)
    rad <- c(radii, radii)
    
    symbols(x,y, circles=rad, inches=FALSE, add=TRUE, bg="bisque4", fg="bisque4")

}



