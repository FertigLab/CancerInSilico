radii = seq(3,ncol(mat),6)
numCells = sum(mat[time,radii]>0)

#Information of cells based on Matrix Values
xcoords = seq(1,(numCells-1)*7,6)
ycoords = xcoords + 1
radii = ycoords + 1
axis_len = radii + 1
axis_ang = axis_len + 1

mn = min(min(mat[,xcoords]),min(mat[,ycoords])) - 2
mx = max(max(mat[,xcoords]),max(mat[,ycoords])) + 2


  AddCircle <- function(x,y,rad) {
    theta <- seq(0,2*pi,length=200)
    lines(x + rad * cos(theta),y + rad * sin(theta))
  }