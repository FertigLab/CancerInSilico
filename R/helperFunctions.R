DrawCircle <- function(x,y,rad) {
    theta <- seq(0,2*pi,length=20)
    points(x + rad * cos(theta),y + rad * sin(theta))

}
