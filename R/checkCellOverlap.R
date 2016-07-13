setGeneric("checkCellOverlap", function(model, time)
  standardGeneric("checkCellOverlap"))

setMethod("checkCellOverlap", "CellModel",
          
	function(model, time) {

		radii = seq(3,ncol(model),6)
        xcoords = radii[model[time,radii]>0] - 2

		x_1 <- model[time,xcoords] + (0.5 * model[time,xcoords+3] - model[time,xcoords+2]) * cos(model[time,xcoords+4])
	    x_2 <- model[time,xcoords] - (0.5 * model[time,xcoords+3] - model[time,xcoords+2]) * cos(model[time,xcoords+4])
	    y_1 <- model[time,xcoords+1] + (0.5 * model[time,xcoords+3] - model[time,xcoords+2]) * sin(model[time,xcoords+4])
	    y_2 <- model[time,xcoords+1] - (0.5 * model[time,xcoords+3] - model[time,xcoords+2]) * sin(model[time,xcoords+4])

	    x <- c(x_1,x_2)
	    y <- c(y_1,y_2)
	    rad <- c(model[time,xcoords+2], model[time,xcoords+2])

		for (i in 1:length(x_1)) {

			for (j in setdiff(c(i, i + length(x_1)), 1:length(x))) {

				dist_1 = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
				dist_2 = sqrt((x[i + length(x_1)] - x[j])^2 + (y[i + length(x_1)] - y[j])^2)
				
				if (dist_1 < rad[i] + rad[j] || dist_2 < rad[i + length(x_1)] + rad[j]) {

					return(FALSE)

				}

			}

		}

		return(TRUE)

	}
         
)
