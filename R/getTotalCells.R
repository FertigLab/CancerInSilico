#' \code{getTotalCells} Plots the total number of cells vs time
#'
#' @param model A CellModel
#' @return the size of the cell population over time
#' @examples
#' getTotalCells(runModel(10,100))
#' @export
#' 
setGeneric("getTotalCells", function(model)
    standardGeneric("getTotalCells"))

#' \code{getTotalCells} Plots the total number of cells vs time
#'
#' @param model A CellModel
#' @return the size of the cell population over time
#' @examples
#' getTotalCells(runModel(10,100))
#' @export

setMethod("getTotalCells", "CellModel",

    function(model) {

        total_cells = c()
        radii = seq(3,ncol(model),6)
        for (t in 1:nrow(model)) {
            total_cells[t] = sum(model[t,radii]>0)
        }
		return(total_cells)
    }

)
