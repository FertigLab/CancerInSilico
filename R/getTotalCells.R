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

        for (t in 1:length(model@cells)) {

            radii = seq(3, length(model@cells[[t]]), 6)
            total_cells[t] = sum(model@cells[[t]][radii] > 0)

        }

		return(total_cells)

    }

)
