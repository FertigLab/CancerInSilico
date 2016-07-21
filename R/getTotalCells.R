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
        radii = seq(3,length(model@cells[[length(model@cells)]]),6)
        for (t in 1:length(model@cells)) {
            cellsatt = model@cells[[t]][radii]
            cellsatt = sum(length(cellsatt[!cellsatt %in% NA]))
            total_cells[t] = cellsatt
        }
		return(total_cells)
    }

)
