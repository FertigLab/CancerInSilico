#' \code{plotTotalCells} Plots the total number of cells vs time
#'
#' @param mat A CellModel
#' @return Plots the size of the cell population over time
#' @examples
#' plotTotalCells(runModel(10,100))

setGeneric("plotTotalCells", function(mat)
    standardGeneric("plotTotalCells"))

#' \code{plotTotalCells} Plots the total number of cells vs time
#'
#' @param mat A CellModel
#' @return Plots the size of the cell population over time
#' @examples
#' plotTotalCells(runModel(10,100))
#' @export

setMethod("plotTotalCells", "CellMatrix",

    function(mat) {

        total_cells = c()
        radii = seq(3,ncol(mat),6)

        for (t in 1:nrow(mat)) {

            total_cells[t] = sum(mat[t,radii]>0)

        }

        plot(total_cells,main="Plot of Total Cells vs Time",xlab = "Time",ylab = "Number of Cells",type="l")

    }

)
