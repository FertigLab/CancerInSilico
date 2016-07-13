#' \code{plotTotalCells} Plots the total number of cells vs time
#'
#' @param model A CellModel
#' @return Plots the size of the cell population over time
#' @examples
#' plotTotalCells(runModel(10,100))

setGeneric("plotTotalCells", function(model)
    standardGeneric("plotTotalCells"))

#' \code{plotTotalCells} Plots the total number of cells vs time
#'
#' @param model A CellModel
#' @return Plots the size of the cell population over time
#' @examples
#' plotTotalCells(runModel(10,100))
#' @export

setMethod("plotTotalCells", "CellModel",

    function(model) {

        total_cells = c()
        radii = seq(3,ncol(model),6)

        for (t in 1:nrow(model)) {

            total_cells[t] = sum(model[t,radii]>0)

        }

        plot(total_cells,main="Plot of Total Cells vs Time",xlab = "Time",ylab = "Number of Cells",type="l")

    }

)
