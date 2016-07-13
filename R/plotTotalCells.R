#' \code{plotTotalCells} Plots the total number of cells vs time
#'
#' @param model A CellModel
#' @return plots the total number of cells over time
#' @examples
#' plotTotalCells(runModel(10,100))
#' @export

setGeneric("plotTotalCells", function(model)
    standardGeneric("plotTotalCells"))

setMethod("plotTotalCells", "CellModel",

    function(model) {

        total_cells = c()
        radii = seq(3,ncol(model),6)

        for (t in 1:nrow(model)) {

            total_cells[t] = sum(model[t,radii]>0)

        }

        plot(total_cells,main="Total Cells",type="l",xlab="Time",ylab="Total Cells")

    }

)
