#' \code{plotTotalCells} Plots the total number of cells vs time
#'
#' @param mat A Cell Matrix
#' @return Plots the size of the cell population over time
#' @examples
#' plotTotalCells(runModel(10,100))
#' @export

setGeneric("plotTotalCells", function(mat)
    standardGeneric("plotTotalCells"))

<<<<<<< HEAD
=======
#' \code{plotTotalCells} Plots the total number of cells vs time
#'
#'
#' @param mat A Cell Matrix
#' @return Plots the size of the cell population over time
#' @export

>>>>>>> 1f571cc86bdbc85a7bd5f5e6d24cf291a11c89a3
setMethod("plotTotalCells", "CellMatrix",

    function(mat) {

        total_cells = c()
        radii = seq(3,ncol(mat),6)

        for (t in 1:nrow(mat)) {

            total_cells[t] = sum(mat[t,radii]>0)

        }

        plot(total_cells,type="l")

    }

)
