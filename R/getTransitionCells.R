#' \code{getTransitionCells}
#'
#' @param model A CellModel
#' @param time time of interest
#' @param phase phase the cell is transitioning into (S or M)
#' @return the indices of cells making the transition
#' @export

setGeneric("getTransitionCells", function(model,time,phase)
    standardGeneric("getTransitionCells"))

setMethod("getTransitionCells", "CellModel",
    function(model,time,phase) {

        time_window = 0.05 # arbitrary for now

        if (phase == 'S') {

            cur_rad <- getRadii(model, t)
            next_rad <- getRadii(model, t + time_window)
            return (which(next_rad > sqrt(3/2) & cur_rad < sqrt(3/2)))

        } else if (phase == 'M') {

            cur_axis <- getAxisLength(model, t)
            next_axis <- getAxisLength(model, t + time_window)
            return (which(next_axis < cur_axis))

        } else {

            stop("getTransitionCells: invalid phase parameter")

        }

    }

)

            
