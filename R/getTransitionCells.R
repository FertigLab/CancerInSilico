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

        time_window = model@parameters

        if (phase == 'S') {

            
