#' \code{simulateGrowthFactor} Simulate Gene Expression Data (Average)
#'
#' @param model A CellModel
#' @param pathway A vector of gene names
#' @return the size of the cell population over time
#' @export


setGeneric("simulateGrowthFactor", function(model,pathway)
    standardGeneric("simulateGrowthFactor"))

setMethod("simulateGeneExpGroup", "CellModel",
          
          function(model,pathway) {
              
              
              
          }
)