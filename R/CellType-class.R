#### class definition ####

#' @title CellType
#' @description An S4 class to represent the type of cell (and its consequent behavior) 
#'
#' @slot mType the type of cell as a character
#' @export

setClass("CellType", representation(
                        mType = "character" ))

#### getters (parameters) ####
#' nothing at the moment

#### getters (cell type data) ####

#' gets the cell type
#'
#' @return cell type

getCellType <- function(model) {

    return (model@mType)

}

)
