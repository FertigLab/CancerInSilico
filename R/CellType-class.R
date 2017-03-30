#### class definition ####

#' @title CellType
#' @description An S4 class to the properties of a cell type
#'
#' @slot type the name of the cell type
#' @slot size the relative size (volume) of the cell
#' @slot cycleLength function that returns sample from distribution of 
#'    cycle lengths
#' @slot inheritCycleLength bool - inherit cycle length from parent
#' @export

setClass('CellType', representation(
                        name = 'character',
                        size = 'numeric',
                        cycleLength = 'function',
                        inheritCycleLength = 'logical'))

newCellType <- function(name, size = 1, cycleLength = NULL,
inheritCycleLength = FALSE)
{
    newType <- new('CellType')
    newType@type <- type
    newType@size <- size
    newType@cycleLength <- cycleLength
    newType@inheritCycleLength <- inheritCycleLength
    
    if (cycleLengthDist == NULL)
    {
        newType@cycleLength <- function() {return(48)}
    }

    if (inheritCycleLength == NULL)
    {
        newType@inheritCycleLength <- FALSE
    } 

    return (newType)
}


