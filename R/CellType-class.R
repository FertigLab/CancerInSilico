#### class definition ####

#' @title CellType
#' @description An S4 class to the properties of a cell type
#'
#' @slot type the name of the cell type
#' @slot size the relative size (volume) of the cell
#' @slot cycleLenth function that returns the average (target) cycle
#'      length of this cell
#' @slot drugEffect function that returns the effect a drug has on the
#'      the cell properties
#' @slot inheritCycleLength bool - inherit cycle length from parent
#' @slot inheritDrugEffect bool - inherit drug effect from parent
#' @export

setClass('CellType', representation(
                        type = 'character',
                        size = 'numeric',
                        cycleLength = 'function',
                        drugEffect = 'function',
                        inheritCycleLength = 'logical',
                        inheritDrugEffect = 'function'))

newCellType <- function(type, size = 1, cycleLength = NULL,
drugEffect = NULL, inheritCycleLength = FALSE, inheritDrugEffect = NULL)
{

    newType <- new('CellType')
    newType@type <- type
    newType@size <- size
    newType@cycleLength <- cycleLength
    newType@drugEffect <- drugEffect
    newType@inheritCycleLength <- inheritCycleLength
    newType@inheritDrugEffect <- inheritDrugEffect
    
    if (cycleLengthDist == NULL)
    {
        newType@cycleLength <- function() {return(48)}
    }

    if (drugEffect == NULL)
    {
        ## cell is list of info about cell
        newType@drugEffect <- function(drug, cell) {return(cell)}
    }

    if (inheritCycleLength == NULL)
    {
        newType@inheritCycleLength <- FALSE
    } 

    if (inheritDrugEffect == NULL)
    {
        newType@inheritDrugEffect <- function(drug) {return(FALSE)}
    }

    return (newType)

}


