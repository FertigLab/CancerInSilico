#### class definition ####

#' @title Drug
#' @description An S4 class to describe the properties of a drug
#'
#' @slot name name of drug
#' @slot timeAdded the time at which this drug is added to the simulation
#' @export

setClass('Drug', representation(
                    name = 'character',
                    timeAdded = 'numeric',
                    cycleLengthEffect = 'function'))

newDrug <- function(name, timeAdded)
{
    newDrug <- new('Drug')
    newDrug@name <- name
    newDrug@timeAdded <- timeAdded
    newDrug@cycleLengthEffect <- cycleLengthEffect

    if (cycleLengthEffect == NULL)
    {
        cycleLengthEffect <- function(type, cycleLength, phase)
        {
            return (cycleLength)
        }
    }
}



