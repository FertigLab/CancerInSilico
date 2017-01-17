#### class definition ####

#' @title Drug
#' @description An S4 class to describe the properties of a drug
#'
#' @slot name name of drug
#' @slot timeAdded the time at which this drug is added to the simulation
#' @export

setClass('Drug', representation(
                        name = 'character',
                        timeAdded = 'numeric'))

newDrug <- function(name, timeAdded)
{

    return (new('Drug', name = name, timeAdded = timeAdded))

}



