#### class definition ####

#' @title Drug
#' @description An S4 class to describe the properties of a drug
#'
#' @slot mName name of drug
#' @export

setClass('Drug', representation(
                        name = 'character',
                        timeAdded = 'numeric'))

newDrug <- function(name, timeAdded) {

    return (new('Drug', name = name, timeAdded = timeAdded))

}



