library(methods)

#' @title Drug
#' @description describes the properties of a drug
#'
#' @slot name name of drug
#' @slot timeAdded the time at which this drug is added to the simulation
#' @slot cycleLengthEffect effect this drug has on cell cycle length
#' @export
setClass('Drug', slots = c(
    name = 'character',
    timeAdded = 'numeric',
    cycleLengthEffect = 'function'
))

setMethod('initialize', 'Drug',
    function(.Object, ...)
    {
        if (is.null(list(...)$name)) stop('missing name')
        .Object@timeAdded <- 0
        .Object@cycleLengthEffect <- function(type, len) {return(len)}
        .Object <- callNextMethod(.Object, ...)
        .Object
    }
)

setValidity('Drug',
    function(object)
    {
        if (length(object@name) == 0)
            "missing 'name'"
        else if (length(object@timeAdded) == 0)
            "missing 'timeAdded'"
        else if (is.null(object@cycleLengthEffect(0,0)))   
            "missing 'cycleLengthEffect'"
        else if (object@timeAdded < 0)
            "'timeAdded' cannot be negative"
        else
            TRUE
    }
)


