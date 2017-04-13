#' @include class-CellModel.R class-OffLatticeModel.R
NULL

library(methods)

################ Class Definition ################

#' @title DrasdoHohmeModel
#' @description Implementation of an off-latice cell-based model
#'  based on the work in Drasdo, Hohme (2003)
#'
#' @slot nG number of monte carlo steps between each growth trial
#' @slot epsilon constant that controls the probability trails are accepted
#' @slot delta controls distance over which short range interactions occur
#' @export
setClass('DrasdoHohmeModel', contains = 'OffLatticeModel', slot = c(
    nG = 'numeric',
    epsilon = 'numeric',
    delta = 'numeric'
))

setMethod('initialize', 'DrasdoHohmeModel',
    function(.Object, nG = 28, epsilon = 10.0, delta = 0.2, ...)
    {
        # get supplied parameters
        .Object@nG <- nG
        .Object@epsilon <- epsilon
        .Object@delta <- delta

        # get min cycle length
        types <- list(...)$cellTypes
        if (is.null(types)) {types <- c(new('CellType', name = 'TEMP'))}
        minCycle <- min(sapply(types, slot, name = 'minCycle'))

        # calculate time increment
        t1 <- delta / (4 * nG * (4 - sqrt(2)))
        t2 <- delta * (minCycle - 1) / (8 * nG * (sqrt(2) - 1))
        .Object@timeIncrement <- min(t1, t2)

        # calculate off lattice parameters
        .Object@maxTranslation <- delta / 2
        .Object@maxRotation <- acos((16 + delta^2 - 4 * delta) / 16)
        .Object@maxDeformation <- 2 * min(t1,t2) * nG * (4 - sqrt(2))

        # finish intialization, return object
        .Object <- callNextMethod(.Object, ...)
        validObject(.Object)
        return(.Object)
    }
)

setValidity('DrasdoHohmeModel',
    function(object)
    {
        if (length(object@nG) == 0)
            "missing 'nG'"
        else if (length(object@epsilon) == 0)
            "missing 'epsilon'"
        else if (length(object@delta) == 0)
            "missing 'delta'"
        else if (object@nG < 1)
            "'nG' must be at least one"
        else if (object@epsilon <= 0)
            "'epsilon' must be greater than zero"
        else if (object@delta <= 0)
            "'delta' must be greater than zero"
    }
)

##################### Methods ####################

#setMethod('run', signature('DrasdoHohmeModel'),
#    function(model)
#    {
#        model@cells <- NULL
#        model@cells <- cppRunModel(model, 'DrasdoHohme')
#        return(model)
#    }
#)


