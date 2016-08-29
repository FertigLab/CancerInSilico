#' \code{getParameters} get a named list of parameters in the model
#'
#' @param model A CellModel
#' @param fullDist whether or not to return full distribution of cycle length
#' @return a named list of parameters in the model
#' @examples
#' getParameters(runCancerSim(1,1))
#' @export
#' 

setGeneric("getParameters", function(model, fullDist=FALSE)
    standardGeneric("getParameters"))

#' \code{getParameters} get a named list of parameters in the model
#'
#' @param model A CellModel
#' @param fullDist whether or not to return full distribution of cycle length
#' @return a named list of parameters in the model
#' @examples
#' getParameters(runCancerSim(1,1))
#' @rdname CellModel-class
#' @export

setMethod("getParameters", "CellModel",

    function(model, fullDist=FALSE) {

        retDist = mean(.cycleLengthDist(model))

        if (fullDist) {
 
         retDist = .cycleLengthDist(model)

        }

        ret_val = list(

            initialNum = .initialNumCells(model),
            runTime = .runTime(model),           
            initialDensity = .initialDensity(model),           
            inheritGrowth = .inheritGrowth(model),           
            outputIncrement = .outputIncrement(model),           
            randSeed = .randSeed(model),           
            epsilon = .epsilon(model),           
            nG = .nG(model),
            timeIncrement = .timeIncrement(model),
            cycleLengthDist = retDist

        )           

        return(ret_val)

    }

)
