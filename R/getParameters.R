#' \code{getParameters} get a named list of parameters in the model
#'
#' @param model A CellModel
#' @param fullDist whether or not to return full distribution of cycle length
#' @return a named list of parameters in the model
#' @examples
#' getParameters(runModel(1,1))
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
#' getParameters(runModel(1,1))
#' @export
#' 

setMethod("getParameters", "CellModel",

    function(model, fullDist=FALSE) {

        len = length(model@parameters)
        retDist = mean(model@parameters[10:len])
        if (fullDist) {
          retDist = model@parameters[10:len]
        }
        ret_val = list(
            initialNum=model@parameters[1],
            runTime=model@parameters[2],           
            density=model@parameters[3],           
            inheritGrowth=model@parameters[4],           
            outputIncrement=model@parameters[5],           
            randSeed=model@parameters[6],           
            epsilon=model@parameters[7],           
            nG=model@parameters[8],
            cycleLengthDist=retDist
        )           

        return(ret_val)

    }

)
