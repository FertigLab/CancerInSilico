#' \code{getParameters} get a named list of parameters in the model
#'
#' @param model A CellModel
#' @return a named list of parameters in the model
#' @examples
#' getParameters(runModel(10,100))
#' @export

setGeneric("getParameters", function(model)
    standardGeneric("getParameters"))

#' \code{getParameters} get a named list of parameters in the model
#'
#' @param model A CellModel
#' @return a named list of parameters in the model
#' @examples
#' getParameters(runModel(10,100))
#' @export

setMethod("getParameters", "CellModel",

    function(model) {

        ret_val = list(
            initialNum=model@parameters[1],
            runTime=model@parameters[2],           
            density=model@parameters[3],           
            cycleTimeDist=model@parameters[4],           
            inheritGrowth=model@parameters[5],           
            timeIncrement=model@parameters[6],           
            outputIncrement=model@parameters[7],           
            randSeed=model@parameters[8],           
            epsilon=model@parameters[9]
        )           

        return(ret_val)

    }

)
