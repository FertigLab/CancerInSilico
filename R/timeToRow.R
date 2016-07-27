setGeneric("timeToRow", function(model, time)
  standardGeneric("timeToRow"))

setMethod("timeToRow", "CellModel",
          
	function(model, time) {

        return(floor(time / model@parameters[9]) + 1)

	}
         
)
