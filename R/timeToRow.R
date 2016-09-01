setGeneric("timeToRow", function(model, time)
  standardGeneric("timeToRow"))

setMethod("timeToRow", "CellModel",
          
	function(model, time) {

        return(floor(time / .timeIncrement(model)) + 1)

	}
         
)
