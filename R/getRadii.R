setGeneric("getRadii", function(model, time)
  standardGeneric("getRadii"))

setMethod("getRadii", "CellModel",
          
	function(model, time) {

        indices <- seq(3,length(model@cells[[timetoRow(model,t)]]),6)
        return(radii <- model@cells[[timetoRow(model,t)]][radii])

	}
         
)
