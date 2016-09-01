setGeneric("getRadii", function(model, time)
  standardGeneric("getRadii"))

setMethod("getRadii", "CellModel",
          
	function(model, time) {

        indices <- seq(3,length(model@cells[[timetoRow(model,t)]]),6)
        return(model@m_cells[[timetoRow(model,t)]][indices])

	}
         
)

setGeneric("getAxisLength", function(model, time)
  standardGeneric("getAxisLength"))

setMethod("getAxisLength", "CellModel",
          
	function(model, time) {

        indices <- seq(4,length(model@cells[[timetoRow(model,t)]]),6)
        return(model@m_cells[[timetoRow(model,t)]][indices])

	}
         
)
