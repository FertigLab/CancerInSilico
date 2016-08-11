#' \code{simulateDrugEffect} Simulate Gene Expression Data (Per Cell)
#'
#' @param model A CellModel
#' @param pathway A list of pathways, Format:(GtoM, GtoS, Prox)
#' @return the size of the cell population over time
#' @export

setGeneric("simulateDrugEffect", function(model,pathway)
    standardGeneric("simulateDrugEffect"))

setMethod("simulateDrugEffect", "CellModel",
          
        function(model,pathway) {
            numdgenes = rexp(length(pathway),1/3)
            drmatrix = matrix(0,model@parameters[2],length(numdgenes))
            
            probs = runif(model@parameters[2] * length(numdgenes),0,1)
            mn = min(probs)
            mx = max(probs)
            
            drug = getDrugEffect(FUN = funtion(x)(return((x - mn)/(mx - mn))),cycleLengthDist = model@paramCycleLengthDist)
            for(t in drug[1]:model@parameters[2]){
                xcoords = seq(1,length(model@cells[[timeToRow(model,t)]]),6)
                numcells = sum(model@cells[[timeToRow(model,t)]][xcoords] > 0)
                dist = rnorm(length(numdgenes),drug[2],sqrt(sum()))
            }
            test = runif(model@parameters[2] * length(pathway),0,1)
            
        }
)