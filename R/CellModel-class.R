## Class Definitions ##

#' @title CellModel
#' @description An S4 class that contains the output of a cell simulation,
#'      along with all the parameters used to generate the simulation
#' @slot cells A list object that describes the state of the cells
#'      at each time. The state representation depends on the type of
#'      model run, and is accessed by the function designed for each
#'      model type. 
#' @slot params A list object containing all parameters used in the model
#' @export
setClass('CellModel', representation(
    cells = "list",
    params = "list")
)

setClass('OffLatticeModel', contains='CellModel')
setClass('DrasdoHohmeModel', contains='OffLatticeModel')

## Constructor ##

createCellModel <- function(params, ...)
{
    type <- params[['modelType']]
    if (type == 'DrasdoHohme')
    {
        model <- new("DrasdoHohmeModel", cells = NULL, params = NULL)
        model <- processParameters(model, params, ...)
    }
    else
    {
        stop('invalid model type')
    }
    return (model)
}

## Generics ##

setGeneric('run',
    function(model)
    {
        standardGeneric('run')
    }
)

setGeneric('processParameters',
    function(model, params, ...),
    {
        standardGeneric('processParameters')
    }
)

setGeneric('timeToRow',
    function(model, time)
    {
        standardGeneric('timeToRow')
    }
)

setGeneric('getColumn',
    function(model, time, col)
    {
        standardGeneric('getColumn')
    }
)

setGeneric('getCellPhase',
    function(model, id)
    {
        standardGeneric('getCellPhase')
    }
)

setGeneric('getCoordinates',
    function(model, time)
    {
        standardGeneric('getCoordinates')
    }
)

setGeneric('getRadii',
    function(model, time)
    {
        standardGeneric('getRadii')
    }
)

setGeneric('getAxisLength',
    function(model, time)
    {
        standardGeneric('getAxisLength')
    }
)

setGeneric('getAxisAngle',
    function(model, time)
    {
        standardGeneric('getAxisAngle')
    }
)

setGeneric('getGrowthRates',
    function(model, time)
    {
        standardGeneric('getGrowthRates')
    }
)

setGeneric('getCellTypes',
    function(model, time)
    {
        standardGeneric('getCellTypes')
    }
)

setGeneric('getCycleLengths',
    function(model, time)
    {
        standardGeneric('getCycleLengths')
    }
)

setGeneric('getNumberOfCells',
    function(model, time)
    {
        standardGeneric('getNumberOfCells')
    }
)

setGeneric('getNumberOfNeighbors',
    function(model, time, index, radius)
    {
        standardGeneric('getNumberOfNeighbors')
    }
)

setGeneric('getDensity',
    function(model, time)
    {
        standardGeneric('getDensity')
    }
)

setGeneric('plotCells',
    function(model, time)
    {
        standardGeneric('plotCells')
    }
)

setGeneric('interactivePlot',
    function(model)
    {
        standardGeneric('interactivePlot')
    }
)

setGeneric('cellSummary',
    function(model, time)
    {
        standardGeneric('cellSummary')
    }
)

