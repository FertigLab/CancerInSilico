setMethod('processParameters',
    signature('OffLatticeModel'),
    function(model),
    {
        # call function on base class
        base <- new('CellModel', cells = NULL, params = model@params)
        model@params <- processParameters(base)
        
        # check off lattice parameters
        
    }
)

setGeneric('getCellPhase',
    function(model, id)
    {
        standardGeneric('getCellPhase')
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

setGeneric('getNumberOfNeighbors',
    function(model, time, index, radius)
    {
        standardGeneric('getNumberOfNeighbors')
    }
)



