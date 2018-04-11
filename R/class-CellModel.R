library(methods)

################ Class Definition ################

#' CellModel
#' @description The top-level CellModel class. All other cell model
#' classes inherit from this in some way
#' @slot cells A list object that describes the state of the cells
#' at each time. The state representation depends on the type of
#' model run, and is accessed by the function designed for each
#' model type. 
#' @slot initialNum number of cells at time 0
#' @slot runTime number of model hours to run the simulation
#' @slot density initial density of cells
#' @slot boundary keep cells within circular boundary
#' @slot syncCycles start all cells in the beginning of interphase
#' @slot randSeed random seed used for both R and C++ functions
#' @slot outputIncrement how often simulation info is displayed
#' @slot recordIncrement how often cell info is recorded (controls size
#' of resulting CellModel object
#' @slot timeIncrement controls how fine the model timestep is
#' @slot cellTypes list of CellType objects used in the model
#' @slot cellTypeInitFreq initial frequency of cell types among cells
#' @slot drugs list of Drug objects used in the model
#' @export
setClass('CellModel', contains = 'VIRTUAL', slots = c(
    cells = 'list',
    initialNum = 'numeric',
    runTime = 'numeric',
    density = 'numeric',
    boundary = 'numeric',
    syncCycles = 'logical',
    randSeed = 'numeric',
    outputIncrement = 'numeric',
    recordIncrement = 'numeric',
    timeIncrement = 'numeric',
    cellTypes = 'list',
    cellTypeInitFreq = 'numeric',
    drugs = 'list'
))

#' Constructor for CellModel
#' @param .Object CellModel object
#' @param initialNum initial number of cells
#' @param runTime run time of the model in hours
#' @param density initial density of the cell population
#' @param boundary impose a physical boundary on the cells
#' @param syncCycles synchronization all cells to the same point in the cycle
#' @param randSeed random seed
#' @param outputIncrement how often (model hours) to print simulation status
#' @param recordIncrement how oftern (model hours) to record cell information
#' @param timeIncrement internal time step (model hours) used by the model
#' @param cellTypes list of CellType objects
#' @param cellTypeInitFreq initial proportions of all cell types
#' @param drugs list of Drug objects
#' @param ... model specific parameters
#' @return initialized cell model object
#' @importFrom methods callNextMethod
setMethod('initialize', 'CellModel',
    function(.Object, initialNum, runTime, density,
    boundary = 1, syncCycles = FALSE, randSeed = 0, 
    outputIncrement = 4, recordIncrement = 0.1, timeIncrement = 0.001,
    cellTypes = c(new('CellType', name='DEFAULT')), cellTypeInitFreq = c(1),
    drugs = list(), ...)
    {
        # store parameters, don't overwrite existing values
        if (!length(.Object@initialNum))
            .Object@initialNum <- initialNum
        if (!length(.Object@runTime))
            .Object@runTime <- runTime
        if (!length(.Object@density))
            .Object@density <- density
        if (!length(.Object@boundary))
            .Object@boundary <- boundary
        if (!length(.Object@syncCycles))
            .Object@syncCycles <- syncCycles
        if (!length(.Object@randSeed))
            .Object@randSeed <- randSeed
        if (!length(.Object@outputIncrement))
            .Object@outputIncrement <- outputIncrement
        if (!length(.Object@recordIncrement))
            .Object@recordIncrement <- recordIncrement
        if (!length(.Object@timeIncrement))
            .Object@timeIncrement <- timeIncrement
        if (!length(.Object@cellTypes))
            .Object@cellTypes <- cellTypes
        if (!length(.Object@cellTypeInitFreq))
            .Object@cellTypeInitFreq <- cellTypeInitFreq
        if (!length(.Object@drugs))
            .Object@drugs <- drugs

        # adjust time/record increment to divide runtime evenly
        .Object@timeIncrement <- .Object@runTime /
            ceiling(.Object@runTime / .Object@timeIncrement)
        .Object@recordIncrement <- .Object@runTime /
            ceiling(.Object@runTime / .Object@recordIncrement)

        # finish intialization, return object
        .Object <- callNextMethod(.Object, ...)
        return(.Object)
    }
)

setValidity('CellModel',
    function(object)
    {
        isValidCellType <- function(ct) {is(ct, 'CellType') &
            is.logical(validObject(ct, test=TRUE))}

        isValidDrug <- function(d) {is(d, 'Drug') &
            is.logical(validObject(d, test=TRUE))}
            
        if (length(object@initialNum) == 0)
            "missing 'initialNum'"
        else if (length(object@runTime) == 0)
            "missing 'runTime'"
        else if (length(object@density) == 0)
            "missing 'density'"
        else if (length(object@boundary) == 0)
            "missing 'boundary'"
        else if (length(object@syncCycles) == 0)
            "missing 'syncCycles'"
        else if (length(object@randSeed) == 0)
            "missing 'randSeed'"
        else if (length(object@outputIncrement) == 0)
            "missing 'outputIncrement'"
        else if (length(object@recordIncrement) == 0)
            "missing 'recordIncrement'"
        else if (length(object@timeIncrement) == 0)
            "missing 'timeIncrement'"
        else if (length(object@cellTypes) == 0)
            "missing 'cellTypes'"
        else if (length(object@cellTypeInitFreq) == 0)
            "missing 'cellTypeInitFreq'"
        else if (object@initialNum < 1)
            "'initialNum' must be at least 1"
        else if (object@runTime <= 0)
            "'runTime must be greater than zero"
        else if (object@density <= 0)
            "'density' must greater than zero"
        else if (object@density >= 0.5)
            "'density' too high for initial seeding process"
        else if (object@outputIncrement <= 0)
            "'outputIncrement' must be greater than zero"
        else if (object@recordIncrement <= 0)
            "'recordIncrement' must be greater than zero"
        else if (object@timeIncrement <= 0)
            "'timeIncrement' must be greater than zero"
        else if (prod(sapply(object@cellTypes, isValidCellType)) == 0)
            "not all cell types are valid"
        else if (length(object@cellTypes)!=length(object@cellTypeInitFreq))
            "cell type frequency size != cell type size"
        else if (sum(object@cellTypeInitFreq) != 1)
            "'cellTypeInitFreq' does not sum to 1"
        else if (length(object@drugs) > 0)
            if (prod(sapply(object@drugs, isValidDrug)) == 0)
                "not all drugs are valid"
        else if (length(object@drugs) > 64)
            "can\'t have more than 64 drugs in the simulation"
        else
            TRUE
    }
)
        
##################### Generics ###################

#' run a cell model
#' @export
#' @docType methods
#' @rdname run-methods
#'
#' @param model cell model object
#' @return cell model object with simulation info
#' @examples
#' data(SampleModels)
#' run(modDefault)
setGeneric('run', function(model)
    {standardGeneric('run')})

#' get number of cells in the model at a given time
#' @export
#' @docType methods
#' @rdname getNumberOfCells-methods
#'
#' @param model cell model object
#' @param time hour of the model to query
#' @return number of cells
#' @examples
#' data(SampleModels)
#' getNumberOfCells(modDefault, modDefault@runTime)
setGeneric('getNumberOfCells', function(model, time)
    {standardGeneric('getNumberOfCells')})

#' get density of the cell population at a given time
#' @export
#' @docType methods
#' @rdname getDensity-methods
#'
#' @param model cell model object
#' @param time hour of the model to query
#' @return density
#' @examples
#' data(SampleModels)
#' getDensity(modDefault, modDefault@runTime)
setGeneric('getDensity', function(model, time)
    {standardGeneric('getDensity')})

#' get phase of a cell at a given time
#' @export
#' @docType methods
#' @rdname getCellPhase-methods
#'
#' @param model cell model object
#' @param time hour of the model to query
#' @param cell id of cell to query
#' @return cell phase
#' @examples
#' data(SampleModels)
#' getCellPhase(modDefault, modDefault@runTime, 1)
setGeneric('getCellPhase', function(model, time, cell)
    {standardGeneric('getCellPhase')})

#' get type of a cell at a given time
#' @export
#' @docType methods
#' @rdname getCellType-methods
#'
#' @param model cell model object
#' @param time hour of the model to query
#' @param cell id of cell to query
#' @return cell type
#' @examples
#' data(SampleModels)
#' getCellType(modDefault, modDefault@runTime, 1)
setGeneric('getCellType', function(model, time, cell)
    {standardGeneric('getCellType')})

#' get cycle length of a cell at a given time
#' @export
#' @docType methods
#' @rdname getCycleLength-methods
#'
#' @param model cell model object
#' @param time hour of the model to query
#' @param cell id of cell to query
#' @return cycle length in hours
#' @examples
#' data(SampleModels)
#' getCycleLength(modDefault, modDefault@runTime, 1)
setGeneric('getCycleLength', function(model, time, cell)
    {standardGeneric('getCycleLength')})

#' get neighborhood density around a cell at a given time
#' @export
#' @docType methods
#' @rdname getLocalDensity-methods
#'
#' @param model cell model object
#' @param time hour of the model to query
#' @param cell id of cell to query
#' @param radius distance to search for neighboring cells
#' @return density
#' @examples
#' data(SampleModels)
#' getLocalDensity(modDefault, modDefault@runTime, 1, 3.3)
setGeneric('getLocalDensity', function(model, time, cell, radius)
    {standardGeneric('getLocalDensity')})

#' get distance between two cells
#' @export
#' @docType methods
#' @rdname getCellDistance-methods
#'
#' @param model cell model object
#' @param time hour of the model to query
#' @param cellA id of cell to query
#' @param cellB id of cell to query
#' @return distance between cellA and cellB
#' @examples
#' data(SampleModels)
#' getCellDistance(modDefault, modDefault@runTime, 1, 2)
setGeneric('getCellDistance', function(model, time, cellA, cellB)
    {standardGeneric('getCellDistance')})

#' summary of cell model at a given time
#' @export
#' @docType methods
#' @rdname cellSummary-methods
#'
#' @param model cell model object
#' @param time hour of the model to query
#' @return string containing summary of model
#' @examples
#' data(SampleModels)
#' cellSummary(modDefault, modDefault@runTime)
setGeneric('cellSummary', function(model, time)
    {standardGeneric('cellSummary')})

#' plot cell population at a given time
#' @export
#' @docType methods
#' @rdname plotCells-methods
#'
#' @param model cell model object
#' @param time hour of the model to query
#' @return plot
#' @examples
#' data(SampleModels)
#' plotCells(modDefault, modDefault@runTime)
setGeneric('plotCells', function(model, time)
    {standardGeneric('plotCells')})

#' plot the cell population and interactively scroll through time points
#' @export
#' @docType methods
#' @rdname interactivePlot-methods
#'
#' @param model cell model object
#' @return plot
setGeneric('interactivePlot', function(model)
    {standardGeneric('interactivePlot')})

##################### Methods ####################

#' @rdname interactivePlot-methods
#' @aliases interactivePlot
setMethod('interactivePlot', signature('CellModel'),
    function(model)
    {
        helpMSG <- paste('Basic Commands:\n',
            'b ARG = back ARG timesteps (default ARG = 1)\n',
            'n ARG = forward ARG timesteps (default ARG = 1)\n',
            't ARG - jump to timestep ARG (default ARG = 1)\n',
            'd ARG - change default ARG for other commands\n',
            's = summary of cells\n q = quit\n h = basic command help\n')

        time <- 0 # start time
        defaultArg <- 1 # default command line arg value    
        quit <- FALSE; # option to quit the plot

        # main loop for interactive plot
        while (!quit)
        {
            # correct for invalid time
            if (time > model@runTime)
                time <- model@runTime
            else if (time < 0 || !is.numeric(time) || is.na(time))
                time <- 0

            time <-model@recordIncrement*ceiling(time/model@recordIncrement)

            plotCells(model, time) # plot cells at current time
            read <- readline()     # get keyboard input
            place <- unlist(gregexpr(" ",read))[1] # separate "command arg"
            if (place == -1) {place <- nchar(read)} # no arg given

            # parse command
            cmd <- gsub(" ", "", substring(read, 1, place))
            arg <- suppressWarnings(as.numeric(gsub(" ", "", 
                substring(read, place, nchar(read)))))
            
            cmds <- c("n","b","t","d","s","q","h") # possible commands
            if (is.na(arg)) {arg <- defaultArg}   # set default if no arg

            # check if valid command
            if ((cmd %in% cmds))
            {
                # switch on which command was provided
                switch (match(cmd, cmds),
                    {time = time + arg},    # 'n' - increase time by arg
                    {time = time - arg},    # 'b' - decrease time by arg
                    {time = arg},           # 't' - go to time arg
                    {defaultArg = arg},    # 'd' - set default arg
                    {cat(cellSummary(model, time))}, # 's' - cell summary
                    {quit <- TRUE},         # 'q' - quit display
                    {cat(helpMSG)}        # 'h' - display help command
                )
            }
            else
            {
                cat("Enter a valid command. Type \"h\" for further help\n")
            }
        }
    }
)

#' @rdname cellSummary-methods
#' @aliases cellSummary
setMethod('cellSummary', signature(model='CellModel'),
    function(model, time)
    {
        return (paste('Cell Density: ', round(getDensity(model, time), 2),
            '\nNumber of Cells: ', getNumberOfCells(model, time), '\n'))
    }
)

#' @importFrom methods showDefault
setMethod('show', signature('CellModel'),
    function(object)
    {
        object@cells <- list()
        showDefault(object)
    }
)
