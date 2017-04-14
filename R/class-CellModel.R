library(methods)

################ Class Definition ################

#' @title CellModel
#' @description The top-level CellModel class. All other cell model
#'  classes inherit from this in some way
#' @slot cells A list object that describes the state of the cells
#'  at each time. The state representation depends on the type of
#'  model run, and is accessed by the function designed for each
#'  model type. 
#' @slot initialNum number of cells at time 0
#' @slot runTime number of model hours to run the simulation
#' @slot density initial density of cells
#' @slot boundary keep cells within circular boundary
#' @slot syncCycles start all cells in the beginning of interphase
#' @slot randSeed random seed used for both R and C++ functions
#' @slot outputIncrement how often simulation info is displayed
#' @slot recordIncrement how often cell info is recorded (controls size
#'  of resulting CellModel object
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

setMethod('initialize', 'CellModel',
    function(.Object, initialNum, runTime, density,
    boundary = 1, syncCycles = FALSE, randSeed = 0, 
    outputIncrement = 4, recordIncrement = 0.1, timeIncrement = 0.001,
    cellTypes = c(new('CellType', name='DEFAULT')), cellTypeInitFreq = c(1),
    ...)
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
        else if (object@density > 0.7)
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
        else
            TRUE
    }
)
        
##################### Generics ###################

#' @export
setGeneric('run', function(model)
    {standardGeneric('run')})

#' @export
setGeneric('getCellPhases', function(model, time)
    {standardGeneric('getCellPhases')})

#' @export
setGeneric('getCellTypes', function(model, time)
    {standardGeneric('getCellTypes')})

#' @export
setGeneric('getCycleLengths', function(model, time)
    {standardGeneric('getCycleLengths')})

#' @export
setGeneric('getNumberOfCells', function(model, time)
    {standardGeneric('getNumberOfCells')})

#' @export
setGeneric('getDensity', function(model, time)
    {standardGeneric('getDensity')})

#' @export
setGeneric('timeToRow', function(model, time)
    {standardGeneric('timeToRow')})

#' @export
setGeneric('getColumn', function(model, time, col)
    {standardGeneric('getColumn')})

#' @export
setGeneric('cellSummary', function(model, time)
    {standardGeneric('cellSummary')})

#' @export
setGeneric('plotCells', function(model, time)
    {standardGeneric('plotCells')})

#' @export
setGeneric('interactivePlot', function(model)
    {standardGeneric('interactivePlot')})

##################### Methods ####################

setMethod('interactivePlot', signature('CellModel'),
    function(model)
    {
        # help message
        helpMSG <- paste('Basic Commands:\n',
            'b ARG = back ARG timesteps (default ARG = 1)\n',
            'n ARG = forward ARG timesteps (default ARG = 1)\n',
            't ARG - jump to timestep ARG (default ARG = 1)\n',
            'd ARG - change default ARG for other commands\n',
            's = summary of cells\nq = quit\nh = basic command help\n')

        default_arg <- 1 # default arg value    
        quit <- FALSE; # option to quit the plot

        # main loop for interactive plot
        while (!quit)
        {
            # correct for invalid time
            if (time > model@runTime)
            {
                time <- model@runTime
            }
            else if (time < 0 || !is.numeric(time) || is.na(time))
            {
                time <- 0
            }

            plotCells(model, time) # plot cells at current time
            read <- readline()     # get keyboard input
            place <- unlist(gregexpr(" ", read))[1] # separate "command arg"
            if (place == -1) {place <- nchar(read)} # no arg given

            # parse command
            cmd <- gsub(" ", "", substring(read, 1, place))
            arg <- suppressWarnings(as.numeric(gsub(" ", "", 
                substring(read, place, nchar(read)))))
            
            cmds <- c("n","b","t","d","s","q","h") # possible commands
            if (is.na(arg)) {arg <- default_arg}   # set default if no arg

            # check if valid command
            if ((cmd %in% cmds))
            {
                # switch on which command was provided
                switch (match(cmd, cmds),
                    {time = time + arg},    # 'n' - increase time by arg
                    {time = time - arg},    # 'b' - decrease time by arg
                    {time = arg},           # 't' - go to time arg
                    {default_arg = arg},    # 'd' - set default arg
                    {print(cellSummary(model, time))}, # 's' - cell summary
                    {quit <- TRUE},         # 'q' - quit display
                    {print(helpMSG)}        # 'h' - display help command
                )
            }
            else
            {
                cat("Enter a valid command. Type \"h\" for further help\n")
            }
        }
    }
)

setMethod('cellSummary', signature('CellModel'),
    function(model, time)
    {
        return (paste('Cell Density: ', getDensity(model, time),
            '\nNumber of Cells: ', getNumerOfCells(model, time), '\n'))
    }
)
