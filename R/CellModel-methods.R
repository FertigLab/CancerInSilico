setMethod('run',
    signature('CellModel'),
    function(model)
    {
        model@cells <- NULL
        model@cells <- cppRunModel(model@params)
        return (model)
    }
)

setMethod('processParameters',
    signature('CellModel'),
    function(model, params, ...)
    {
        # for readability
        p <- model@params

        # check for missing parameters
        if (is.null(p[['initialNum']])) stop('missing initialNum')
        if (is.null(p[['runTime']])) stop('missing runTime')
        if (is.null(p[['density']])) stop('missing density')
        if (is.null(p[['boundary']])) stop('missing boundary')
        if (is.null(p[['syncCycles']])) stop('missing syncCycles')

        # set default parameters
        if (is.null(p[['outputIncrement']])) p[['outputIncrement']] <- 6
        if (is.null(p[['recordIncrement']])) p[['recpordIncrement']] <- 0.25

        # check parameter values
        if (p[['initialNum']] < 1) stop('invalid initialNum')
        if (p[['runTime']] <= 0) stop('invalid runTime')
        if (p[['density']] <= 0) stop('invalid density')
        if (p[['density']] > 0.7) stop('density too high, seeding too slow')
        if (p[['outputIncrement']] <= 0) stop('invalid outputIncrement')
        if (p[['recordIncrement']] <= 0) stop('invalid recordIncrement')

        # check cell types
        if (is.null(p[['cellTypes']]))
        {
            stop('no cell types provided')
        }
        else
        {
            if (length(p[['cellTypes']]) != length(p[['cellTypeInitFreq']]))
            {
                stop('cell type frequency size != cell type size')
            }

            if (sum(p[['cellTypeInitFreq']]) != 1)
            {
                stop('cell type frequency doesn\'t sum to 1')
            }
        }

        # check drugs
        if (!is.null(p[['drugs']]))
        {
            isDrug <- function(d) {return (class(d)[1] == 'Drug')}
            if (prod(sapply(p[['drugs']], isDrug)) == 0)
            {
                stop('not all drugs are of class "Drug"')
            }
        }

        # return parameters
        return (p)
    }
)

setMethod('interactivePlot',
    signature('CellModel'),
    function(model)
    {
        # help message
        helpMSG <- paste('Basic Commands:\n',
            'b ARG = back ARG timesteps (default ARG = 1)\n',
            'n ARG = forward ARG timesteps (default ARG = 1)\n',
            't ARG - jump to timestep ARG (default ARG = 1)\n',
            'd ARG - change default ARG for other commands\n',
            's = summary of cells\nq = quit\nh = basic command help\n')

        # default arg value    
        default_arg <- 1

        # option to quit the plot
        quit <- FALSE;

        # main loop for interactive plot
        while (!quit)
        {
            # correct for invalid time
            if (time > model@params[['runTime']])
            {
                time <- model@params[['runTime']]
            }
            else if (time < 0 || !is.numeric(time) || is.na(time))
            {
                time <- 0
            }

            # call internal function which plots cells at current time
            plotCells(model, time)
    
            # get keyboard input
            read <- readline()

            # find separator between command and arg
            place <- unlist(gregexpr(" ", read))[1]

            if (place == -1) {place <- nchar(read)}

            # parse command
            cmd <- gsub(" ", "", substring(read, 1, place))
            arg <- suppressWarnings(as.numeric(gsub(" ", "", 
                substring(read, place, nchar(read)))))
            
            # list of possible commands
            cmds <- c("n","b","t","d","s","q","h")

            # set to default if no command is present
            if (is.na(arg)) {arg <- default_arg}

            # check if valid command
            if ((cmd %in% cmds))
            {
                # switch on which command was provided
                switch (match(cmd, cmds),
                    # 'n' increase time by arg
                    {time = time + arg},
                    # 'b' decrease time by arg
                    {time = time - arg},
                    # 't' go to time arg
                    {time = arg},
                    # 'd' set default arg
                    {default_arg = arg},
                    # 's' get cell summary
                    {
                        print(cellSummary(model, time))
                    },
                    # 'q' quit display
                    {
                        quit <- TRUE
                    },
                    # display help command
                    {
                        print(helpMSG)
                    }
                )
            }
            else
            {
                cat("Enter a valid command. Type \"h\" for further help\n")
            }
        }
    }
)

setMethod('cellSummary',
    signature('CellModel'),
    function(model, time)
    {
        return (paste('Cell Density: ', getDensity(model, time),
            '\nNumber of Cells: ', getNumerOfCells(model, time), '\n'))
    }
)

