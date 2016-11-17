#' \code{simulatePathway}
#'
#' Simulates gene expression data for a given pathway
#'
#' @param model A CellModel object
#' @param pathway List of genes and range of expression values (fields:
#'      'genes', 'max', 'min')
#' @param type Cellular function related to pathway ('S', 'M', 'PROX',
#'      'GROWTH')
#' @param sampFreq Time (hrs) between each sample
#' @param sampSize Size of sample for single cell data
#' @param singleCell T/F: generate single cell data
#' @param timeWindow Length of window examined to see if relevant cell
#'      state changed (for 'S' and 'M' pathways)
#' @param downReg T/F: pathway is down regulated by cell activity (type)
#' @return Gene expression matrix for given pathway 
#'
simulatePathway <- function(args) {

    # if not single cell data, sample size is "1" (only return mean)
    if (!args[['singleCell']]) { args[['nCells']] <- 1 }

    # find closest, valid, sampling frequency
    args[['sampFreq']] <- .recordIncrement(args[['model']]) *
            ceiling(args[['sampFreq']] / .recordIncrement(args[['model']]))

    # vector of times to get gene expression for
    times <- seq(0, .runTime(args[['model']]) - 
                    args[['timeWindow']], args[['sampFreq']])

    # create return matrix
    gsMatrix <- matrix(nrow = length(pathway[["genes"]]),
                       ncol = length(times) * args[['sampSize']])
    colnames(gsMatrix) <- rep(times, each = args[['sampSize']])
    rownames(gsMatrix) <- pathway[["genes"]]

    # loop through each time
    for (i in 1:length(times)) {

        # get a vector of all cells
        cells <- 1:getNumberOfCells(args[['model']], times[i])

        # if doing single cell, get sample of cells
        if (args[['singleCell']]) {
    
            cells <- sort(sample(cells, args[['sampSize']]))

        }

        # get scaling factor for gene expression
        scalingFactor <- getScalingFactor(args[['model']], cells, times[i],
            timeWindow, type)

        # rows to replace
        cols <- (args[['sampSize']] * (i - 1) + 1):(args[['sampSize']] * i)

        # get gene expresssion
        gsMatrix[,cols] <- getExpression(scalingFactor, args)

    }

    # return the gene expression matrix
    return (gsMatrix)

}

# get scaling factor related to a cellular function
getScalingFactor <- function(model, cells, t, timeWindow, type) {

    # switch on type of cellular function
    switch (type,

        # pathway is related to the G to S transition in the cell cycle
        GtoS =  { return (scalingFactorGtoS(model, cells, t, timeWindow)) },

        # pathway is related to the G to M transition in the cell cycle
        GtoM = { return (scalingFactorGtoM(model, cells, t, timeWindow)) },

        # pathway is related to the growth rate of a cell
        Growth = { return (scalingFactorGROWTH(model, cells, t)) },

        # pathway is related to number of neighboring cells
        Prox = { return (scalingFactorPROX(model, cells, t)) }

    )

}

# get gene expression, given a pathway and a scaling factor
getExpression <- function(scalingFactor, args) {

    # get range of expression values
    range <- args[['pathway']][["max"]] - args[['pathway']][["min"]]

    # multiply pathway by scaling factor
    if (args[['singleCell']]) {

        # multiply each cells value by each gene in the pathway
        return (range %*% t(scalingFactor) + args[['pathway']][["min"]])

    } else {

        # average the gene expression across cells
        return (mean(scalingFactor) * range + args[['pathway']][["min"]])

    }

}

# calculate gene expression for G to S related pathway
scalingFactorGtoS <- function(model, cells, time, time_window) {

    # get the axis lengths of each cell at current time 
    cur_rad <- getRadii(model, time)[cells]

    # get the axis lengths of each cell at end of time window
    next_rad <- getRadii(model, time + time_window)[cells]

    # return logical vector (0's and 1's), 1 if cell made the G to S
    # transition in this window. G to S transition is defined by a cell
    # crossing the halfway point of its growth
    return (next_rad > sqrt(3/2) & cur_rad < sqrt(3/2))

}

# calculate gene expression G to M related pathway 
scalingFactorGtoM <- function(model, cells, time, time_window) {

    # get the axis lengths of each cell at current time 
    cur_ax <- getAxisLength(model, time)[cells]

    # get the axis lengths of each cell at end of time window
    next_ax <- getAxisLength(model, time + time_window)[cells]

    # return logical vector (0's and 1's), 1 if cell made the
    # G to M transition in this window
    return (next_ax < cur_ax)

}

# calculate gene expression for a pathway effected by neighboring cells
scalingFactorPROX <- function(model, cells, time) {

    # create empty vector
    exp <- c()

    # loop through each cell
    for (i in cells) {

        # expression = number of neighboring cells (max 6) divided by 6
        exp[i] <- getNumNeighbors(model, time, i) / 6

    }

    # return expression vector
    return (exp)

}

# calculate gene expression for a pathway effected by growth rate
scalingFactorGROWTH <- function(model, cells, time) {

    # get the cycle lengths of each cell
    cycle_len <- getCycleLengths(model,time)[cells]

    # this function maps cycle length to [0,1] in a smooth way
    return (1 - 1 / (1 + exp(-0.25 * (cycle_len - 16))))

}

# Determine the pathways from which to simulate gene expression data
getPathways <- function(args) {

    # possible pathway types
    validTypes <- c('GtoM', 'GtoS', 'Prox', 'Growth')
  
    # get pathway names
    pwy_names <- names(args[['pathways']])

    # if no user pathways provided
    if (is.null(args[['pathways']])) {

        # use default pathways from the package 
        validPwys <- inSilicoPathways

    } else {

        # keep only valid pathways
        validPwys <- args[['pathways']][pwy_names %in% validTypes]
    
    }

    # check that at least 1 pathway is valid
    if (length(validPwys) == 0) {

        stop(paste('No valid pathways provided. Valid choices are: ',
                 paste(validTypes, collapse=', ')))

    # warn user if some of the provided pathways were not valid
    } else if (length(validPwys) < length(pwy_names)) {

        warning(paste('The following pathways are invalid and will not ',
            'be simulated: ', pwy_names[!(pwy_names %in% names(validPwys))]
        ))

    }

    return (validPwys)

}

# Determine the range of gene expression values to use for each pathway
setPathwayExpressionRange <- function(args) {

    # if no reference data set provided, sample values
    if (is.null(args[['ReferenceDataSet']])) {
    
        # iterate through each pathway
        for (pwy in names(args[['pathways']])) {

            # get number of genes in pathway
            pwy_len <- length(args[['pathways']][[pwy]][["genes"]])

            # sample min and max (min + diff) values from exponential
            args[['pathways']][[pwy]][["min"]] <- rexp(pwy_len, 1) + 3
            args[['pathways']][[pwy]][["max"]] <- 
                    args[['pathways']][[pwy]][["min"]] +
                        rexp(pwy_len, 1.9) + 0.75

        }

    } else {

        # get the names of genes in the pathways
        gene_names <- unique(unname(unlist(args[['pathways']])))

        # check that the reference dataset is valid
        checkReferenceDataSet(args[['ReferenceDataSet']], gene_names)

        # iterate through each pathway
        for (pwy in names(args[['pathways']])) {

            # get min/max genes in pathway for each sample
            D_path <- args[['ReferenceDataSet']][
                            args[['pathways']][[pwy]][["genes"]], ]

            # get min/max expression values           
            args[['pathways']][[pwy]][["min"]] <-
                                    unname(apply(D_path, 1, min))
            args[['pathways']][[pwy]][["max"]] <-
                                    unname(apply(D_path, 1, max))

        }

    }

    return (args[['pathways']])

}

# verify reference data set
checkReferenceDataSet <- function(ReferenceDataSet, genes) {

    # data set contains all neccesary genes
    if (length(setdiff(genes, row.names(ReferenceDataSet))) > 0) {

        stop('ReferenceDataSet does not contain all neccesary genes')

    }

    # no NA values for pathway genes (only need to check median)
    medians <- unname(apply(ReferenceDataSet[genes,], 2, median))
    if (sum(is.na(medians)) > 0) {

        stop('ReferenceDataSet has NA values for pathway genes')

    }
    
    # expression data is log transformed
    if (max(ReferenceDataSet[!is.na(ReferenceDataSet)]) > 50) {

        warning(paste('ReferenceDataSet has a maximum value of',
                      max(ReferenceDataSet),
                      'check that the data is log transformed'))

    }

    # data is positive  
    if (min(ReferenceDataSet[!is.na(ReferenceDataSet)]) < 0) {
    
        stop('ReferenceDataset should be strictly non-negative.')

    }
  
}


