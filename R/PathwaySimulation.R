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
simulatePathway <- function(model, pathway, type, sampFreq = 1,
sampSize, singleCell = FALSE, timeWindow = 1, downReg = FALSE) {

    # if not single cell data, sample size is "1" (only return mean)
    if (!singleCell) { sampSize <- 1 }

    # find closest, valid, sampling frequency
    sampFreq <- .recordIncrement(model) *
                    ceiling(sampFreq / .recordIncrement(model))

    # vector of times to get gene expression for
    times <- seq(0, .runTime(model) - timeWindow, sampFreq)

    # create return matrix
    gsMatrix <- matrix(nrow = length(pathway[["min"]]),
                       ncol = length(times) * sampSize)
    colnames(gsMatrix) <- rep(times, each = sampSize)
    rownames(gsMatrix) <- pathway[["genes"]]

    # loop through each time
    for (i in 1:length(times)) {

        # get a vector of all cells
        cells <- 1:getNumberOfCells(model, times[i])

        # if doing single cell, get sample of cells
        if (singleCell) {
    
            cells <- sort(sample(cells, sampSize))

        }

        # get scaling factor for gene expression
        scalingFactor <- getScalingFactor(model, cells, times[i],
            timeWindow, type)

        # rows to replace
        cols <- (sampSize * (i - 1) + 1):(sampSize * i)

        # get gene expresssion
        gsMatrix[,cols] <- getExpression(scalingFactor, pathway, singleCell)

    }

    # return the gene expression matrix
    return (gsMatrix)

}

# get scaling factor related to a cellular function
getScalingFactor <- function(model, cells, t, timeWindow, type) {

    # switch on type of cellular function
    switch (type,

        # pathway is related to the G to S transition in the cell cycle
        S =  { return (scalingFactorGtoS(model, cells, t, timeWindow)) },

        # pathway is related to the G to M transition in the cell cycle
        M = { return (scalingFactorGtoM(model, cells, t, timeWindow)) },

        # pathway is related to the growth rate of a cell
        GROWTH = { return (scalingFactorGROWTH(model, cells, t)) },

        # pathway is related to number of neighboring cells
        PROX = { return (scalingFactorPROX(model, cells, t)) }

    )

}

# get gene expression, given a pathway and a scaling factor
getExpression <- function(scalingFactor, pathway, singleCell) {

    # get range of expression values
    range <- pathway[["max"]] - pathway[["min"]]

    # multiply pathway by scaling factor
    if (singleCell) {

        # multiply each cells value by each gene in the pathway
        return (range %*% t(scalingFactor) + pathway[["min"]])

    } else {

        # average the gene expression across cells
        return (mean(scalingFactor) * range + pathway[["min"]])

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
getPathways <- function(pathways = NULL) {

    # possible pathway types
    validTypes <- c('GtoM', 'GtoS', 'Prox', 'Growth')
  
    # if no user pathways provided
    if (is.null(pathways)) {

        # use default pathways from the package 
        validPwys <- inSilicoPathways

    } else {

        # keep only valid pathways
        validPwys <- pathways[names(pathways) %in% validTypes]
    
    }

    # check that at least 1 pathway is valid
    if (length(validPwys) == 0) {

        stop(paste('No valid pathways provided. Valid choices are: ',
                 paste(validTypes, collapse=', ')))

    # warn user if some of the provided pathways were not valid
    } else if (length(validPwys) < length(pathways)) {

        warning(paste('The following pathways are invalid and will not ',
                'be simulated: ',
                names(pathways)[!(names(pathways) %in%
                                names(validPwys))]))

    }

    return (validPwys)

}

# Determine the range of gene expression values to use for each pathway
getPathwayExpressionRange <- function(ReferenceDataset = NULL,
lambda = 1/3, pathways) {

    # if no reference data set provided, sample values
    if (is.null(ReferenceDataset)) {
    
        # iterate through each pathway
        for (pwy in names(pathways)) {

            # get number of genes in pathway
            pwy_len <- length(pathways[[pwy]][["genes"]])

            # sample min and max (min + diff) values from exponential
            pathways[[pwy]][["min"]] <- rexp(pwy_len, 1) + 3
            pathways[[pwy]][["max"]] <- pathways[[pwy]][["min"]] +
                                        rexp(pwy_len, 1.9) + 0.75

        }

    } else {

        # check that the reference dataset is valid
        checkReferenceDataset(ReferenceDataset)

        # get the names of genes in the data set
        gene_names <- row.names(ReferenceDataset)

        # iterate through each pathway
        for (pwy in names(pathways)) {

            # find pathway genes in data set
            genes <- intersect(gene_names, pathways[[pwy]][["genes"]])

            # get median expression for pathway genes
            D_path <- ReferenceDataset[genes]
            med_exp <- unname(apply(D_path, 2, median))

            # remove genes with NA values
            genes <- genes[!is.na(unname(med_exp[genes]))]

            # if no genes remain, quit with error
            if (length(genes) == 0) {

                stop(paste("The following pathway had genes without any",
                            "values in the reference data set: ",
                            pwy))

            }

            # get min/max expression values by taking the expression
            # values from the sample with min/max median expression
            pathways[[pwy]][["min"]] <- unname(D_path[,which.min(med_exp)])
            pathways[[pwy]][["max"]] <- unname(D_path[,which.max(med_exp)])

        }

    }

}

# verify reference data set
checkReferenceDataset <- function(ReferenceDataset) {

    # make sure expression data is log transformed
    if (max(ReferenceDataset) > 50) {

        warning(paste('ReferenceDataset has a maximum value of',
                      max(ReferenceDataset),
                      'check that the data is log transformed'))

    }

    # make sure data is positive  
    if (min(ReferenceDataset) < 0) {
    
        warning('ReferenceDataset should be strictly non-negative.')

    }
  
}


