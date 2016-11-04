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
#' @export
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
    gsMatrix <- matrix(0,length(times),length(pathway[["max"]]) * sampSize)
    colnames(gsMatrix) <- rep(pathway[["genes"]], sampSize)
    rownames(gsMatrix) <- times

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

        # get gene expresssion
        gsMatrix[i,] <- getExpression(scalingFactor, pathway, singleCell)

    }

    # return the gene expression matrix
    return (gsMatrix)

}

# get scaling factor related to a cellular function
getScalingFactor <- function(model, cells, t, timeWindow, type) {

    # switch on type of cellular function
    switch (type,

        # pathway is related to the G to S transition in the cell cycle
        S =  { return (getGtoSexpression(model, cells, t, timeWindow)) },

        # pathway is related to the G to M transition in the cell cycle
        M = { return (getGtoMexpression(model, cells, t, timeWindow)) },

        # pathway is related to the growth rate of a cell
        GROWTH = { return (getGROWTHexpression(model, cells, t)) },

        # pathway is related to number of neighboring cells
        PROX = { return (getPROXexpression(model, cells, t)) }

    )

}

# get gene expression, given a pathway and a scaling factor
getExpression <- function(scalingFactor, pathway, singleCell) {

    # get range of expression values
    range <- pathway[["max"]] - pathway[["min"]]

    # multiply pathway by scaling factor
    if (singleCell) {

        # multiply each cells value by each gene in the pathway
        return (c(t(scalingFactor %*% t(range) + pathway[["min"]])))

    } else {

        # average the gene expression across cells
        return (mean(scalingFactor) * range + pathway[["min"]])

    }

}

# calculate gene expression for G to S related pathway
getGtoSexpression <- function(model, cells, time, time_window) {

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
getGtoMexpression <- function(model, cells, time, time_window) {

    # get the axis lengths of each cell at current time 
    cur_ax <- getAxisLength(model, time)[cells]

    # get the axis lengths of each cell at end of time window
    next_ax <- getAxisLength(model, time + time_window)[cells]

    # return logical vector (0's and 1's), 1 if cell made the
    # G to M transition in this window
    return (next_ax < cur_ax)

}

# calculate gene expression for a pathway effected by neighboring cells
getPROXexpression <- function(model, cells, time) {

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
getGROWTHexpression <- function(model, cells, time) {

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

        warning(paste('The following pathways are invalid and will not be simulated: ',
                        names(pathways)[!(names(pathways) %in% names(validPwys))]))

    }

    return (validPwys)

}

#' \code{getPathwayExpressionValues} Determine the gene expression values to use for each gene in each pathway.
#' @param pathway A list of genes associated with pathways.
#' @param ReferenceDataSet Optional; Reference gene expression dataset to use to calculate the gene expression values for genes in the pathway to use as mean values in the simulation. Defaults values randomly selected from an exponential distribution with parameter \code{lambda}. If specified by the user, \code{row.names} of the ReferenceDataset must match gene names in the \code{pathway} argument. In this case, gene expression values will be set as the values for genes in the pathway of the sample with max median value of raw RNA expression and expression Z-score for pathway genes. 
#' @param lambda Optional; Parameter of the exponential distribution used to determine the maximum expression value of each simulated gene in the pathway. Defaults to 1/3. Not used if values are determined from a dataset in \code{ReferenceDataSet}.
#' @return list of numeric objects for the maximum expression value of each gene in each pathway. 
#' 

# Determine the range of gene expression values to use for each pathway
getPathwayExpressionRange <- function(ReferenceDataset = NULL, lambda = 1/3, pathways) {

    # if no reference data set provided, sample values from exponential distribution
    if (is.null(ReferenceDataset)) {
    
        # iterate through each pathway
        for (pwy in names(pathways)) {

            # get number of genes in pathway
            pwy_len <- length(pathways[[pwy]][["genes"]])

            # sample min and max (min + difference) values from exponential distribution
            pathways[[pwy]][["min"]] <- rexp(pwy_len, 1) + 3
            pathways[[pwy]][["max"]] <- pathways[[pwy]][["min"]] + rexp(pwy_len, 1.9) + 0.75

        }

    } else {
    
        D <- sweep(ReferenceDataset,1,apply(ReferenceDataset,1,max),FUN="/")
    
    # limit pathways to genes contained in the reference dataset
    pathways <- lapply(pathways, intersect, row.names(ReferenceDataset))
    
    # check that the reference dataset is a valid, log transformed dataset
    checkReferenceDataset(ReferenceDataset,pathways)
    
    if (any(sapply(pathways,length)==0)) {
      stop(paste('The following pathways do not have any genes in the reference dataset: ',
                 paste(names(which(sapply(pathways,length)==0)), collapse=",")))
    }
    
    # for each pathway return the gene expression values in the reference dataset for the sample
    # that has the maximum median gene expression for all pathway genes
    return(lapply(pathways,function(x){ReferenceDataset[x,names(which.max(apply(D[x,],2,median)))]}))

  }
}


