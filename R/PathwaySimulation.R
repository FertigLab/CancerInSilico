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
    sampSize <- 1

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




