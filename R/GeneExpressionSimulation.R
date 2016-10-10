#' \code{simulatePathway} Simulates gene expression data for a given pathway
#'
#' @param model A CellModel
#' @param pathway The pathway to simulate expression for
#' @param type the type of cellular function the pathway is related too
#' @param sampFreq how often to sample the cells and record gene expression data
#' @param sampSize how many cells to sample in the case of single-cell expression data
#' @param singleCell whether or not to simulate single-cell expression data
#' @return gene expression matrix for given pathway 
#' @export
#'

simulatePathway <- function(model, pathway, type, sampFreq = 1, sampSize = 1, singleCell=FALSE) {

    # time window to search for events related to gene expression
    time_window = 1 

    # average cycle time used for determining if cell is growing 'fast' or 'slow'
    mean_cycle_len = 16

    # vector of times to get gene expression for
    times <- seq(0, .runTime(model) - time_window, sampFreq)

    # create return matrix
    gsMatrix <- matrix(0,length(times),length(pathway) * sampSize)
    colnames(gsMatrix) <- rep(names(pathway), sampSize)
    rownames(gsMatrix) <- times

    # loop through each time
    for (t in times) {

        # get a vector of all cells
        cells <- 1:getNumberOfCells(model, t)

        # if doing single cell ...
        if (singleCell) {
    
            # get a sample of cells
            cells <- sample(cells, sampSize)

        }

        # switch on type of gene expression
        switch (type,

            # pathway is related to the G to S transition in the cell cycle
            S =  { single_cell_exp <- getGtoSexpression(model, cells, t, time_window) },

            # pathway is related to the G to M transition in the cell cycle
            M = { single_cell_exp <- getGtoMexpression(model, cells, t, time_window) },

            # pathway is related to the growth rate of a cell
            GROWTH = { single_cell_exp <- getGROWTHexpression(model, cells, t) },

            # pathway is related to number of neighboring cells
            PROX = { single_cell_exp <- getPROXexpression(model, cells, t) }

        )

        # if single cell expression is desired
        if (singleCell) {

            # multiply each cells value by each gene in the pathway
            gsMatrix[t,] = c(t(pathway %*% t(single_cell_exp)))

        } else {

            # average the gene expression across cells
            gsMatrix[t,] = pathway * mean(single_cell_exp)

        }

    }

    # return the gene expression matrix
    return (gsMatrix)

}

#' \code{getGtoSexpression} calculate gene expression for a pathway effected by the G to S transition
#'
#' @param model A CellModel
#' @param cells the indices of cells to calculate expression for
#' @param time model time in hours
#' @param time_window the window of time (model hours) to check if a cell made the transition
#' @return gene expression for each cell 

getGtoSexpression(model, cells, time, time_window) {

    # get the axis lengths of each cell at current time 
    cur_rad <- getRadii(model, time)[cells]

    # get the axis lengths of each cell at end of time window
    next_rad <- getRadii(model, time + time_window)[cells]

    # return logical vector (0's and 1's), 1 if cell made the G to S transition in this window
    # G to S transition is defined by a cell crossing the halfway point of its growth
    return (next_rad > sqrt(3/2) & cur_rad < sqrt(3/2))

}

#' \code{getGtoMexpression} calculate gene expression for a pathway effected by the G to M transition
#'
#' @param model A CellModel
#' @param cells the indices of cells to calculate expression for
#' @param time model time in hours
#' @param time_window the window of time (model hours) to check if a cell made the transition
#' @return gene expression for each cell 

getGtoMexpression(model, cells, time, time_window) {

    # get the axis lengths of each cell at current time 
    cur_ax <- getAxisLength(model, time)[cells]

    # get the axis lengths of each cell at end of time window
    next_ax <- getAxisLength(model, time + time_window)[cells]

    # return logical vector (0's and 1's), 1 if cell made the G to M transition in this window
    return (next_ax < cur_ax)

}

#' \code{getPROXexpression} calculate gene expression for a pathway effected by neighboring cells
#'
#' @param model A CellModel
#' @param cells the indices of cells to calculate expression for
#' @param time model time in hours
#' @return gene expression for each cell 

getPROXexpression(model, cells, time) {

    # create empty vector
    exp <- c()

    # loop through each cell
    for (i in cells) {

        # expression = number of neighboring cells (max 6) divided by 6
        exp[i] <- getNumNeighbors(model, time, i, 2 * sqrt(3) - 0.01) / 6

    }

    # return expression vector
    return (exp)

}

#' \code{getGROWTHexpression} calculate gene expression for a pathway effected by the growth rate of the cells
#'
#' @param model A CellModel
#' @param cells the indices of cells to calculate expression for
#' @param time model time in hours
#' @return gene expression for each cell 

getGROWTHexpression(model, cells, time) {

    # get the cycle lengths of each cell
    cycle_len <- getCycleLengths(model,time)[cells]

    # this function maps cycle length to [0,1] in a smooth way
    return (1 - 1 / (1 + exp(-0.25 * (cycle_len - mean_cycle_len))))

}




