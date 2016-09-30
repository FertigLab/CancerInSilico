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

            S =  {

                cur_rad <- getRadii(model, t)[cells]
                next_rad <- getRadii(model, t + time_window)[cells]
                single_cell_exp <- next_rad > sqrt(3/2) & cur_rad < sqrt(3/2)

            }, M = {

                cur_ax <- getAxisLength(model, t)[cells]
                next_ax <- getAxisLength(model, t + time_window)[cells]
                single_cell_exp <- next_ax < cur_ax

            }, GROWTH = {

                cur_rad <- getRadii(model, t)[cells]
                next_rad <- getRadii(model, t + time_window)[cells]
                not_growing <- cur_rad > next_rad
                single_cell_exp <- (next_rad - cur_rad) / cur_rad
                single_cell_exp[not_growing] = 0

            }, PROX = {

                single_cell_exp <- c()
                for (i in cells) {

                    single_cell_exp[i] <- getNumNeighbors(model, t, i, 2 * sqrt(3) - 0.01) / 6

                }

            }

        )

        if (singleCell) {

            gsMatrix[t,] = c(t(pathway %*% t(single_cell_exp)))

        } else {

            gsMatrix[t,] = pathway * mean(single_cell_exp)

        }

    }

    return (gsMatrix)

}


