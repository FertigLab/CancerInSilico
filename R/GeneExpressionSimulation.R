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

    time_window = 1 # arbitrary for now
    mean_cycle_len = 16 # avg cycle length in hours
    times <- seq(0, .runTime(model) - time_window, sampFreq)

    gsMatrix <- matrix(0,length(times),length(pathway) * sampSize)
    colnames(gsMatrix) <- rep(names(pathway), sampSize)
    rownames(gsMatrix) <- times

    for (t in times) {

        cells <- 1 + (seq(1,getNumberOfCells(model, t), 6) - 1) / 6
        if (singleCell) {

            cells <- sample(cells, sampSize)

        }

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

                cycle_len <- getCycleLengths(model,t)[cells]
                single_cell_exp <- 1 - 1 / (1 + exp(-0.25 * (cycle_len - mean_cycle_len)))
    
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

}


