#' @title CellModel
#' @description An S4 class that contains the output of a cell simulation,
#'  along with all the parameters used to generate the simulation
#' @slot cells A list object that describes the state of the cells
#'  at each time. The state representation depends on the type of
#'  model run, and is accessed by the function designed for each
#'  model type. 
#' @slot params A list object containing all parameters used in the model
#' @export
setClass('DrasdoHohmeModel', contains = 'OffLatticeModel', representation(
    nG = 'numeric',
    epsilon = 'numeric',
    delta = 'numeric')
)


setMethod('processParameters',
    signature('DrasdoHohmeModel'),
    function(model, ...)
    {
        # for readability
        p <- model@params

        # get extra parameters
        if (is.null(p[['nG']])) {p[['nG']] <- list(...)$nG}
        if (is.null(p[['epsilon']])) {p[['epsilon']] <- list(...)$epsilon}
        if (is.null(p[['delta']])) {p[['delta']] <- list(...)$delta}

        # check drasdo parameters
        if (is.null(p[['nG']])) stop('missing nG')
        if (is.null(p[['epsilon']])) stop('missing epsilon')
        if (is.null(p[['delta']])) stop('missing delta')

        if (p[['nG']] <= 0) stop('invalid nG')
        if (p[['epsilon']] <= 0) stop('invalid epsilon')
        if (p[['delta']] <= 0) stop('invalid delta')

        # call function on top level base class
        base <- new('CellModel', params = p)
        p <- processParameters(base)

        # get minimum cycle length
        minCycle <- p[['cellTypes']][[1]]@minCycle
        for (i in 1:length(p[['cellTypes']]))
        {
            minCycle <- min(minCycle, p[['cellTypes']][[i]]@minCycle)
        }

        # get off lattice parameters
        e <- p[['epsilon']]
        d <- p[['delta']]
        ng <- p[['nG']]

        t1 <- d / (4 * ng * (4 - sqrt(2)))
        t2 <- d * (minCycle - 1) / (8 * ng * (sqrt(2) - 1))
        p[['timeIncrement']] <- min(t1, t2)

        p[['maxTranslation']] <- d / 2
        p[['maxRotate']] <- acos((16 + d^2 - 4 * d) / 16)
        p[['maxDeform']] <- 2 * p[['timeIncrement']] * ng * (4 - sqrt(2))

        # call function on direct base class
        base <- new('OffLatticeModel', params = p)
        p <- processParameters(base)

        # return parameters
        return (p)        
    }
)
