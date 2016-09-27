#' \code{interactivePlot} Plots a CellModel and allows the user to control the plot with various commands
#'
#' @param model A CellModel
#' @param time The timestep at which to plot the matrix. Must be below
#'      the specified max amount of timesteps
#' @return Plot a visual representation of cells that 
#' takes in command-line-like inputs as well
#' see "plotCellsAtTime and plotInteractive" in vignette for
#' command prompt inputs.
#' @export
#' @rdname CellModel-class

interactivePlot <- function(model, time = 0) {

    defnum = 1
    while (time <= .runTime(model)) {
      
        plotCellsAtTime(model,time)

        read = readline()
        place = unlist(gregexpr(" ",read))[1]

        if (place == -1) {

	        place = nchar(read)

        }

        cmd = gsub(" ","",substring(read,1,place))
        arg_num = suppressWarnings(as.numeric(gsub(" ","",substring(read,place,nchar(read)))))
        
        cmds <- c("n","b","t","i","s","q","h")

        if (is.na(arg_num)) {

            arg_num = defnum

        }

		if ((cmd %in% cmds)) {
	
			switch (match(cmd,cmds),
        
                {time = time + arg_num},
                {time = time - arg_num},
                {time = arg_num},
		        {defnum = arg_num},
                {
                    cat("Cell Density = ", getDensity(model,time), "\n")
					cat("Number of Cells = ", length(getCycleLengthDistribution(model,time)), "\n")
                },
                {
                    graphics.off()
                    break
                },
                {
                    cat("Basic Commands:\n")
                    cat("b ARG = back ARG timesteps (default ARG = 1)\n")
                    cat("n ARG = forward ARG timesteps (default ARG = 1)\n")
                    cat("t ARG - jump to timestep ARG (default ARG = 1)\n")
                    cat("i ARG - change default ARG for other commands\n")
                    cat("s = summary of cells\nq = quit\nh = basic command help\n")
                }

			)

        } else {

            cat("Enter a valid command. Type \"h\" for further help\n")

        }

    }

}

#'\code{getDrugEffect} 
#'
#' @param cycleLengthSeq A sequence spanning the range of cycle lengths
#' @return A list of vectors specifying the distribution of drug effects depending on the growth rate of the cell
#' @examples
#' getDrugEffect(0.3, seq(2,12,0.1))
#' @export

getDrugEffect <- function(FUN = function(x) {0}, ...) {

    cycleLengthDist <- list(...)$cycleLengthDist
    cycleLengthSeq <- list(...)$cycleLengthSeq

    if (is.null(cycleLengthDist)) {

        if (is.null(cycleLengthSeq)) {

            stop("must specifiy either cycleLengthDist or cycleLengthSeq")

        }

    } else if (is.null(cycleLengthSeq)) {

        cycleLengthSeq <- seq(min(cycleLengthDist), max(cycleLengthDist), 0.1)

    } 

    ret_list <- vector("list", length(cycleLengthSeq))

    for (i in 1:length(ret_list)) {

        ret_list[[i]] = c(cycleLengthSeq[i], FUN(cycleLengthSeq[i]))

    }

    return (ret_list)

}
