#' \code{plotInteractive} Plots a CellModel at a certain point in time
#'
#' @param model A CellModel
#' @param time The timestep at which to plot the matrix. Must be below
#'      the specified max amount of timesteps
#' @return Plot a visual representation of cells that 
#' takes in command-line-like inputs as well
#' see "plotCellsAtTime and plotInteractive" in vignette for
#' command prompt inputs.
#' @export

setGeneric("plotInteractive", function(model,time = 1)
  standardGeneric("plotInteractive"))

#' \code{plotInteractive} Plots a CellModel at a certain point in time
#'
#' @param model A CellModel
#' @param time The timestep at which to plot the matrix. Must be below
#'      the specified max amount of timesteps
#' @return Plot a visual representation of cells that 
#' takes in command-line-like inputs as well
#' see "plotCellsAtTime and plotInteractive" in vignette for
#' command prompt inputs.
#' @export


setMethod("plotInteractive", "CellModel",
          
    function(model, time = 1) {

    	while (time <= nrow(model)) {
          
			plotCellsAtTime(model,time)

			read = readline()
			place = unlist(gregexpr(" ",read))[1]

			if (place == -1) {

				place = nchar(read)

			}

			cmd = gsub(" ","",substring(read,1,place))
			arg_num = suppressWarnings(as.numeric(gsub(" ","",substring(read,place,nchar(read)))))
			cmds <- c("n","b","t","s","q","h")

			if (is.na(arg_num)) {

				arg_num <- 1

			}

			if ((cmd %in% cmds)) {
		
				switch (match(cmd,cmds),

					{time = time + arg_num},
					{time = time - arg_num},
					{time = arg_num},
					{
						cat("Cell Density = ", getDensity(model,time), "\n")
						cat("Number of Cells = ", length(getGrowthRateDistribution(model,time)), "\n")
											
					},
					{
						break
					},
					{
						
						cat("Basic Commands: \n b = back one timestep \n n = forward one timestep \n s = summary of cells \n q = quit \"console\"\n h = basic command help\n")
					}

				)

			} else {

				cat("Enter a valid command. Type \"h\" for further help\n")

			}
		}
	}
 
)
