#### class definition ####

#' @title CellModel
#' @description An S4 class to represent the output of a cell-based model
#'
#' @slot mCells A list object where each row of the list describes the state of all the cells in the model at a given time. Each cell is described over 6 columns: [1] x-coordinate, [2] y-coordinate, [3] radius, [4] axis length, [5] axis angle, [6] growth rate. For instance, the x-coordinates of the first 3 cells will be in columns 1,7,13.
#' @slot mInitialNumCells the initial number of cells in the model
#' @slot mRunTime the total run time (hours) of the model
#' @slot mInitialDensity the density the cells were seeded at
#' @slot mInheritGrowth whether or not cells inherit growth rates from their parent
#' @slot mOutputIncrement the frequency of print statements during the run
#' @slot mRandSeed the random seed 
#' @slot mEpsilon model specific parameter 
#' @slot mNG model specific parameter
#' @slot mTimeIncrement amount of time elapsed in each model step
#' @slot mCycleLengthDist initial distribution of cell-cycle lengths 
#' @export

setClass("CellModel", representation(
                        mCells = "list",
                        mInitialNumCells = "numeric",
                        mRunTime = "numeric",
                        mInitialDensity = "numeric",
                        mInheritGrowth = "logical",
                        mOutputIncrement = "numeric",
                        mRandSeed = "numeric",
                        mEpsilon = "numeric",
                        mNG = "numeric",
                        mTimeIncrement = "numeric",
                        mCycleLengthDist = "numeric" ))

#### getters (parameters) ####

.initialNumCells <- function(model) { return (model@mInitialNumCells) }

.runTime <- function(model) { return (model@mRunTime) }

.initialDensity <- function(model) { return (model@mInitialDensity) }

.inheritGrowth <- function(model) { return (model@mInheritGrowth) }

.outputIncrement <- function(model) { return (model@mOutputIncrement) }

.randSeed <- function(model) { return (model@mRandSeed) }

.epsilon <- function(model) { return (model@mEpsilon) }

.nG <- function(model) { return (model@mNG) }

.timeIncrement <- function(model) { return (model@mTimeIncrement) }

.cycleLengthDist <- function(model) { return (model@mCycleLengthDist) }

#### getters (cell data) ####

#' \code{getCoordinates} get a two dimensional matrix of all the cell coordinates
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return an N X 2 matrix of cell coordinates at time

getCoordinates <- function(model, time) {

    # find the row corresponding to the given time 
    row <- timeToRow(model, time)

    # get the sequence of indices that contain x-coordinates
    indices <- seq(1,length(model@mCells[[row]]),6)
    
    # create return matrix, col 1 = x-coord & col 2 = y-coord
    ret_mat <- matrix(nrow = length(indices), ncol = 2)

    # for each cell (each index), store the x and y coordinate
    for (i in 1:length(indices)) {

        # store x-coordinate
        ret_mat[i,1] = model@mCells[[row]][indices[i]]

        # store y-coordinate (x-coord index + 1)
        ret_mat[i,2] = model@mCells[[row]][indices[i]+1]

    }

    # return the coordinate matrix
    return (ret_mat)
         
}

#' \code{getRadii} get the radius of each cell
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return vector containing the radius of each cell at time

getRadii <- function(model, time) {

    # find the row corresponding to the given time
    row <- timeToRow(model, time)

    # get the sequence of indices that contain the cell radius (starts at 3)
    indices <- seq(3,length(model@mCells[[row]]),6)

    # return the values at these indices
    return(model@mCells[[row]][indices])

}

#' \code{getAxisLength} get the axis length of each cell
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return vector containing the axis length of each cell at time
         
getAxisLength <- function(model, time) {

    # find the row corresponding to the given time
    row <- timeToRow(model, time)

    # get the sequence of indices that contain the axis length (starts at 4)
    indices <- seq(4,length(model@mCells[[row]]),6)

    # return the values at these indices
    return(model@mCells[[row]][indices])

}

#' \code{getAxisAngle} get the axis angle of each cell
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return vector containing the axis angle of each cell at time

getAxisAngle <- function(model, time) {

    # find the row corresponding to the given time
    row <- timeToRow(model, time)

    # get the sequence of indices that contain the axis angle (starts at 5)
    indices <- seq(5,length(model@mCells[[row]]),6)

    # return the values at these indices
    return(model@mCells[[row]][indices])

}

#' \code{getGrowthRates} get the model growth rates of each cell
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return vector containing the growth rate of each cell at time

getGrowthRates <- function(model, time) {

    # find the row corresponding to the given time
    row <- timeToRow(model, time)

    # get the sequence of indices that contain the cell growth rate (starts at 6)
    indices <- seq(6,length(model@mCells[[row]]),6)

    # return the values at these indices
    return(model@mCells[[row]][indices])

}

#### exported functions for this class ####

#' \code{getCycleLengths} return the cycle lengths of each cells at time
#'
#' @param model a CellModel object
#' @param time time in model hours
#' @return the cycle lengths of each cell at time
#' @examples
#' getCycleLengths(runCancerSim(1,1), 1)
#' @export

getCycleLengths <- function(model, time) {

    # get the raw model growth rates
    gr_rates <- getGrowthRates(model, time)

    # convert the raw growth rates to cell cycle time in hours
    return (1 + 2 * (sqrt(2) - 1) * .timeIncrement(model) * .nG(model) / gr_rates)

}

#' \code{getNumberOfCells} get the number of cells alive
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return the number of cells at this time
#' @examples 
#' getNumberOfCells(runCancerSim(1,1), 1)
#' @export

getNumberOfCells <-    function(model, time) {

    return (sum(getRadii(model, time) > 0))   

}

#' \code{getDensity} gets the density of cells at a given time
#'
#' @param model A Cell Model
#' @param time time in model hours
#' @return The density of cells at that time (not quite the same as confluency)
#' @examples
#' getDensity(runCancerSim(1,1),1)
#' @export

getDensity <- function(model,time) {

    # get coordinates of each cell
    coords <- getCoordinates(model, time)

    #get radius of each cell
    radii <- getRadii(model, time)
    
    # farthest distance from (0,0) of cell
    d <- max(sqrt(coords[,1]**2 + coords[,2]**2) + radii)

    # ratio of area of all cells and the disk which contains them       
    return(sum(radii ** 2) / (d ^ 2))

}


#' \code{getParameters} get a named list of parameters in the model
#'
#' @param model A CellModel
#' @param fullDist [bool] return full distribution of cycle length
#' @return a named list of parameters in the model
#' @examples
#' getParameters(runCancerSim(1,1))
#' @export

getParameters <- function(model, fullDist=FALSE) {

    # get the average cell cycle time
    retDist = mean(.cycleLengthDist(model))

    # if the full distribution is requested, store the entire distribution
    if (fullDist) {

        retDist = .cycleLengthDist(model)

    }

    # create a named list of all the parameters
    ret_val = list(

        initialNum = .initialNumCells(model),
        runTime = .runTime(model),           
        initialDensity = .initialDensity(model),           
        inheritGrowth = .inheritGrowth(model),           
        outputIncrement = .outputIncrement(model),           
        randSeed = .randSeed(model),           
        epsilon = .epsilon(model),           
        nG = .nG(model),
        timeIncrement = .timeIncrement(model),
        cycleLengthDist = retDist

    )           

    # return named list
    return(ret_val)

}

#' \code{interactivePlot} plots a CellModel and allows the user to control the plot with various commands
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return plot a visual representation of cells that takes in command-line-like inputs, type 'h' for help and a list of all available commands
#' @export

interactivePlot <- function(model, time = 0) {

    # default arg value    
    default_arg = 1

    # while time is valid
    while (time <= .runTime(model) && time >= 0) {
      
        # call internal function which plots cells at current time
        plotCells(model,time)

        # get keyboard input
        read = readline()

        # find separator between command and arg
        place = unlist(gregexpr(" ",read))[1]

        if (place == -1) {

            place = nchar(read)

        }

        # parse command
        cmd = gsub(" ","",substring(read,1,place))
        arg = suppressWarnings(as.numeric(gsub(" ","",substring(read,place,nchar(read)))))
        
        # list of possible commands
        cmds <- c("n","b","t","i","s","q","h")

        # set to default if no command is present
        if (is.na(arg)) {

            arg = default_arg

        }

        # check if valid command
        if ((cmd %in% cmds)) {
    
            # switch on which command was provided
            switch (match(cmd,cmds),
        
                # 'n' increase time by arg
                {time = time + arg},

                # 'b' decrease time by arg
                {time = time - arg},

                # 't' go to time arg
                {time = arg},

                # 'i' set default arg
                {default_arg = arg},

                # 's' get cell summary
                {
                    cat("Cell Density = ", getDensity(model,time), "\n")
                    cat("Number of Cells = ", getNumberOfCells(model, time), "\n")
                },

                # 'q' quit display
                {
                    graphics.off()
                    break
                },

                # display help command
                {
                    cat("Basic Commands:\n")
                    cat("b ARG = back ARG timesteps (default ARG = 1)\n")
                    cat("n ARG = forward ARG timesteps (default ARG = 1)\n")
                    cat("t ARG - jump to timestep ARG (default ARG = 1)\n")
                    cat("i ARG - change default ARG for other commands\n")
                    cat("s = summary of cells\nq = quit\nh = basic command help\n")
                }

            )

        # invalid command                
        } else {

            cat("Enter a valid command. Type \"h\" for further help\n")

        }

    }

}

#### helper (non-exported) functions ####

#' \code{timeToRow} return the correct row in the mCells list corresponding to a given time
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return corresponding row in mCells list

timeToRow <- function(model, time) {

    return (floor(time / .timeIncrement(model)) + 1)

}

#' \code{plotCell} plots a CellModel at a given time
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return plot a visual representation of cells 

plotCells <- function(model,time)  {

    # get all the cell information
    coords <- getCoordinates(model, time)
    radii <- getRadii(model, time)
    axis_len <- getAxisLength(model, time)
    axis_ang <- getAxisAngle(model, time)

    # find a square that contains all cells
    mn <- min(coords) - 2
    mx <- max(coords) + 2
    
    # create the plot template
    plot(c(mn,mx),c(mn,mx),main=paste("Plot of CellModel At Time",time),xlab = "",ylab="",type="n",asp=1)
          
    # get all (x,y) pairs for each of the two circles that make up a cell
    x_1 <- coords[,1] + (0.5 * axis_len - radii) * cos(axis_ang)
    x_2 <- coords[,1] - (0.5 * axis_len - radii) * cos(axis_ang)
    y_1 <- coords[,2] + (0.5 * axis_len - radii) * sin(axis_ang)
    y_2 <- coords[,2] - (0.5 * axis_len - radii) * sin(axis_ang)

    # combine all coordinate pairs along with the radii for each cell    
    x <- c(x_1,x_2)
    y <- c(y_1,y_2)
    rad <- c(radii, radii)
    
    # plot the cells
    symbols(x,y, circles=rad, inches=FALSE, add=TRUE, bg="bisque4", fg="bisque4")

}

#' \code{show} display summary of CellModel class
#'
#' @param object A CellModel Object
#' @examples show(runCancerSim(1,1))
#' @return shows all available functions and parameters of model
#' @export

setMethod("show", "CellModel",

    function (object) {

        cat("model parameters:")
        print(getParameters(object))
        cat("available functions:\n")
        cat("interactivePlot\n")
        cat("getParameters\n")
        cat("getDensity\n")
        cat("getCycleLengths\n")
        cat("getNumberOfCells\n")

    }

)

