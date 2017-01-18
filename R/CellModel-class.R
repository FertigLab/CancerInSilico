#### class definition ####

#' @title CellModel
#' @description An S4 class that contains the output of a cell simulation,
#'      along with all the parameters used to generate the simulation
#' @slot mCells A list object where each row of the list describes the
#'      state of all the cells in the model at a given time. Each cell is
#'      described over 6 columns: [1] x-coordinate, [2] y-coordinate, [3]
#'      radius, [4] axis length, [5] axis angle, [6] growth rate. For
#'      instance, the x-coordinates of the first 3 cells will be in columns
#'      1,7,13.
#' @slot mParameters A list object containing all parameters used to run
#'      the model
#' @export

setClass("CellModel", representation(cells = "list", params = "list"))

#### constructor ####

createCellModel <- function(params, output) {

    return (new("CellModel", cells = output, params = params))

}

#### getters ####

# helper function to find corresponding row of a given time
timeToRow <- function(model, time) {

    if (time == model@params[['runTime']]) {

        return (length(model@cells))

    } else {

        return (floor(time / model@params[['recordIncrement']] + 1))

    }

}


# get column of data for each cell
getColumn <- function(model, time, col) {

    # find row corresponding to time
    row <- timeToRow(model, time)

    # get the sequence of indices that correspond to column
    indices <- seq(col,length(model@cells[[row]]), 6)

    # return the values at these indices
    return(model@cells[[row]][indices])

}    

# return an N X 2 matrix of cell coordinates at time
getCoordinates <- function(model, time) {

    # get x and y coordinates
    xCoord <- getColumn(model, time, 1)
    yCoord <- getColumn(model, time, 2)

    # return matrix with coordinates
    return (matrix(c(xCoord, yCoord), nrow = length(xCoord)))

}

# get a vector containing the radius of each cell at time
getRadii <- function(model, time) {

    return (getColumn(model, time, 3))

}

# get a vector containing the radius of each cell at time
getAxisLength <- function(model, time) {

    return (getColumn(model, time, 4))

}

# get a vector containing the radius of each cell at time
getAxisAngle <- function(model, time) {

    return (getColumn(model, time, 5))

}

# get a vector containing the radius of each cell at time
getGrowthRates <- function(model, time) {

    return (getColumn(model, time, 6))

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
    grRates <- getGrowthRates(model, time)

    # convert the raw growth rates to cell cycle time in hours
    return (1 + 2 * (sqrt(2) - 1) * model@params[["timeIncrement"]] *
        model@params[["nG"]] / grRates)

}

#' \code{getNumberOfCells} get the number of cells alive
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return the number of cells at this time
#' @examples 
#' getNumberOfCells(runCancerSim(1,1), 1)
#' @export

getNumberOfCells <- function(model, time) {

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

    # get radius of each cell
    radii <- getRadii(model, time)

    # not sure about how R treats non-logical values, so this is safe
    if (model@params[['boundary']] != FALSE) {

        # ratio of area of all cells and the disk which contains them       
        return(sum(radii ** 2) / (model@params[['boundary']] ^ 2))
        
    } else {

        # get coordinates of each cell
        coords <- getCoordinates(model, time)

        # farthest distance from (0,0) of cell
        d <- max(sqrt(coords[,1]**2 + coords[,2]**2) + radii)

        # ratio of area of all cells and the disk which contains them       
        return(sum(radii ** 2) / (d ^ 2))

    }

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
    while (time <= model@params[['runTime']] && time >= 0) {
      
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
                    cat(paste("s = summary of cells\nq = quit\n",
                                "h = basic command help\n"))
                }

            )

        # invalid command                
        } else {

            cat("Enter a valid command. Type \"h\" for further help\n")

        }

    }

}

#' \code{getCellTypes} return the list of cell types
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return the list of cell types for the model

getCellTypes <- function(model, time) {

    return (model@params[['cellTypeInitFreq']])

}

#### helper (non-exported) functions ####

#' \code{getNumNeighbors} return the number of neighboring cells
#'
#' @param model A CellModel
#' @param time time in model hours
#' @param index index of cell
#' @return the number of neighbors around the cell at index

getNumNeighbors <- function(model, time, index) {

    # initialize number of neighbors
    num <- 0 

    # get coordinates of all cells
    coords <- getCoordinates(model, time)
    
    # search through all cells besides the one at index
    for (i in setdiff(1:nrow(coords), index)) {

        # if radius is within distance of inner ring of neighbors
        if ((coords[i,1] - coords[index,1])^2 + (coords[i,2] -
            coords[index,2])^2 < 12.25) {

            # add to neighbor count
            num <- num + 1

        }

    }

    # return neighbor count
    return (num)

}

#' \code{plotCell} plots a CellModel at a given time
#'
#' @param model A CellModel
#' @param time time in model hours
#' @return plot a visual representation of cells 
#' @examples plotCells(runCancerSim(10,1), 1)
#' @export

plotCells <- function(model,time,drawBoundary = TRUE)  {

    # get all the cell information
    coords <- getCoordinates(model, time)
    radii <- getRadii(model, time)
    axis_len <- getAxisLength(model, time)
    axis_ang <- getAxisAngle(model, time)

    if (model@params[['boundary']]) {

        mn <- -model@params[['boundary']] - 2
        mx <- model@params[['boundary']] + 2

    } else {

        mn <- min(coords) - 2
        mx <- max(coords) + 2

    }

     
    # create the plot template
    plot(c(mn, mx), c(mn, mx), main = paste("Plot of CellModel At Time",
        time), xlab = "", ylab = "", type = "n", asp = 1)
          
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
    symbols(x,y, circles=rad, inches=FALSE, add=TRUE, bg="bisque4",
        fg="bisque4")

    # draw boundary
    if (drawBoundary) {

        symbols(0,0,circles= model@params[['boundary']], inches=FALSE, add=TRUE,
                    lwd = 2)

    }

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
        print(object@params)
        cat("available functions:\n")
        cat("interactivePlot\n")
        cat("plotCells\n")
        cat("getDensity\n")
        cat("getCycleLengths\n")
        cat("getNumberOfCells\n")

    }

)
