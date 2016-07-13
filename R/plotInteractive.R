#' \code{plotInteractive} Plots a CellModel at a certain point in time
#'
#' @param mat A CellModel object
#' @param time The timestep at which to plot the matrix. Must be below
#'      the specified max amount of timesteps
#' @return Plot a visual representation of cells that 
#' takes in command-line-like inputs as well
#' see "plotCellsAtTime and plotInteractive" in vignette for
#' command prompt inputs.
#' @export

setGeneric("plotInteractive", function(mat,time)
  standardGeneric("plotInteractive"))

setMethod("plotInteractive", "CellMatrix",
          
          function(mat,time){
            while(time <= nrow(mat)) {
              radii = seq(3,ncol(mat),6)
              grrate = seq(6,ncol(mat),6)
              plotCellsAtTime(mat,time)
              read = readline()
              count = 1
              place = unlist(gregexpr(" ",read))[1]
              if(place == -1){
                place = nchar(read)
              }
              cmd = gsub(" ","",substring(read,1,place))
              argunum = suppressWarnings(as.numeric(gsub(" ","",substring(read,place,nchar(read)))))
              arguchr = gsub(" ","",substring(read,place,nchar(read)))
              testcase = suppressWarnings(as.numeric("abcdefghijklmnopqrstuvwxyz"))
              cmds <- c("n","b","t","summ","q","h")
              if((cmd %in% cmds)==TRUE){
                if(suppressWarnings(identical(argunum,testcase) == FALSE)){
                  count = argunum
                }
                if(match(cmd,cmds) == 1){
                  time = time + count
                }
                else if(match(cmd,cmds) == 2){
                  time = time - count
                }
                else if(match(cmd,cmds) == 3){
                  time = count
                }
                else if(match(cmd,cmds) == 4){
                  vari = sum(mat[time,radii]>0)
                  if(arguchr == "d"){
                    vari = getDensity(mat,time)
                    cat(paste("Density of Cells at time",time,"=",vari))
                  }
                  else if(arguchr == "g"){
                    vari = (sum(mat[time,grrate]))/(sum(mat[time,radii]))
                    cat(paste("Average Growth Rate of Cells at time",time,"=",vari))
                  }
                  else{
                    cat(paste("Number of Cells at time",time,"=",vari))
                  }
                }
                else if(match(cmd,cmds) == 5){
                  graphics.off()
                  break
                }
                else if(match(cmd,cmds) == 6){
                  cat("Basic Commands: \n b = back one timestep \n n = forward one timestep \n summ = summary of cells \n q = quit \"console\"\n h = basic command help")
                }
              }
              else{
                cat("Enter a valid command. Type \"h\" for further help")
              }
              if(time > nrow(data) | time <= 0){
                print("Time out of bound")
                graphics.off()
                break
              }
            }
          }
)