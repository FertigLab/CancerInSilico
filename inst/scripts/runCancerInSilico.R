library(CancerInSilico)

args <- commandArgs(TRUE)
numOfRuns <- as.integer(args[1])
arrayNum <- as.integer(args[2])
jobName <- args[3]

#### Set Parameters ####

initialNum <- rep(80, numOfRuns) 
runTime <- rep(168, numOfRuns)
density <- rep(0.2, numOfRuns)
## cycleLengthDist <- rep(48, numOfRuns)
cycleLengthDist <- replicate(numOfRuns, 20 + rexp(1000, 1/5))
drugEffect <- replicate(numOfRuns, function(x) {return(1)})
# drugEffect <- c(function(x) {return(runif(1,0,1))})
inheritGrowth <- c(FALSE, TRUE)
outputIncrement <- rep(6, numOfRuns)
recordIncrement <- rep(0.25, numOfRuns)
randSeed <- rep(0, numOfRuns)
modelType <- rep("DrasdoHohme2003", numOfRuns)
drugTime <- rep(0.0, numOfRuns)
boundary <- rep(TRUE, numOfRuns)
nG <- rep(8,numOfRuns)
epsilon <- rep(10, numOfRuns)
syncCycles <- rep(TRUE, numOfRuns)

########################

output <- runCancerSim (initialNum=initialNum[arrayNum],
            runTime=runTime[arrayNum],
            density=density[arrayNum],
            cycleLengthDist=cycleLengthDist[,arrayNum],
            inheritGrowth=inheritGrowth[arrayNum],
            drugEffect=drugEffect[[arrayNum]],
            outputIncrement=outputIncrement[arrayNum],
            recordIncrement=recordIncrement[arrayNum],
            randSeed=randSeed[arrayNum],
            modelType=modelType[arrayNum],
            drugTime=drugTime[arrayNum],
            boundary=boundary[arrayNum],
            syncCycles=syncCycles[arrayNum],
            nG=nG[arrayNum],
            epsilon=epsilon[arrayNum]
            )


saveRDS(output, paste("output_", jobName, "_", arrayNum, ".rds", sep=""))
