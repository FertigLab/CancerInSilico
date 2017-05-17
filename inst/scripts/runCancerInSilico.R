library('CancerInSilico')
library(methods)

# max arrayNum = 210

args <- commandArgs(TRUE)
arrayNum <- as.integer(args[1])
jobName <- args[2]

#### Set Defaults ####

initialNum <- 80
runTime <- 168
density <- 0.2
boundary <- 1
syncCycles <- FALSE
randSeed <- 0
outputIncrement <- 4
recordIncrement <- 0.25
timeIncrement <- 0.001
cellTypes <- c(new('CellType', name='DEFAULT'))
cellTypeInitFreq <- c(1)
drugs <- list()
maxDeformation <- 0.1
maxTranslation <- 0.1
maxRotation <- 0.3
nG <- 28
epsilon <- 10.0
delta <- 0.2

#### Set Custom Values ####

allInitFreq_typeA <- seq(0,1,0.05)
allCellTypes_typeB <- lapply(seq(12,48,4), function(l) new('CellType', name='B', minCycle=l, cycleLength=function() l))

dim <- c(length(allInitFreq_typeA), length(allCellTypes_typeB))
indexArray <- array(1:prod(dim), dim)
index <- which(indexArray==arrayNum, arr.ind=TRUE)

initFreq_typeA <- allInitFreq_typeA[index[1]]
ctB <- allCellTypes_typeB[index[2]]
ctA <- new('CellType', name='A', cycleLength=function() {return(24)}, minCycle=24)

# Set changed params

cellTypes <- c(ctA, ctB)
cellTypeInitFreq <- c(initFreq_typeA, 1 - initFreq_typeA)

# Sanity check for distribution and cycle lengths

print(paste("Cell Type dist:", "A:", cellTypeInitFreq[1], "B:", cellTypeInitFreq[2]))
print(paste("Cell Type B minCycle:", slot(cellTypes[2][[1]], "minCycle")))


#### Run Simulation ####

output <- runCellSimulation(initialNum=initialNum,
                            runTime=runTime,
                            density=density,
                            boundary=boundary,
                            syncCycles=syncCycles,
                            randSeed=randSeed,
                            modelType=modelType,
                            outputIncrement=outputIncrement,
                            recordIncrement=recordIncrement,
                            timeIncrement=timeIncrement,
                            cellTypes=cellTypes,
                            cellTypeInitFreq=cellTypeInitFreq,
                            drugs=drugs,
                            maxDeformation=maxDeformation,
                            maxTranslation=maxTranslation,
                            maxRotation=maxRotation,
                            nG=nG,
                            epsilon=epsilon,
                            delta=delta
                            )       

save(output, file=paste("output_", jobName, "_", arrayNum, ".RData", sep=""))
