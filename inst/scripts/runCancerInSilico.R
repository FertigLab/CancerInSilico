library(CancerInSilico)

args <- commandArgs(TRUE)
arrayNum <- as.integer(args[1])
jobName <- args[2]

#### Set Defaults ####

boundary <- 1
syncCycles <- FALSE
randSeed <- 0, 
outputIncrement <- 4
recordIncrement <- 0.1
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


#### Run Simulation ####

output <- runCellSimulation(initialNum=initialNum,
                            runTime=runTime,
                            density=density,
                            boundary=boundary,
                            syncCycles=syncCycles,
                            randSeed=randSeed,
                            modelType=modelType,
                            outputIncrement=outputIncrement,
                            recordIncrement=timeIncrement,
                            timeIncrement=timeIncrement,
                            cellTypes=cellTypes,
                            cellTypeInitFreq,
                            drugs,
                            maxDeformation=maxDeformation,
                            maxTranslation=maxTranslation,
                            maxRotation=maxRotation,
                            nG=nG,
                            epsilon=epsilon,
                            delta=delta
                            )       

save(output, paste("output_", jobName, "_", arrayNum, ".RData", sep=""))
