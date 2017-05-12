library('CancerInSilico')
library(methods)

args <- commandArgs(TRUE)
arrayNum <- as.integer(args[1])
jobName <- args[2]

#### Set Defaults ####

boundary <- 1
syncCycles <- FALSE
randSeed <- 0 
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

typeA_dist <- seq(from = 0, to = 1, length.out = 11)
gr_AtoB_rat_list <- seq(from = 0, to = 2, length.out = 9)

# Assuming Indexing 1 to 99

mtypeA_dist <- typeA_dist[floor((arrayNum- 1)/9) + 1]
mgRate <- gr_AtoB_rat_list[((arrayNum - 1) %% 9) + 1]

# Create 2 cell types : A, B

ctA <- new('CellType', name='A', cycleLength <- function() {return(48)})
ctB <- new('CellType', name='B', cycleLength <- function() {return(48*mgRate)})

# Set changed params

cellTypes <- c(ctA, ctB)
cellTypeInitFreq <- c(mtypeA_dist, 1 - mtypeA_dist)

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
                            cellTypeInitFreq=cellTypeInitFreq,
                            drugs=drugs,
                            maxDeformation=maxDeformation,
                            maxTranslation=maxTranslation,
                            maxRotation=maxRotation,
                            nG=nG,
                            epsilon=epsilon,
                            delta=delta
                            )       

save(output, paste("output_", jobName, "_", arrayNum, ".RData", sep=""))
