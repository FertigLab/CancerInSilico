library(CancerInSilico)
library(methods)

## generate sample parameters

## generate sample models

simpleModel <- runCellSimulation(10,10,0.1)
longRun <- runCellSimulation(5,100,0.1)
largeRun <- runCellSimulation(1000,1,0.1)
highDensity <- runCellSimulation(100,10,0.4)

save(simpleModel, longRun, largeRun, highDensity, 
    file = 'SampleModels.RData')
