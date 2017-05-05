library(CancerInSilico)
library(methods)

## generate sample models

modDefault <- runCellSimulation(10,10,0.1)
modLongRun <- runCellSimulation(5,100,0.1)
modLargeRun <- runCellSimulation(1000,1,0.1)
modHighDensity <- runCellSimulation(100,10,0.4)

save(modDefault, modLongRun, modLargeRun, modHighDensity,
    file = 'SampleModels.RData')

## generate sample pathways

growthExp <- function(model, cell, time)
{
    cycLength <- getCycleLength(model, time)[cell]
    return (1 / exp(cycLength / 24))
}

mitosisExp <- function(model, cell, time)
{
    window <- c(max(time - 1, 0), min(time + 1, model@runTime))

    a1 <- getAxisLength(model, window[1])[cell]
    a2 <- getAxisLength(model, window[2])[cell]

    return(ifelse(a2 < a1, 1, 0))
}

SPhaseExp <- function(model, cell, time)
{


}

contactInhibitionExp <- function(model, cell, time)
{


}

neighborsExp <- function(model, cell, time)
{


}

getGenes <- function(str) paste(str, letters[1:26], '_')

pwyGrowth <- new('Pathway', genes = getGenes('growth'), minExpression = 3,
    maxExpression = 5, expressionScale = growthExp)

pwyMitosis <- new('Pathway', genes = getGenes('mitosis'), minExpression = 0,
    maxExpression = 2, expressionScale = mitosisExp)

pwySPhase <- new('Pathway', genes = getGenes('Sphase'), minExpression = 0,
    maxExpression = 2, expressionScale = mitosisExp)

pwyContactInhibition <- new('Pathway', genes=getGenes('contact_inhibition'),
    minExpression = 0, maxExpression = 2, expressionScale = mitosisExp)

pwyNeighbors <- new('Pathway', genes = getGenes('neighbors'),
    minExpression = 0, maxExpression = 2, expressionScale = mitosisExp)

save(pwyGrowth, pwyMitosis, file = 'SamplePathways.RData')
