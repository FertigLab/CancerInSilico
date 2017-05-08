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
    window <- c(max(time - 1, 0), min(time + 1, model@runTime))

    r1 <- getRadii(model, window[1])[cell]
    r2 <- getRadii(model, window[2])[cell]

    type <- getCellTypes(model, time)[cell]

    return (ifelse(r1 < sqrt(1.5 * type@size) & r2 > sqrt(1.5 * type@size),
        1, 0))
}

contactInhibitionExp <- function(model, cell, time)
{
    return(getContactInhibition(model, time)[cell])
}

neighborsExp <- function(model, cell, time)
{
    num <- getNumberOfNeighbors(model, time, cell, 5)
    return(num / 6)
}

getGenes <- function(str) paste(str, letters[1:26], '_')

pwyGrowth <- new('Pathway', genes = getGenes('growth'), minExpression = 3,
    maxExpression = 5, expressionScale = growthExp)

pwyMitosis <- new('Pathway', genes = getGenes('mitosis'), minExpression = 0,
    maxExpression = 2, expressionScale = mitosisExp)

pwySPhase <- new('Pathway', genes = getGenes('Sphase'), minExpression = 0,
    maxExpression = 2, expressionScale = SPhaseExp)

pwyContactInhibition <- new('Pathway', genes=getGenes('contact_inhibition'),
    minExpression = 0, maxExpression = 2,
    expressionScale = contactInhibitionExp)

pwyNeighbors <- new('Pathway', genes = getGenes('neighbors'),
    minExpression = 0, maxExpression = 2, expressionScale = neighborsExp)

save(pwyGrowth, pwyMitosis, pwySPhase, pwyContactInhibition, pwyNeighbors,
    file = 'SamplePathways.RData')



