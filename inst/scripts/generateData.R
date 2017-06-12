library(CancerInSilico)
library(methods)

## Full run examples

#modDefault <- runCellSimulation(10,10,0.1)
#modLongRun <- runCellSimulation(5,100,0.1)
#modLargeRun <- runCellSimulation(1000,1,0.1)
#modHighDensity <- runCellSimulation(100,10,0.4)

## CellType examples

c1 <- new('CellType', name='DEFAULT', size=1, minCycle=48,
    cycleLength=function() 48)

c2 <- new('CellType', name='DOUBLE_SIZE', size=2, minCycle=48,
    cycleLength=function() 48)

c3 <- new('CellType', name='SHORT_CYCLE', size=1, minCycle=6,
    cycleLength=function() 6)

c4 <- new('CellType', name='RANDOM_CYCLE', size=1, minCycle=10,
    cycleLength=function() runif(1,10,20))

modCellTypes <- new('DrasdoHohmeModel', initialNum=1, runTime=1,
    density=0.1, cellTypes=c(c1,c2,c3,c4), cellTypeInitFreq=rep(0.25,4))

## Drug examples

d1 <- new('Drug', name='NO_EFFECT', timeAdded=0, cycleLengthEffect=
    function(type,len) len)

d2 <- new('Drug', name='HALF_CYCLE_LENGTH', timeAdded=0, cycleLengthEffect=
    function(type,len) len/2)

d3 <- new('Drug', name='HALF_DEFAULT_TYPE', timeAdded=0, cycleLengthEffect=
    function(type,len) ifelse(type=='DEFAULT', len/2, len))

d4 <- new('Drug', name='ADD_LATE', timeAdded=6, cycleLengthEffect=
    function(type,len) len)

modDrugs <- new('DrasdoHohmeModel', initialNum=100, runTime=4, density=0.4,
    randSeed=0, outputIncrement=4, recordIncrement=0.1, syncCycle=FALSE,
    cellTypes=c(c1,c2), cellTypeInitFreq=c(0.5,0.5), drugs=c(d1,d2,d3,d4))

save(modDefault, modLongRun, modLargeRun, modHighDensity, modCellTypes,
    modDrugs, file = 'SampleModels.RData')

## Pathway examples

growthExp <- function(model, cell, time)
{
    cycLength <- getCycleLengths(model, time)[cell]
    return(1 / exp(cycLength / 24))
}

mitosisExp <- function(model, cell, time)
{
    window <- c(max(time - 2, 0), min(time + 2, model@runTime))

    a1 <- getAxisLength(model, window[1])[cell]
    a2 <- getAxisLength(model, window[2])[cell]
    if (is.na(a1)) a1 <- 0


    return(ifelse(a2 < a1, 1, 0))
}

SPhaseExp <- function(model, cell, time)
{
    window <- c(max(time - 1, 0), min(time + 1, model@runTime))

    r1 <- getRadii(model, window[1])[cell]
    r2 <- getRadii(model, window[2])[cell]

    typeNdx <- getCellTypes(model, time)[cell]
    type <- model@cellTypes[[typeNdx]]

    return(ifelse(r1 < sqrt(1.5 * type@size) & r2 > sqrt(1.5 * type@size),
        1, 0))
}

contactInhibitionExp <- function(model, cell, time)
{
    return(getContactInhibition(model, time)[cell])
}

getGenes <- function(str) paste(str, letters[1:26], sep='_')
getBound <- function(mu) runif(26,mu,mu+2)

pwyGrowth <- new('Pathway', genes = getGenes('growth'), minExpression = getBound(1),
    maxExpression = getBound(5), expressionScale = growthExp)

pwyMitosis <- new('Pathway', genes = getGenes('mitosis'), minExpression = getBound(3),
    maxExpression = getBound(8), expressionScale = mitosisExp)

pwySPhase <- new('Pathway', genes = getGenes('Sphase'), minExpression = getBound(0),
    maxExpression = getBound(3), expressionScale = SPhaseExp)

pwyContactInhibition <- new('Pathway', genes=getGenes('contact_inhibition'),
    minExpression = getBound(1), maxExpression = getBound(3),
    expressionScale = contactInhibitionExp)

save(pwyGrowth, pwyMitosis, pwySPhase, pwyContactInhibition, pwyNeighbors,
    file = 'SamplePathways.RData')
