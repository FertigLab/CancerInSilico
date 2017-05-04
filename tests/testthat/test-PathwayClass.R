context('Testing Pathway Class')

test_that('Pathway Initialization',
{  
    expect_error(new('Pathway'), 'missing genes')
    expect_error(new('Pathway', genes=c('a'), minExpression=0,
        maxExpression=2), 'missing expressionScale')
    expect_error(new('Pathway', genes=c('a'),
        expressionScale = function(x,y,z) x), 'missing \'maxExpression\'')
    expect_error(new('Pathway', genes=c('a'), maxExpression=0,
        expressionScale = function(x,y,z) x), 'missing \'minExpression\'')
    expect_error(new('Pathway', genes = c('a'), minExpression = 0,
        maxExpression=2, expressionScale = function(x,y,z) x), NA)   
})

growthExp <- function(model, cell, time)
{
    cycLength <- getCycleLengths(model, time)[cell]
    return (1/exp(cycLength/20))
}

mitosisExp <- function(model, cell, time)
{
    window <- c(max(time - 1,0), min(time + 1, model@runTime))

    a1 <-getAxisLength(model, window[1])[cell]
    a2 <-getAxisLength(model, window[2])[cell]

    if (a2 < a1) return(1)
    else return(0)
}

test_that('Pathway Expression Function',
{
    pwyGrowth <- new('Pathway', genes=letters[1:26], minExpression=3,
        maxExpression=5, expressionScale = growthExp)

    pwyMitosis <- new('Pathway', genes=LETTERS[1:26], minExpression=0,
        maxExpression=2, expressionScale = mitosisExp)

#    print(getNumberOfCells(simpleModel, 5))
#    pwyGrowth@expressionScale(simpleModel, 5, 5)
})

test_that('Full Expression Simulation',
{

})
