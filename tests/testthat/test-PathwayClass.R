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

test_that('Pathway Expression Function',
{
    expect_equal(pwyMitosis@expressionScale(modDefault, 9, 5), 1)
    expect_equal(pwyGrowth@expressionScale(modDefault, 5, 5), 0.14,
        tolerance=0.01)
})

test_that('Pathway Expression Simulation',
{
    meanExp <- simulateExpression(pwyMitosis, modDefault, 1)
    expect_true(all(meanExp[,3] < meanExp[,4]))

    meanExpSS <- simulateExpression(pwyMitosis, modDefault, 1, TRUE, 10)
    expect_true(all(meanExpSS[,'c9_t2'] < meanExpSS[,'c9_t3']))
})
