context('Testing Pathway Class')

test_that('Pathway Initialization',
{  
    expect_error(new('Pathway'), 'missing genes')
    expect_error(new('Pathway', genes=c('a'), minExpression=0,
        maxExpression=2), 'missing expressionScale')
    expect_error(new('Pathway', genes = c('a'), minExpression = 0,
        maxExpression=2, expressionScale = function(x,y,z) x), NA)   
})

test_that('Pathway Expression Function',
{
    expect_equal(pwyMitosis@expressionScale(modDefault, 9, 5), 1)
    expect_equal(pwyGrowth@expressionScale(modDefault, 5, 5), 0.607,
        tolerance=0.01)
})

test_that('Calibrate Pathways',
{  
    pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
    expect_equal(pwyMitosis@minExpression[1], 4.78, tolerance=0.01)
    expect_equal(pwyMitosis@maxExpression[1], 5.55, tolerance=0.01)
})

