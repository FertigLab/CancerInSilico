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

test_that('Pathway Expression Simulation',
{
    pwyMitosis@minExpression <- runif(length(pwyMitosis@genes), 2,4)
    pwyMitosis@maxExpression <- runif(length(pwyMitosis@genes), 6,8)
    activity <- simulatePathwayActivity(pwyMitosis, modDefault, 1, 10)
    meanExp <- simulatePathwayExpression(pwyMitosis, activity)
    expect_true(all(meanExp[,2] < meanExp[,3]))

    activitySS <- simulatePathwayActivity(pwyMitosis, modDefault, 1, 10, T)
    meanExpSS <- simulatePathwayExpression(pwyMitosis, activitySS)
    expect_true(all(meanExpSS[,'c9_t1'] < meanExpSS[,'c9_t2']))
})

test_that('Calibrate Pathways',
{  
    pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
    expect_equal(pwyMitosis@minExpression[1], 4.78, tolerance=0.01)
    expect_equal(pwyMitosis@maxExpression[1], 5.55, tolerance=0.01)
})

