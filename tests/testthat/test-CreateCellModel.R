context('Testing CellModel Creation')

test_that('cell model initialization - needed parameters',
{  
    expect_error(new('DrasdoHohmeModel'), paste('argument "initialNum" is',
        'missing, with no default'))
    expect_error(new('DrasdoHohmeModel', initialNum = 1), paste('argument',
        '"runTime" is missing, with no default'))
    expect_error(new('DrasdoHohmeModel', initialNum = 1, runTime = 1),
        'argument "density" is missing, with no default')
    expect_error(new('DrasdoHohmeModel', initialNum = 1, runTime = 1,
        density = 0.02), NA)
})

test_that('cell model initialization - validObject is called',
{
    expect_error(new('DrasdoHohmeModel', initialNum = 1, runTime = 1,
        density = 0.75), "'density' too high for initial seeding process")
})

test_that('object validity check',
{
    testMod <- new('DrasdoHohmeModel', initialNum = 1, runTime = 1,
        density = 0.01)

    expect_error(validObject(testMod), NA)
    
    testMod@timeIncrement <- -1
    expect_error(validObject(testMod), paste("'timeIncrement' must be",
        "greater than zero"))
    testMod@timeIncrement <- 1

    testMod@nG <- -1
    expect_error(validObject(testMod), "'nG' must be at least one")
})

test_that('default arguments',
{
    testMod <- new('DrasdoHohmeModel', initialNum = 1, runTime = 1,
        density = 0.01)

    expect_equal(testMod@initialNum, 1)
    expect_equal(testMod@runTime, 1)
    expect_equal(testMod@density, 0.01)
    expect_equal(testMod@boundary, 1)
    expect_equal(testMod@syncCycles, FALSE)
    expect_equal(testMod@randSeed, 0)
    expect_equal(testMod@outputIncrement, 4)
    expect_equal(testMod@recordIncrement, 0.1)
    expect_equal(testMod@timeIncrement, 0.00589, tolerance = 1e-3)
    expect_equal(testMod@cellTypes[[1]]@name, 'DEFAULT')
    expect_equal(testMod@cellTypes[[1]]@size, 1)
    expect_equal(testMod@cellTypes[[1]]@minCycle, 24)
    expect_false(is.null(testMod@cellTypes[[1]]@cycleLength()))
    expect_equal(testMod@cellTypeInitFreq[1], 1)
    expect_equal(testMod@maxTranslation, 0.1)
    expect_equal(testMod@maxRotation, 0.05, tolerance = 1e-3)
    expect_equal(testMod@nG, 28)
    expect_equal(testMod@epsilon, 10)
    expect_equal(testMod@delta, 0.2)
})

test_that('drasdo parameter calculations',
{
    testMod <- new('DrasdoHohmeModel', initialNum = 1, runTime = 1,
        density = 0.01, nG = 10, delta = 0.4, maxTranslation = 1,
        maxRotation = 1, timeIncrement = 1)
    
    expect_equal(testMod@nG, 10)
    expect_equal(testMod@delta, 0.4)
    expect_equal(testMod@maxTranslation, 0.2)
    expect_equal(testMod@maxRotation, 0.1, tolerance = 1e-3)
    expect_equal(testMod@timeIncrement, 0.0333, tolerance = 1e-3)
})



