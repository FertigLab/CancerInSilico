context("Testing Gene Expression Data Simulation for Pathways")

test_that("getScalingFactor", {

    # test G to S scaling factor
    expect_equal(getScalingFactor(GE_testmod, 1:2, 0, 1, 'S'),
                c(TRUE, FALSE))
    expect_equal(getScalingFactor(GE_testmod, 1:2, 1, 1, 'S'),
                c(FALSE, FALSE))
    
    # test G to M scaling factor
    expect_equal(getScalingFactor(GE_testmod, 1:2, 0, 1, 'M'),
                c(FALSE, FALSE))
    expect_equal(getScalingFactor(GE_testmod, 1:2, 1, 1, 'M'),
                c(TRUE, FALSE))

    # test PROX scaling factor
    expect_equal(getScalingFactor(GE_testmod, 1:2, 0, 1, 'PROX'),
                c(1 / 6, 1 / 6))
    expect_equal(getScalingFactor(GE_testmod, 1:2, 1, 1, 'PROX'),
                c(1 / 6, 1 / 6))

    # test GROWTH scaling factor
    expect_equal(getScalingFactor(GE_testmod, 1:2, 0, 1, 'GROWTH'),
                c(0.966, 0.974), tolerance = 1e-3)
    expect_equal(getScalingFactor(GE_testmod, 1:2, 1, 1, 'GROWTH'),
                c(0.966, 0.974), tolerance = 1e-3)

})

test_that("simulatePathway - pooled", {

    ## create pathway
    pathway <- list()    
    pathway[["genes"]] <- letters[1:5]
    pathway[["min"]] <- rexp(5,1/4)
    pathway[["max"]] <- pathway[["min"]] + 4
    
    # for readability
    base_exp <- pathway[["min"]]
    
    ## test G to S  
    gs <- unname(simulatePathway(GE_testmod, pathway, 'S', sampFreq = 1,
                          sampSize = 1, singleCell = FALSE))
    expect_equal(gs[,1], 4 * 0.5 + base_exp)
    expect_equal(gs[,2], 4 * 0.0 + base_exp)

    ## test G to M
    gs <- unname(simulatePathway(GE_testmod, pathway, 'M', sampFreq = 1,
                          sampSize = 1, singleCell = FALSE))
    expect_equal(gs[,1], 4 * 0.0 + base_exp)
    expect_equal(gs[,2], 4 * 0.5 + base_exp)

    ## test PROX
    gs <- unname(simulatePathway(GE_testmod, pathway, 'PROX', sampFreq = 1,
                          sampSize = 1, singleCell = FALSE))
    expect_equal(gs[,1], 4 / 6 + base_exp)
    expect_equal(gs[,2], 4 / 6 + base_exp)

    ## test GROWTH
    gs <- unname(simulatePathway(GE_testmod, pathway, 'GROWTH', sampFreq=1,
                          sampSize = 1, singleCell = FALSE))
    expect_equal(gs[,1], 4 / 1.03 + base_exp, tolerance=1e-3)
    expect_equal(gs[,2], 4 / 1.03 + base_exp, tolerance=1e-3)

})

test_that("simulatePathway - single cell", {

    ## create pathway
    pathway <- list()    
    pathway[["genes"]] <- letters[1:5]
    pathway[["min"]] <- rexp(5,1/4)
    pathway[["max"]] <- pathway[["min"]] + 4
    
    # for readability
    base_exp <- pathway[["min"]]

    ## test G to S  
    gs <- unname(simulatePathway(GE_testmod, pathway, 'S', sampFreq = 1,
                          sampSize = 2, singleCell = TRUE))
    expect_equal(gs[,1], 4 * 1 + base_exp)
    expect_equal(gs[,2], 4 * 0 + base_exp)
    expect_equal(gs[,3], 4 * 0 + base_exp)
    expect_equal(gs[,4], 4 * 0 + base_exp)

    ## test G to M
    gs <- unname(simulatePathway(GE_testmod, pathway, 'M', sampFreq = 1,
                          sampSize = 2, singleCell = TRUE))
    expect_equal(gs[,1], 4 * 0 + base_exp)
    expect_equal(gs[,2], 4 * 0 + base_exp)
    expect_equal(gs[,3], 4 * 1 + base_exp)
    expect_equal(gs[,4], 4 * 0 + base_exp)

    ## test PROX
    gs <- unname(simulatePathway(GE_testmod, pathway, 'PROX', sampFreq = 1,
                          sampSize = 2, singleCell = TRUE))
    expect_equal(gs[,1], 4 / 6 + base_exp)
    expect_equal(gs[,2], 4 / 6 + base_exp)
    expect_equal(gs[,3], 4 / 6 + base_exp)
    expect_equal(gs[,4], 4 / 6 + base_exp)

    ## test GROWTH
    gs <- unname(simulatePathway(GE_testmod, pathway, 'GROWTH', sampFreq=1,
                          sampSize = 2, singleCell = TRUE))
    expect_equal(gs[,1], 4 / 1.035 + base_exp, tolerance = 1e-3)
    expect_equal(gs[,2], 4 / 1.027 + base_exp, tolerance = 1e-3)
    expect_equal(gs[,3], 4 / 1.035 + base_exp, tolerance = 1e-3)
    expect_equal(gs[,4], 4 / 1.027 + base_exp, tolerance = 1e-3)

})

test_that("getPathways - default value", {

    pwys <- getPathways()
    expect_equal(pwys, inSilicoPathways)

})

test_that("getPathways - custom pathway", {

    # list of pathways
    pwys <- list()

    # invalid pathway
    pwys[["invalid"]] <- inSilicoPathways[["GtoM"]]

    # should quit with error
    expect_error(getPathways(pwys))
    
    # add valid pathway to pwys
    pwys[["Prox"]] <- inSilicoPathways[["Prox"]]

    # should only throw a warning now
    expect_warning(getPathways(pwys))

    # rename invalid pathway
    names(pwys) <- c("GtoM", "Prox")

    # no warnings or errors now
    expect_error(getPathways(pwys), regexp = NA)
    expect_warning(getPathways(pwys), regexp = NA)

})

test_that("checkReferenceDataset", {

    # create valid dataset
    data <- replicate(5, runif(5,10,20)) 

    # no warnings or errors
    expect_error(checkReferenceDataset(data), regexp = NA)
    expect_warning(checkReferenceDataset(data), regexp = NA)

    # add negative entry
    data[1,1] <- -1

    # warning for negative number
    expect_warning(checkReferenceDataset(data))

    # add large entry
    data[1,1] <- 51

    # warning for non log transformed data
    expect_warning(checkReferenceDataset(data))

})








