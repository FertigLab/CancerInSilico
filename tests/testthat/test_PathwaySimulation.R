context("Testing functions in PathwaySimulation.R")

test_that("getScalingFactor", {

    # test G to S scaling factor
    expect_equal(getScalingFactor(PS_test_model, 1:2, 0, 1, 'GtoS'),
                c(TRUE, FALSE))
    expect_equal(getScalingFactor(PS_test_model, 1:2, 1, 1, 'GtoS'),
                c(FALSE, FALSE))
    
    # test G to M scaling factor
    expect_equal(getScalingFactor(PS_test_model, 1:2, 0, 1, 'GtoM'),
                c(FALSE, FALSE))
    expect_equal(getScalingFactor(PS_test_model, 1:2, 1, 1, 'GtoM'),
                c(TRUE, FALSE))

    # test PROX scaling factor
    expect_equal(getScalingFactor(PS_test_model, 1:2, 0, 1, 'Prox'),
                c(1 / 6, 1 / 6))
    expect_equal(getScalingFactor(PS_test_model, 1:2, 1, 1, 'Prox'),
                c(1 / 6, 1 / 6))

    # test GROWTH scaling factor
    expect_equal(getScalingFactor(PS_test_model, 1:2, 0, 1, 'Growth'),
                c(0.966, 0.974), tolerance = 1e-3)
    expect_equal(getScalingFactor(PS_test_model, 1:2, 1, 1, 'Growth'),
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
    gs <- unname(simulatePathway(PS_test_model, pathway, 'GtoS',
                    sampFreq = 1, sampSize = 1, singleCell = FALSE))
    expect_equal(gs[,1], 4 * 0.5 + base_exp)
    expect_equal(gs[,2], 4 * 0.0 + base_exp)

    ## test G to M
    gs <- unname(simulatePathway(PS_test_model, pathway, 'GtoM',
                    sampFreq = 1, sampSize = 1, singleCell = FALSE))
    expect_equal(gs[,1], 4 * 0.0 + base_exp)
    expect_equal(gs[,2], 4 * 0.5 + base_exp)

    ## test PROX
    gs <- unname(simulatePathway(PS_test_model, pathway, 'Prox',
                    sampFreq = 1, sampSize = 1, singleCell = FALSE))
    expect_equal(gs[,1], 4 / 6 + base_exp)
    expect_equal(gs[,2], 4 / 6 + base_exp)

    ## test GROWTH
    gs <- unname(simulatePathway(PS_test_model, pathway, 'Growth',
                    sampFreq=1, sampSize = 1, singleCell = FALSE))
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
    gs <- unname(simulatePathway(PS_test_model, pathway, 'GtoS',
                    sampFreq = 1, sampSize = 2, singleCell = TRUE))
    expect_equal(gs[,1], 4 * 1 + base_exp)
    expect_equal(gs[,2], 4 * 0 + base_exp)
    expect_equal(gs[,3], 4 * 0 + base_exp)
    expect_equal(gs[,4], 4 * 0 + base_exp)

    ## test G to M
    gs <- unname(simulatePathway(PS_test_model, pathway, 'GtoM',
                    sampFreq = 1, sampSize = 2, singleCell = TRUE))
    expect_equal(gs[,1], 4 * 0 + base_exp)
    expect_equal(gs[,2], 4 * 0 + base_exp)
    expect_equal(gs[,3], 4 * 1 + base_exp)
    expect_equal(gs[,4], 4 * 0 + base_exp)

    ## test PROX
    gs <- unname(simulatePathway(PS_test_model, pathway, 'Prox',
                    sampFreq = 1, sampSize = 2, singleCell = TRUE))
    expect_equal(gs[,1], 4 / 6 + base_exp)
    expect_equal(gs[,2], 4 / 6 + base_exp)
    expect_equal(gs[,3], 4 / 6 + base_exp)
    expect_equal(gs[,4], 4 / 6 + base_exp)

    ## test GROWTH
    gs <- unname(simulatePathway(PS_test_model, pathway, 'Growth',
                    sampFreq=1, sampSize = 2, singleCell = TRUE))
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

test_that("checkReferenceDataSet", {

    # create dummy gene names
    genes <- letters[1:10]

    # create valid dataset
    data <- replicate(5, runif(15,10,20)) 
    row.names(data) <- letters[1:15]

    # no warnings or errors
    expect_error(checkReferenceDataSet(data, genes), regexp = NA)
    expect_warning(checkReferenceDataSet(data, genes), regexp = NA)

    # error for gene not in data set
    invalid_genes <- letters[1:16]
    expect_error(checkReferenceDataSet(data, invalid_genes))

    # no error for NA values for non-pathway gene
    data[12,] <- rep(NA, 5)
    expect_error(checkReferenceDataSet(data, genes), regexp = NA)
    
    # error for NA values for pathway gene
    data[8,] <- rep(NA, 5)
    expect_error(checkReferenceDataSet(data, genes))
    data[8,] <- rep(10, 5)

    # error for negative number
    data[1,1] <- -1
    expect_error(checkReferenceDataSet(data, genes))

    # warning for non log transformed data
    data[1,1] <- 51
    expect_warning(checkReferenceDataSet(data, genes))

})

test_that("getPathwayExpressionRange", {

    # create valid dataset
    data <- replicate(5, runif(15,10,20)) 
    row.names(data) <- letters[1:15]

    # create pathways
    pwys <- list()
    pwys[["A"]][["genes"]] <- letters[1:10]
    pwys[["B"]][["genes"]] <- letters[5:15]    

    # set range without using reference data set
    pwys <- setPathwayExpressionRange(pathways = pwys)

    # check that min < max for each gene
    expect_equal(sum(pwys[["A"]][["min"]] > pwys[["A"]][["max"]]), 0)
    expect_equal(sum(pwys[["B"]][["min"]] > pwys[["B"]][["max"]]), 0)

    # check length of min/max
    expect_equal(length(pwys[["A"]][["min"]]), 10)
    expect_equal(length(pwys[["A"]][["max"]]), 10)
    expect_equal(length(pwys[["B"]][["min"]]), 11)
    expect_equal(length(pwys[["B"]][["max"]]), 11)
    
    # reset pathways
    pwys <- list()
    pwys[["A"]][["genes"]] <- letters[1:10]
    pwys[["B"]][["genes"]] <- letters[5:15]    

    # set range using reference data set
    expect_error(pwys <- setPathwayExpressionRange(pathways = pwys,
                            ReferenceDataSet = data), regexp = NA)

    # check that min < max for each gene
    expect_equal(sum(pwys[["A"]][["min"]] > pwys[["A"]][["max"]]), 0)
    expect_equal(sum(pwys[["B"]][["min"]] > pwys[["B"]][["max"]]), 0)

    # check length of min/max
    expect_equal(length(pwys[["A"]][["min"]]), 10)
    expect_equal(length(pwys[["A"]][["max"]]), 10)
    expect_equal(length(pwys[["B"]][["min"]]), 11)
    expect_equal(length(pwys[["B"]][["max"]]), 11)

})







