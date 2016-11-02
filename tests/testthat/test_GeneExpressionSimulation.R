context("Testing Gene Expression Data Simulation")

test_that("pathway simulation - pooled", {

    ## create pathway
    pathway <- rexp(5, 1/20)
    names(pathway) <- letters[1:5]
    
    ## test G to S  
    gs <- simulatePathway(GE_testmod, pathway, 'S', sampFreq = 1,
                          sampSize = 1, singleCell = FALSE)
    expect_equal(gs[1,], pathway * 0.5)
    expect_equal(gs[2,], pathway * 0)

    ## test G to M
    gs <- simulatePathway(GE_testmod, pathway, 'M', sampFreq = 1,
                          sampSize = 1, singleCell = FALSE)
    expect_equal(gs[1,], pathway * 0)
    expect_equal(gs[2,], pathway * 0.5)

    ## test PROX
    gs <- simulatePathway(GE_testmod, pathway, 'PROX', sampFreq = 1,
                          sampSize = 1, singleCell = FALSE)
    expect_equal(gs[1,], pathway / 6)
    expect_equal(gs[2,], pathway / 6)

    ## test GROWTH
    gs <- simulatePathway(GE_testmod, pathway, 'GROWTH', sampFreq = 1,
                          sampSize = 1, singleCell = FALSE)
    expect_equal(gs[1,], pathway / 1.030972, tolerance=1e-3)
    expect_equal(gs[2,], pathway / 1.030972, tolerance=1e-3)

})

test_that("pathway simulation - single cell", {

    ## create pathway
    pathway <- rexp(5, 1/20)
    names(pathway) <- letters[1:5]

    ## test G to S  
    gs <- simulatePathway(GE_testmod, pathway, 'S', sampFreq = 1,
                          sampSize = 2, singleCell = TRUE)
    expect_equal(gs[1,], c(pathway * 1, pathway * 0))
    expect_equal(gs[2,], c(pathway * 0, pathway * 0))

    ## test G to M
    gs <- simulatePathway(GE_testmod, pathway, 'M', sampFreq = 1,
                          sampSize = 2, singleCell = TRUE)
    expect_equal(gs[1,], c(pathway * 0, pathway * 0))
    expect_equal(gs[2,], c(pathway * 1, pathway * 0))

    ## test PROX
    gs <- simulatePathway(GE_testmod, pathway, 'PROX', sampFreq = 1,
                          sampSize = 2, singleCell = TRUE)
    expect_equal(gs[1,], c(pathway / 6, pathway / 6))
    expect_equal(gs[2,], c(pathway / 6, pathway / 6))

    ## test GROWTH
    gs <- simulatePathway(GE_testmod, pathway, 'GROWTH', sampFreq = 1,
                          sampSize = 2, singleCell = TRUE)
    expect_equal(gs[1,], c(pathway / 1.035101, pathway / 1.026876),
                         tolerance=1e-3)
    expect_equal(gs[2,], c(pathway / 1.035101, pathway / 1.026876),
                         tolerance=1e-3)

})

