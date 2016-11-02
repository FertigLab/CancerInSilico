context("Testing Gene Expression Data Simulation")

test_that("pathway simulation - pooled", {

    ## create pathway

    set.seed(0)
    pathway <- rexp(5, 20)
    names(pathway) <- letters[1:5]
    
    ## test G to S  
    gs <- simulatePathway(GE_testmod, pathway, 'S')
    expect_equal(gs[1,], pathway * 0.5)
    expect_equal(gs[2,], pathway * 0)

    ## test G to M
    gs <- simulatePathway(GE_testmod, pathway, 'M')
    expect_equal(gs[1,], pathway * 0)
    expect_equal(gs[2,], pathway * 0.5)

    ## test PROX
    gs <- simulatePathway(GE_testmod, pathway, 'PROX')
    expect_equal(gs[1,], pathway / 6)
    expect_equal(gs[2,], pathway / 6)

    ## test GROWTH
    gs <- simulatePathway(GE_testmod, pathway, 'GROWTH')
    expect_equal(gs[1,], pathway / 1.030972, tolerance=1e-3)
    expect_equal(gs[2,], pathway / 1.030972, tolerance=1e-3)

})

