context("Model Output Testing")

test_that("no NA values in output", {
  
	output <- runCancerSim(5,1,0.1)
	expect_equal(sum(is.na(output@mCells)), 0)

})

test_that("CellModel getters", {

    mod <- runCancerSim(5,4,0.1, cycleLengthDist = 12)
    
    expect_equal(getCycleLengths(mod, 2)[1], 12)
    expect_equal(getParameters(mod)$runTime, 4)
    expect_equal(getParameters(mod)$initialNum, 5)
    expect_equal(getParameters(mod)$initialDensity, 0.1)
    show(mod)

})

test_that("getDrugEffect", {

    cyc_seq <- seq(8,16,1)
    de <- getDrugEffect(cycleLengthSeq = cyc_seq)
    expect_equal(length(de), 9)
    expect_equal(de[[9]][2], 0)
    expect_equal(length(de[[1]]), 2)

    f <- function(x) {return(0.4)}
    de <- getDrugEffect(f, cycleLengthSeq = cyc_seq)
    expect_equal(length(de), 9)
    expect_equal(de[[9]][2], 0.4)
    expect_equal(length(de[[1]]), 2)

    f <- function(x) {
        ret <- rnorm(100, (x - 8) / 16, 0.2)
        ret[ret > 1] <- 1
        ret[ret < 0] <- 0
        return (ret)

    }
    de <- getDrugEffect(f, cycleLengthSeq = cyc_seq)
    expect_equal(length(de), 9)
    expect_true(min(unlist(de)) >= 0)
    expect_true(max(de[[9]][2:101]) <= 1)
    expect_equal(length(de[[1]]), 101)

})

test_that("simulatePathway", {

    # cell models to test
    m1 <- runCancerSim(1,10)
    m2 <- runCancerSim(2,20) # goes to size 38
    m3 <- readRDS("testObj.rds") # how can I access this?

    # pathway
    p1 <- seq(0, 9, 1) # tests expressions of value integers 0-9
    zeros <- seq(0, 0, length.out = 10) # vector of all zeros, length 10

    # get gene expression matrix for m1, type PROX
    geMatrix_m1_PROX <- simulatePathway(m1, p1, type = PROX)
    expect_equal(geMatrix_PROX[0,], zeros) # no neighbors at t = 0
    expect_equal(geMatrix_PROX[1,], zeros) # no neighbors at t = 1
    for (i in 2:9) {
        expect_equal(geMatrix_PROX[i,], p1/6) # one neighbor at t = 2-9
    }
    # this 
    pathway <- rnorm(10,10,1)
    names(pathway) <- letters[1:10]

    # get gene expression mwwwatrix for m2, type PROX
    geMatrix_m2_PROX <- simulatePathway(m2, p1, type = PROX)

    # test gene expression matrix for m2
    for (j in 1:60) {
        # test gene expression at time i
   #     currExp_m2_PROX <- geMatrix_m2_PROX[j,]
   #     expect_equal(?,currExp_m2_PROX) # HOW TO CALCULATE EXPECTED?
    }


})
