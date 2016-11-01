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
    m1 <- runCancerSim(2,10)
    m3 <- readRDS("testObj.rds") # how can I access this?

    # pathway of length 10 with random values
    p1<- rnorm(10,10,1)
    names(pathway) <- letters[1:10]

    # get gene expression matrix for m1, type S
    geMatrix_m1_S <- simulatePathway(m2, p1, type = S)

    # check correct pathways for times 1-10
    for (t in 1:10) {
        # check correct num neighbors for indices 1-10
        for (i in 1:10) {
            cur_rad <- getRadii(m1, t)
            next_rad <- getRadii(m1, t + .1)
            transition <- ifelse(next_rad > sqrt(3/2) & cur_rad < sqrt(3/2), 1, 0)
            expect_equal(geMatrix_S[t,i], p1 * transition)
        }
    }

    # get gene expression matrix for m1, type M
    geMatrix_m1_M <- simulatePathway(m2, p1, type = M)

    # check correct pathways for times 1-10
    for (t in 1:10) {
        # check correct num neighbors for indices 1-10
        for (i in 1:10) {
            cur_ax <- getAxisLengthLength(m1, t)
            next_ax <- getAxisLength(m1, t + .1)
            transition <- ifelse(next_ax < cur_ax, 1, 0)
            expect_equal(geMatrix_S[t,i], p1 * transition)
        }
    }

    # get gene expression matrix for m1, type GROWTH
    geMatrix_m1_GROWTH <- simulatePathway(m2, p1, type = GROWTH)

    # check correct pathways for times 1-10
    for (t in 1:10) {
        # check correct num neighbors for indices 1-10
        for (i in 1:10) {
            # smooth out function
            cycle_len <- getCycleLengths(model,time)
            expect_equal(geMatrix_GROWTH[t,i], p1 * cycle_len)
        }
    }

    # get gene expression matrix for m1, type PROX
    geMatrix_m1_PROX <- simulatePathway(m1, p1, type = PROX)
 
    # check correct pathways for times 1-10
    for (t in 1:10) {
        # check correct num neighbors for indices 1-10
        for (i in 1:10) {
            neighbors <- getNumNeighbors(m1, t, i)
            expect_equal(geMatrix_PROX[t,i], p1 * (neighbors/6))
        }
    }
})
