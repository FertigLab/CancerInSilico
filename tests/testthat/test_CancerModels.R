context("Model Output Testing")

test_that("no NA values in output", {
  
	output <- runCancerSim(5,3,0.1, nG = 4)
	expect_equal(sum(is.na(output@mCells)), 0)

})

test_that("CellModel getters", {

    expect_equal(getCycleLengths(CM_test_model, 1)[1], 12)
    expect_equal(getParameters(CM_test_model)$runTime, 10)
    expect_equal(getParameters(CM_test_model)$initialNum, 10)
    expect_equal(getParameters(CM_test_model)$initialDensity, 0.1)

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

