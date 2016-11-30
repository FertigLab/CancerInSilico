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

