context("Model Output Testing")

test_that("no NA values in output", {
  
	output <- runModel(100,250,0.1)
	expect_equal(sum(is.na(output)), 0)
	expect_true(checkCellOverlap(output,250))

})
