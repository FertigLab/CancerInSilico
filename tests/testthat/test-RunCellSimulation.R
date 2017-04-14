context("Testing functions in RunCellSimulation.R")

test_that("no NA values in output",
{  
	output <- runCellSimulation(5, 3, 0.1, nG = 4)
	expect_equal(sum(is.na(output@cells)), 0)
})
