context("Testing functions in RunCellSimulation.R")

test_that("no NA values in output",
{  
	output <- inSilicoCellModel(3, 100, 0.1, nG = 4,outputIncrement=100)
	expect_equal(sum(is.na(output@cells)), 0)

	output <- inSilicoCellModel(1000, 1, 0.1, nG = 4)
	expect_equal(sum(is.na(output@cells)), 0)
})
