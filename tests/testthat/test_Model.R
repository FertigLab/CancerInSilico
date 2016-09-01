context("Model Output Testing")

test_that("no NA values in output", {
  
	output <- runCancerSim(5,1,0.1)
	expect_equal(sum(is.na(output@m_cells)), 0)

})
