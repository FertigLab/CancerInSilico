context("Testing functions in CellModel-class.R")

test_that("timeToRow", {
  
    expect_equal(timeToRow(PS_test_model, 0), 1)
    expect_equal(timeToRow(PS_test_model, 2), 3)
    expect_equal(timeToRow(PS_test_model, 1.9), 2)

})

test_that("getColumn", {

    expect_equal(getColumn(PS_test_model, 0, 1), c(0, 2))
    expect_equal(getColumn(PS_test_model, 1, 1), c(0, 2))
    expect_equal(getColumn(PS_test_model, 2, 1), c(1, 2))

    expect_equal(getColumn(PS_test_model, 1, 2), c(0, 2))
    expect_equal(getColumn(PS_test_model, 1, 3), c(1.3, 1))
    expect_equal(getColumn(PS_test_model, 1, 4), c(2.6, 2))
    expect_equal(getColumn(PS_test_model, 1, 5), c(0, 0))
    expect_equal(getColumn(PS_test_model, 1, 6), c(0.01, 0.03))

})

test_that("getCoordinates", {

    expect_equal(getCoordinates(PS_test_model, 0)[1,], c(0,0))
    expect_equal(getCoordinates(PS_test_model, 0)[2,], c(2,2))

})

test_that("getRadii", {

    expect_equal(getRadii(PS_test_model, 1), c(1.3, 1))

})

test_that("getAxisLength", {

    expect_equal(getAxisLength(PS_test_model, 1), c(2.6, 2))

})

test_that("getAxisAngle", {

    expect_equal(getAxisAngle(PS_test_model, 2), c(0, 0))

})

test_that("getGrowthRates", {

    expect_equal(getGrowthRates(PS_test_model, 2), c(0.03, 0.03))

})

test_that("getCycleLengths", {

    expect_equal(getCycleLengths(PS_test_model, 2), c(6.3, 6.3),
        tolerance = 0.01)

})

test_that("getNumberOfCells", {

    expect_equal(getNumberOfCells(PS_test_model, 2), 2)

})

test_that("getDensity", {

    expect_equal(getDensity(PS_test_model, 2), 0.02)

})

test_that("getNumNeighbors", {


    expect_equal(getNumNeighbors(PS_test_model, 2, 1), 1)    

})




