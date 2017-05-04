context('Testing CellModel Access')

test_that('cell model getters',
{
    data(SampleModels)  

    expect_equal(getNumberOfCells(simpleModel,      0),   10)
    expect_equal(getNumberOfCells(longRun,          0),    5)
    expect_equal(getNumberOfCells(largeRun,         0), 1000)
    expect_equal(getNumberOfCells(highDensity,      0),  100)

    expect_equal(getNumberOfCells(simpleModel,     10),   14)
    expect_equal(getNumberOfCells(longRun,        100),   20)
    expect_equal(getNumberOfCells(largeRun,         1), 1037)
    expect_equal(getNumberOfCells(highDensity,     10),  122)

    expect_equal(getCoordinates(simpleModel, 0)[1,1], -8.56, tolerance=0.01)
    expect_equal(getCoordinates(longRun,     0)[1,1], -6.24, tolerance=0.01)
    expect_equal(getCoordinates(largeRun,    0)[1,1],   -68, tolerance=0.01)
    expect_equal(getCoordinates(highDensity, 0)[1,1],  7.98, tolerance=0.01)

    expect_equal(sum(getRadii(simpleModel, 0)), 12.65, tolerance=0.01)
    expect_equal(sum(getRadii(longRun,     0)),  6.08, tolerance=0.01)
    expect_equal(sum(getRadii(largeRun,    0)),  1208, tolerance=0.01)
    expect_equal(sum(getRadii(highDensity, 0)),   123, tolerance=0.01)

    expect_equal(sum(getAxisLength(simpleModel, 0)), 25.30, tolerance=0.01)
    expect_equal(sum(getAxisLength(longRun,     0)),  12.2, tolerance=0.01)
    expect_equal(sum(getAxisLength(largeRun,    0)),  2467, tolerance=0.01)
    expect_equal(sum(getAxisLength(highDensity, 0)),   248, tolerance=0.01)

    expect_equal(sum(getAxisAngle(simpleModel, 0)), 29.5, tolerance=0.01)
    expect_equal(sum(getAxisAngle(longRun,     0)), 16.6, tolerance=0.01)
    expect_equal(sum(getAxisAngle(largeRun,    0)), 3125, tolerance=0.01)
    expect_equal(sum(getAxisAngle(highDensity, 0)),  297, tolerance=0.01)

    expect_equal(mean(getCycleLengths(simpleModel, 0)), 48)
    expect_equal(mean(getCycleLengths(longRun,     0)), 48)
    expect_equal(mean(getCycleLengths(largeRun,    0)), 48)
    expect_equal(mean(getCycleLengths(highDensity, 0)), 48)

    phases <- getCellPhases(largeRun, 0)
    expect_equal(sum(phases == 'I'), 956)
    expect_equal(sum(phases == 'M'),  44)

    expect_equal(getCellTypes(simpleModel, 0)[1], 'DEFAULT')
    expect_equal(getCellTypes(longRun,     0)[1], 'DEFAULT')
    expect_equal(getCellTypes(largeRun,    0)[1], 'DEFAULT')
    expect_equal(getCellTypes(highDensity, 0)[1], 'DEFAULT')
})

test_that('cell model calculations',
{
    data(SampleModels)


})


