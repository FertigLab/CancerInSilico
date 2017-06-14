context('Testing CellModel Access')

test_that('cell model getters',
{
    data(SampleModels)  

    expect_equal(getNumberOfCells(modDefault,     0),   10)
    expect_equal(getNumberOfCells(modLongRun,     0),    5)
    expect_equal(getNumberOfCells(modLargeRun,    0), 1000)
    expect_equal(getNumberOfCells(modHighDensity, 0),  100)

    expect_equal(getNumberOfCells(modDefault, 10), 15)
    expect_equal(getCoordinates(modDefault, 0)[1,1], -6.74, tolerance=0.01)
    expect_equal(sum(getRadius(modDefault, 0)), 12.3, tolerance=0.01)
    expect_equal(sum(getAxisLength(modDefault, 0)), 26.4, tolerance=0.01)
    expect_equal(sum(getAxisAngle(modDefault, 0)), 29.5, tolerance=0.01)
    expect_equal(mean(getCycleLengths(modDefault, 0)), 24)

    phases <- getCellPhases(modLargeRun, 0)
    expect_equal(sum(phases == 'I'), 904)
    expect_equal(sum(phases == 'M'), 96)

    ct <- getCellTypes(modDefault, 0)
    expect_equal(modDefault@cellTypes[[ct[1]]]@name, 'DEFAULT')
})

test_that('cell model calculations',
{
    data(SampleModels)

    expect_equal(getDensity(modDefault, 0), 0.153, tolerance=0.01)
    expect_equal(getDensity(modDefault, 10), 0.205, tolerance=0.01)

    # TODO: test local density
})


