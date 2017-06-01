context('Testing CellModel Access')

test_that('cell model getters',
{
    data(SampleModels)  

    expect_equal(getNumberOfCells(modDefault,     0),   10)
    expect_equal(getNumberOfCells(modLongRun,     0),    5)
    expect_equal(getNumberOfCells(modLargeRun,    0), 1000)
    expect_equal(getNumberOfCells(modHighDensity, 0),  100)

    expect_equal(getNumberOfCells(modDefault,     10),   14)
    expect_equal(getNumberOfCells(modLongRun,    100),   20)
    expect_equal(getNumberOfCells(modLargeRun,     1), 1037)
    expect_equal(getNumberOfCells(modHighDensity, 10),  122)

    expect_equal(getCoordinates(modDefault,    0)[1,1],-8.56,tolerance=0.01)
    expect_equal(getCoordinates(modLongRun,    0)[1,1],-6.24,tolerance=0.01)
    expect_equal(getCoordinates(modLargeRun,   0)[1,1],  -68,tolerance=0.01)
    expect_equal(getCoordinates(modHighDensity,0)[1,1], 7.98,tolerance=0.01)

    expect_equal(sum(getRadii(modDefault,     0)), 12.65, tolerance=0.01)
    expect_equal(sum(getRadii(modLongRun,     0)),  6.08, tolerance=0.01)
    expect_equal(sum(getRadii(modLargeRun,    0)),  1208, tolerance=0.01)
    expect_equal(sum(getRadii(modHighDensity, 0)),   123, tolerance=0.01)

    expect_equal(sum(getAxisLength(modDefault,    0)),25.30, tolerance=0.01)
    expect_equal(sum(getAxisLength(modLongRun,    0)), 12.2, tolerance=0.01)
    expect_equal(sum(getAxisLength(modLargeRun,   0)), 2467, tolerance=0.01)
    expect_equal(sum(getAxisLength(modHighDensity,0)),  248, tolerance=0.01)

    expect_equal(sum(getAxisAngle(modDefault,     0)), 29.5, tolerance=0.01)
    expect_equal(sum(getAxisAngle(modLongRun,     0)), 16.6, tolerance=0.01)
    expect_equal(sum(getAxisAngle(modLargeRun,    0)), 3125, tolerance=0.01)
    expect_equal(sum(getAxisAngle(modHighDensity, 0)),  297, tolerance=0.01)

    expect_equal(mean(getCycleLengths(modDefault,     0)), 48)
    expect_equal(mean(getCycleLengths(modLongRun,     0)), 48)
    expect_equal(mean(getCycleLengths(modLargeRun,    0)), 48)
    expect_equal(mean(getCycleLengths(modHighDensity, 0)), 48)

    phases <- getCellPhases(modLargeRun, 0)
    expect_equal(sum(phases == 'I'), 956)
    expect_equal(sum(phases == 'M'),  44)

    ct <- getCellTypes(modDefault,0)
    expect_equal(modDefault@cellTypes[[ct[1]]]@name, 'DEFAULT')
})

test_that('cell model calculations',
{
    data(SampleModels)

    expect_equal(getDensity(modDefault,     0), 0.1)
    expect_equal(getDensity(modLongRun,     0), 0.1)
    expect_equal(getDensity(modLargeRun,    0), 0.1)
    expect_equal(getDensity(modHighDensity, 0), 0.4)

    expect_equal(getDensity(modDefault,      10), 0.11, tolerance=0.01)
    expect_equal(getDensity(modLongRun,     100), 0.39, tolerance=0.01)
    expect_equal(getDensity(modLargeRun,      1), 0.10, tolerance=0.01)
    expect_equal(getDensity(modHighDensity,  10), 0.46, tolerance=0.01)

    expect_equal(getNumberOfNeighbors(modDefault, 0, 1, 0), 0)
    expect_equal(getNumberOfNeighbors(modDefault, 0, 1,
        2 * modDefault@boundary), 9)

    expect_equal(getNumberOfNeighbors(modLongRun, 0, 1, 0), 0)
    expect_equal(getNumberOfNeighbors(modLongRun, 0, 1,
        2 * modLongRun@boundary), 4)

    expect_equal(getNumberOfNeighbors(modLargeRun, 0, 1, 0), 0)
    expect_equal(getNumberOfNeighbors(modLargeRun, 0, 1,
        2 * modLargeRun@boundary), 999)

    expect_equal(getNumberOfNeighbors(modHighDensity, 0, 1, 0), 0)
    expect_equal(getNumberOfNeighbors(modHighDensity, 0, 1,
        2 * modHighDensity@boundary), 99)
})


