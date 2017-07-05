context('Testing CellModel Access')

test_that('cell model getters',
{
    data(SampleModels)  

    expect_equal(getNumberOfCells(modDefault,     0),   10)
    expect_equal(getNumberOfCells(modLongRun,     0),    5)
    expect_equal(getNumberOfCells(modLargeRun,    0), 1000)
    expect_equal(getNumberOfCells(modHighDensity, 0),  100)

    expect_equal(getNumberOfCells(modDefault, 10), 15)
    expect_equal(getCoordinates(modDefault, 0, 1)[1], -6.74, tolerance=0.01)
    expect_equal(getRadius(modDefault, 0, 1), 1.35, tolerance=0.01)
    expect_equal(getAxisLength(modDefault, 0, 1), 2.7, tolerance=0.01)
    expect_equal(getAxisAngle(modDefault, 0, 1), 3.72, tolerance=0.01)
    expect_equal(getCycleLength(modDefault, 0, 1), 24)

    phases <- sapply(1:1000, getCellPhase, model=modLargeRun, time=0)
    expect_equal(sum(phases == 'I'), 904)
    expect_equal(sum(phases == 'M'), 96)

    ct <- getCellType(modDefault, 0, 1)
    expect_equal(modDefault@cellTypes[[ct]]@name, 'DEFAULT')
})

test_that('cell model calculations',
{
    data(SampleModels)

    expect_equal(getDensity(modDefault, 0), 0.153, tolerance=0.01)
    expect_equal(getDensity(modDefault, 10), 0.205, tolerance=0.01)
})

test_that('local density',
{
    mod <- inSilicoCellModel(2, 1, 0.1, recordIncrement=0.5)
    mod@cells[[1]][1:5] <- c(0, 0, 1, 2, 0)
    mod@cells[[1]][10:14] <- c(2, 2, 1, 2, 0)
    expect_equal(getLocalDensity(mod, 0, 1, 2), 0.01, tol=0.01)
    expect_equal(getLocalDensity(modDefault, 5, 4, 4), 0.34, tol=0.01)
})

