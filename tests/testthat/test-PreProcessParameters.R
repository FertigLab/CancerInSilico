context('Testing CellModel Creation')

test_that('throw error when missing needed parameters',
{  
    cellModel <- new('CellModel')
    offLattice  <- new('OffLatticeModel')
    drasdoHohme <- new('DrasdoHohmeModel')

    # CellModel
    expect_error(processParameters(cellModel), 'missing initialNum')
    cellModel@params[['initialNum']] <- 100
    expect_error(processParameters(cellModel), 'missing runTime')
    cellModel@params[['runTime']] <- 48
    expect_error(processParameters(cellModel), 'missing density')
    cellModel@params[['density']] <- 0.01
    expect_error(processParameters(cellModel), 'missing boundary')
    cellModel@params[['boundary']] <- TRUE
    expect_error(processParameters(cellModel), 'missing syncCycles')
    cellModel@params[['syncCycles']] <- TRUE
    expect_error(processParameters(cellModel), NA)

    # OffLatticeModel
    offLattice@params <- cellModel@params
    expect_error(processParameters(offLattice), 'missing maxDeform')
    offLattice@params[['maxDeform']] <- 1.0
    expect_error(processParameters(offLattice), 'missing maxTranslation')
    offLattice@params[['maxTranslation']] <- 1.0
    expect_error(processParameters(offLattice), 'missing maxRotate')
    offLattice@params[['maxRotate']] <- 1.0
    expect_error(processParameters(offLattice), NA)

    # DrasdoHohmeModel
    drasdoHohme@params <- cellModel@params
    expect_error(processParameters(drasdoHohme), 'missing nG')
    drasdoHohme@params[['nG']] <- 28
    expect_error(processParameters(drasdoHohme), 'missing epsilon')
    drasdoHohme@params[['epsilon']] <- 10.0
    expect_error(processParameters(drasdoHohme), 'missing delta')
    drasdoHohme@params[['delta']] <- 0.5
    expect_error(processParameters(drasdoHohme), NA)
})

test_that('assign parameter default values',
{
    cellModel <- new('CellModel')
    offLattice  <- new('OffLatticeModel')
    drasdoHohme <- new('DrasdoHohmeModel')

    # supply minimum neccesary parameters
    cellModel@params[['initialNum']] <- 100
    cellModel@params[['runTime']] <- 48
    cellModel@params[['density']] <- 0.01
    cellModel@params[['boundary']] <- TRUE
    cellModel@params[['syncCycles']] <- TRUE

    offLattice@params <- cellModel@params
    offLattice@params[['maxDeform']] <- 1.0
    offLattice@params[['maxTranslation']] <- 1.0
    offLattice@params[['maxRotate']] <- 1.0

    drasdoHohme@params <- cellModel@params
    drasdoHohme@params[['nG']] <- 28
    drasdoHohme@params[['epsilon']] <- 10.0
    drasdoHohme@params[['delta']] <- 0.5

    # get default parameters
    cellModel@params <- processParameters(cellModel)
    offLattice@params <- processParameters(offLattice)  
    drasdoHohme@params <- processParameters(drasdoHohme)

    # CellModel
    expect_equal(cellModel@params[['randSeed']], 0)
    expect_equal(cellModel@params[['outputIncrement']], 6)
    expect_equal(cellModel@params[['recordIncrement']], 0.1)
    expect_equal(cellModel@params[['timeIncrement']], 0.01)
    expect_true(is.null(cellModel@params[['drugs']]))
    expect_equal(length(cellModel@params[['cellTypes']]), 1)
    expect_equal(length(cellModel@params[['cellTypesInitFreq']]), 1)
    expect_equal(cellModel@params[['cellTypesInitFreq']][[1]], 1)
    type <- cellModel@params[['cellTypes']][[1]]
    expect_equal(type@name, 'DEFAULT')
    expect_equal(type@size, 1)
    expect_equal(type@minCycle, 48)
    expect_equal(type@cycleLength(), 48)

    # OffLatticeModel
    expect_equal(offLattice@params[['randSeed']], 0)
    expect_equal(offLattice@params[['outputIncrement']], 6)
    expect_equal(offLattice@params[['recordIncrement']], 0.1)
    expect_equal(offLattice@params[['timeIncrement']], 0.01)
    expect_true(is.null(offLattice@params[['drugs']]))
    expect_equal(length(offLattice@params[['cellTypes']]), 1)
    expect_equal(length(offLattice@params[['cellTypesInitFreq']]), 1)
    expect_equal(offLattice@params[['cellTypesInitFreq']][[1]], 1)
    type <- offLattice@params[['cellTypes']][[1]]
    expect_equal(type@name, 'DEFAULT')
    expect_equal(type@size, 1)
    expect_equal(type@minCycle, 48)
    expect_equal(type@cycleLength(), 48)

    # DrasdoHohmeModel
    expect_equal(drasdoHohme@params[['randSeed']], 0)
    expect_equal(drasdoHohme@params[['outputIncrement']], 6)
    expect_equal(drasdoHohme@params[['recordIncrement']], 0.1)
    expect_equal(drasdoHohme@params[['timeIncrement']], 0.00173,
        tolerance=1e-5)
    expect_true(is.null(drasdoHohme@params[['drugs']]))
    expect_equal(length(drasdoHohme@params[['cellTypes']]), 1)
    expect_equal(length(drasdoHohme@params[['cellTypesInitFreq']]), 1)
    expect_equal(drasdoHohme@params[['cellTypesInitFreq']][[1]], 1)
    type <- drasdoHohme@params[['cellTypes']][[1]]
    expect_equal(type@name, 'DEFAULT')
    expect_equal(type@size, 1)
    expect_equal(type@minCycle, 48)
    expect_equal(type@cycleLength(), 48)  
})

test_that('custom parameters - CellModel',
{
    # needed parameters
    cellModel <- new('CellModel')
    cellModel@params[['initialNum']] <- 100
    cellModel@params[['runTime']] <- 48
    cellModel@params[['density']] <- 0.01
    cellModel@params[['boundary']] <- TRUE
    cellModel@params[['syncCycles']] <- TRUE

    # initialNum
    cellModel@params[['initialNum']] <- -2    
    expect_error(processParameters(cellModel), 'invalid initialNum')
    cellModel@params[['initialNum']] <- 100

    # runTime
    cellModel@params[['runTime']] <- -100    
    expect_error(processParameters(cellModel), 'invalid runTime')
    cellModel@params[['runTime']] <- 48

    # density
    cellModel@params[['density']] <- 0    
    expect_error(processParameters(cellModel), 'invalid density')
    cellModel@params[['density']] <- 0.8
    expect_error(processParameters(cellModel), paste('density too high', 
        'for seeding'))
    cellModel@params[['density']] <- 0.01

    # outputIncrement
    cellModel@params[['outputIncrement']] <- 0
    expect_error(processParameters(cellModel), 'invalid outputIncrement')
    cellModel@params[['outputIncrement']] <- 6

    # recordIncrement
    cellModel@params[['recordIncrement']] <- -3
    expect_error(processParameters(cellModel), 'invalid recordIncrement')
    cellModel@params[['recordIncrement']] <- 0.1

    # timeIncrement
    cellModel@params[['timeIncrement']] <- 0
    expect_error(processParameters(cellModel), 'invalid timeIncrement')
    cellModel@params[['timeIncrement']] <- 0.01

    # cellTypes
    type1 <- newCellType('TYPE_1', 1, minCycle = 20, cycleLength = 
        function() {return (runif(1, 20, 50))})
    type2 <- newCellType('TYPE_1', 1, minCycle = 30, cycleLength = 
        function() {return (30 + rexp(1, 1))})

    cellModel@params[['cellTypes']] <- c(type1, type2)
    cellModel@params[['cellTypesInitFreq']] <- c(0.3, 0.7)
    expect_error(processParameters(cellModel), NA)
    cellModel@params[['cellTypesInitFreq']] <- c(0.3, 0.8)
    expect_error(processParameters(cellModel), paste('cell type frequency',
        'doesn\'t sum to 1'))
    cellModel@params[['cellTypesInitFreq']] <- c(0.3, 0.4, 0.3)
    expect_error(processParameters(cellModel), paste('cell type frequency',
        'size != cell type size'))

    cellModel@params[['cellTypes']] <- c(type1, type2, 3.0)
    expect_error(processParameters(cellModel), paste('not all cell types',
        'are valid'))
    cellModel@params[['cellTypes']] <- c(type1, type2)
    cellModel@params[['cellTypesInitFreq']] <- c(0.3, 0.7)

    # drugs
    drug1 <- newDrug('DRUG_1', 0.0, function(t,l,p) {return(l / 2)})
    drug2 <- newDrug('DRUG_2', 24.0)

    cellModel@params[['drugs']] <- c(drug1, drug2)
    expect_error(processParameters(cellModel), NA)

    cellModel@params[['drugs']] <- c(drug1, drug2, 'ABCD')
    expect_error(processParameters(cellModel), 'not all drugs are valid')
})

test_that('custom parameters - OffLatticeModel',
{
    offLattice <- new('OffLatticeModel')
    offLattice@params[['initialNum']] <- 100
    offLattice@params[['runTime']] <- 48
    offLattice@params[['density']] <- 0.01
    offLattice@params[['boundary']] <- TRUE
    offLattice@params[['syncCycles']] <- TRUE
    offLattice@params[['maxDeform']] <- 1.0
    offLattice@params[['maxTranslation']] <- 1.0
    offLattice@params[['maxRotate']] <- 1.0

    offLattice@params[['maxDeform']] <- 0
    expect_error(processParameters(offLattice), 'invalid maxDeform')
    offLattice@params[['maxDeform']] <- 1.0

    offLattice@params[['maxTranslation']] <- -1
    expect_error(processParameters(offLattice), 'invalid maxTranslation')
    offLattice@params[['maxTranslation']] <- 1.0

    offLattice@params[['maxRotate']] <- -5
    expect_error(processParameters(offLattice), 'invalid maxRotate')
})

test_that('custom parameters - DrasdoHohmeModel',
{
    drasdoHohme <- new('OffLatticeModel')
    drasdoHohme@params[['initialNum']] <- 100
    drasdoHohme@params[['runTime']] <- 48
    drasdoHohme@params[['density']] <- 0.01
    drasdoHohme@params[['boundary']] <- TRUE
    drasdoHohme@params[['syncCycles']] <- TRUE
    drasdoHohme@params[['nG']] <- 28
    drasdoHohme@params[['epsilon']] <- 10.0
    drasdoHohme@params[['delta']] <- 0.2

    drasdoHohme@params[['
    expect_error(processParameters(offLattice), 'invalid maxDeform')


})

test_that('ellipsis parameters',
{

})

test_that('CellModel constructor',
{

})
