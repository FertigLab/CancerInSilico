context('Testing Gene Expression Simulation')

test_that('Check Data Set',
{  
    geneNames <- unique(c(pwyGrowth@genes, pwyContactInhibition@genes,
        pwyMitosis@genes, pwySPhase@genes))
    expect_warning(checkDataSet(referenceGeneExpression, geneNames), NA)
})

test_that('Combine Gene Expression',
{  
    mat1 <- matrix(runif(40,0,1), nrow=10, ncol=4)
    mat2 <- matrix(runif(32,2,3), nrow=8, ncol=4)
    mat3 <- matrix(runif(60,4,5), nrow=15, ncol=4)
    row.names(mat1) <- letters[1:10]
    row.names(mat2) <- letters[8:15]
    row.names(mat3) <- letters[12:26]
    colnames(mat1) <- colnames(mat2) <- colnames(mat3) <- LETTERS[1:4]

    mat <- combineGeneExpression(list(mat1, mat2, mat3))
    
    expect_equal(nrow(mat), 26)
    expect_equal(ncol(mat), 4)
    expect_true(all(mat[1:7,] < 1 & mat[1:7,] > 0))
    expect_true(all(mat[8:11,] < 3 & mat[8:11,] > 2))
    expect_true(all(mat[12:26,] < 5 & mat[12:26,] > 4))
})

test_that('Pad Expression Matrix',
{  
    pwyOutput <- lapply(list(pwyGrowth, pwyMitosis), function(p)
    {
        p <- calibratePathway(p, referenceGeneExpression)
        activity <- simulatePathwayActivity(p, modDefault, 1, 10)
        return(simulatePathwayExpression(p, activity))
    })

    mat1 <- combineGeneExpression(pwyOutput)
    expect_equal(nrow(mat1), 546)
    expect_equal(ncol(mat1), 11)

    dist <- function() runif(1,2,14)
    mat2 <- padExpressionMatrix(mat1, 54, dist)
    expect_equal(nrow(mat2), 600)
    expect_equal(ncol(mat2), ncol(mat1))
    expect_true(all(!is.na(mat2)))
})

test_that('Simulate Error - microarray',
{  
    pwyOutput <- lapply(list(pwyGrowth, pwyMitosis), function(p)
    {
        p <- calibratePathway(p, referenceGeneExpression)
        activity <- simulatePathwayActivity(p, modDefault, 1, 10)
        return(simulatePathwayExpression(p, activity))
    })
    mat <- combineGeneExpression(pwyOutput)

    matError <- simulateError(mat, perError=0.1, microArray=TRUE)
    expect_equal(nrow(matError), 546)
    expect_equal(ncol(matError), 11)
})

test_that('Simulate Error - RNA-seq (Neg Bin)',
{  
    pwyOutput <- lapply(list(pwyGrowth, pwyMitosis), function(p)
    {
        p <- calibratePathway(p, referenceGeneExpression)
        activity <- simulatePathwayActivity(p, modDefault, 1, 10)
        return(simulatePathwayExpression(p, activity))
    })
    mat <- combineGeneExpression(pwyOutput)

    matError1 <- simulateError(mat, microArray=FALSE)
    expect_equal(nrow(matError1), 546)
    print(matError1[1:10,])
#    expect_equal(ncol(matError1), 11)

    matError2 <- simulateError(mat, dataSet=referenceGeneExpression, 
        microArray=FALSE)
    expect_equal(nrow(matError2), 546)
#    expect_equal(ncol(matError2), 11)

})

test_that('Simulate Gene Expression - microarray',
{
    pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
    pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
    gs <- inSilicoGeneExpression(modDefault, c(pwyMitosis, pwyGrowth),
    microArray=TRUE, nCells=10)
    expect_equal(nrow(gs$expression), 546)
    expect_equal(ncol(gs$expression), 11)
})

test_that('Simulate Gene Expression - single cell',
{
    pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
    pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
    gs <- inSilicoGeneExpression(modDefault, c(pwyMitosis, pwyGrowth),
    microArray=FALSE, nCells=10)
})


