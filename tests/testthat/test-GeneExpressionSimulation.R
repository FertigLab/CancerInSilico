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
    pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
    pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
    pwyOutput <- lapply(list(pwyGrowth, pwyMitosis), function(p)
        simulatePathwayExpression(p, modDefault, 1))
    mat <- combineGeneExpression(pwyOutput)

    dist <- function() runif(1,2,14)
    mat1 <- padExpressionMatrix(mat, 546, dist)
    expect_equal(mat, mat1)

    mat2 <- padExpressionMatrix(mat, 600, dist)
    expect_equal(nrow(mat2), 600)
    expect_equal(ncol(mat2), ncol(mat))
    expect_true(all(!is.na(mat2)))
})

test_that('Simulate Error - microarray',
{  
    pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
    pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
    pwyOutput <- lapply(list(pwyGrowth, pwyMitosis), function(p)
        simulatePathwayExpression(p, modDefault, 1))
    mat <- combineGeneExpression(pwyOutput)

    matError <- simulateError(mat, perError=0.1, microArray=TRUE)
    
})

test_that('Simulate Error - RNA-seq (Neg Bin)',
{  
    pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
    pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
    pwyOutput <- lapply(list(pwyGrowth, pwyMitosis), function(p)
        simulatePathwayExpression(p, modDefault, 1))
    mat <- combineGeneExpression(pwyOutput)

    matError1 <- simulateError(mat, microArray=FALSE)
#    matError2 <- simulateError(mat, dataSet=referenceGeneExpression, microArray=FALSE)
})

test_that('Simulate Gene Expression',
{
    pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
    pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)
    gs <- inSilicoGeneExpression(modDefault, c(pwyMitosis, pwyGrowth),
    microArray=TRUE)
})

