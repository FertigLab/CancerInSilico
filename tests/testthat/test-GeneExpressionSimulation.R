context('Testing Gene Expression Simulation')

test_that('Check Data Set',
{
    data(ReferenceGeneExpression)  
    data(SamplePathways)
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

test_that('Simulate Gene Expression - microarray',
{
    pwyGrowth <- calibratePathway(pwyGrowth, referenceGeneExpression)
    pwyMitosis <- calibratePathway(pwyMitosis, referenceGeneExpression)

    params <- new('GeneExpressionParams')
    params@RNAseq <- FALSE
    params@singleCell <- FALSE
    params@nCells <- 10

    #gs <- inSilicoGeneExpression(modDefault, c(pwyMitosis, pwyGrowth), params)

    #expect_equal(nrow(gs$expression), 546)
    #expect_equal(ncol(gs$expression), 11)
})

test_that('Simulate Gene Expression - RNA-seq, single cell',
{

})

