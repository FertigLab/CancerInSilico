context("Testing functions in GeneExpressionSimulation.R")

test_that("combineGeneExpression", {

    # create several gene expression matrices
    gs1 <- replicate(25, runif(5,2,8)) 
    gs2 <- replicate(25, runif(10,10,16))
    gs3 <- replicate(25, runif(8,18,24))

    # give them overlapping 'gene' sets    
    row.names(gs1) <- letters[1:5]
    row.names(gs2) <- letters[3:12]
    row.names(gs3) <- letters[11:18]

    # create list of matrices
    gs_list <- list()
    gs_list[[1]] <- gs1
    gs_list[[2]] <- gs2
    gs_list[[3]] <- gs3

    # combine matrices    
    expect_error(combineGeneExpression(gs_list), regexp = NA)
    total_gs <- combineGeneExpression(gs_list)
    
    # check properties of final matrix
    expect_equal(nrow(total_gs), 18)
    expect_equal(ncol(total_gs), 25)

    expect_true(min(total_gs[1:2,]) > 2)
    expect_true(max(total_gs[1:2,]) < 8)

    expect_true(min(total_gs[3:10,]) > 10)
    expect_true(max(total_gs[3:10,]) < 16)

    expect_true(min(total_gs[11:18,]) > 18)
    expect_true(max(total_gs[11:18,]) < 24)

    # throw error if columns uneven
    gs_list[[4]] <- replicate(22, runif(3,0,1))
    expect_error(combineGeneExpression(gs_list))

})

test_that("simulateExpression", {

    



})
