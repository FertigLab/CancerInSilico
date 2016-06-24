#'\code{plotGrowthRateDistribution} Plots a distribution of the growth rates
#'of all the cells
#'
#'
#'@param 

setGeneric("plotGrowthRateDistribution", function(mat, time)
  standardGeneric("plotGrowthRateDistribution"))

setMethod("plotGrowthRateDistribution", "CellMatrix",
          
  function(mat, time) {

    radii = seq(3,ncol(mat),6)
    gr_rates = radii[mat[time,radii]>0] + 3
    plot(density(mat[time,gr_rates]))

  }
          
)