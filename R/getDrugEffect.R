#'\code{getDrugEffect} 
#'
#'
#' @param slwoDown A vector specifying the slow down of the growth rate
#' @param cycleLengthSeq A sequence spanning the range of cycle lengths
#' @return A list of vectors specifying the distribution of drug effects depending on the growth rate of the cell
#' @examples
#' getDrugEffect(0.3, seq(2,12,0.1))
#' @export

getDrugEffect <- function(FUN = function(x) {0}, ...) {

  cycleLengthDist <- list(...)$cycleLengthDist
  cycleLengthSeq <- list(...)$cycleLengthSeq
  
  if (is.null(cycleLengthDist)) {
    
    if (is.null(cycleLengthSeq)) {
    
      stop("must specifiy either cycleLengthDist or cycleLengthSeq")
    
    }

  } else if (is.null(cycleLengthSeq)) {
      
      cycleLengthSeq <- seq(min(cycleLengthDist), max(cycleLengthDist), 0.1)

  } 

  ret_list <- vector("list", length(cycleLengthSeq))

  for (i in 1:length(ret_list)) {

    ret_list[[i]] = c(cycleLengthSeq[i], FUN(cycleLengthSeq[i]))

  }
  
  return (ret_list)

}
