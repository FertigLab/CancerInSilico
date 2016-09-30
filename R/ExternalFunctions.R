#'\code{getDrugEffect} 
#'
#' @param cycleLengthSeq A sequence spanning the range of cycle lengths
#' @return A list of vectors specifying the distribution of drug effects depending on the growth rate of the cell
#' @examples
#' getDrugEffect(0.3, seq(2,12,0.1))
#' @export

getDrugEffect <- function(FUN = function(x) {0}, ...) {

    # store either the cycle lengths themselves or a sequence spanning their range
    cycleLengthDist <- list(...)$cycleLengthDist
    cycleLengthSeq <- list(...)$cycleLengthSeq

    # check which parameters are provided
    if (is.null(cycleLengthDist)) {

        if (is.null(cycleLengthSeq)) {

            # need at least one
            stop("must specifiy either cycleLengthDist or cycleLengthSeq")

        }

    # if only cycle lengths are provided ...
    } else if (is.null(cycleLengthSeq)) {

        # create a sequence that spans their range
        cycleLengthSeq <- seq(min(cycleLengthDist), max(cycleLengthDist), 0.1)

    } 

    # create return list    
    ret_list <- vector("list", length(cycleLengthSeq))

    # for each element in the sequence
    for (i in 1:length(ret_list)) {

        # apply the function to the growth rate, note that FUN could return
        # vector of values which the model then uses to sample from
        ret_list[[i]] = c(cycleLengthSeq[i], FUN(cycleLengthSeq[i]))

    }

    # return the list
    return (ret_list)

}
