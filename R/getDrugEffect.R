#'\code{getDrugEffect} 
#'
#'
#' @param slwoDown A vector specifying the slow down of the growth rate
#' @param cycleLengthSeq A sequence spanning the range of cycle lengths
#' @return A list of vectors specifying the distribution of drug effects depending on the growth rate of the cell
#' @examples
#' getDrugEffect(0.3, seq(2,12,0.1))
#' @export

getDrugEffect <- function(slowDown, cycleLengthSeq) {
## allow user to pass function, raw cycleLengthDist

    ret_list <- vector("list", length(cycleLengthSeq))

    if (length(slowDown) == 1) {

        for (i in 1:length(ret_list)) {

            ret_list[[i]] = c(cycleLengthSeq[i], slowDown)

        }

    } else if (length(slowDown) != length(cycleLengthSeq)) {

        stop()

    } else {

        for (i in 1:length(ret_list)) {

            ret_list[[i]] = c(cycleLengthSeq[i], slowDown[i])

        }

    }

    return (ret_list)

}
