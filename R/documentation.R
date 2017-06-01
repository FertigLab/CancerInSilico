#' CancerInSilico
#'
#' \tabular{ll}{
#' Package: \tab CancerInSilico\cr
#' Type: \tab Package\cr
#' Version: \tab 1.99.0\cr
#' Date: \tab 2017-06-24\cr
#' License: \tab LGPL\cr
#' }
#' @author Maintainer: Elana J. Fertig \email{ejfertig@jhmi.edu},
#'             Thomas D. Sherman \email{tsherma4@jhu.edu}
#' @docType package
#' @name CancerInSilico
#' @useDynLib CancerInSilico
NULL

#' This data contains the list of gene targets of modeled signaling pathways
#' 
#' @docType data
#' @name inSilicoPathways
#' @usage inSilicoPathways
#' @details List of gene targets for GtoM, GtoS, and Growth factor pathways modeled in CancerInSilico. GtoM represents genes associated with the G to M checkpoint of the cell cycle collected from PID pathways PI3KCI_AKT, TNF, TGFBR, and RB1 in MSigDB v5.1. GtoS represented genes associated with the G to S checkpoint of the cell cycle collected from E2F target genes in PID and hallmark genesets from MSigDb v5.1 and exprimentally validated targets in TRANSFAC professional 2014. Growth represents the targets of transcription factors resulting from EGFR activition (STAT1, STAT3, STAT5, Elk-1, Myc, NF-kappaB, and RelA) from the TRANSFAC professional database 2014 as described in Fertig et al, 2012.
#' @references Subramanian, Tamayo, et al, PNAS, 102:155545, 2005; Matys et al, Nuc Ac Res, 34:D108, 2006; and Fertig et al, BMC Genomics, 13:160, 2012.
#' @return List of gene targets for GtoM, GtoS, and Growth factor pathways.
#' 
NULL
