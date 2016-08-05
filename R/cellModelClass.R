#' An S4 class to represent the output of a matrix
#'
#'@export

setClass("CellModel", representation(cells = "list", simple_parameters = "numeric", cycleLengthDist = "numeric", drugEffect = "list")
