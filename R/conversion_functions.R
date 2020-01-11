#' Convert simulated output to MRexperiment object
#'
#' In order to allow investigators to more easily incorporate simulated data,
#' this package converts the raw output into an MRexperiment object used in the
#' \code{\link{metagenomeSeq}} package.
#'
#' @param obj
#'
#' @examples
#' gen_norm_microbiome()
#'
#@export
#'
#' @return
#' An MRexperiment object
#'
#' @importFrom metagenomeSeq newMRexperiment
simulate2MRexperiment <- function(obj){

}

#' Convert simulated output to phyloseq object
#'
#' This function will convert simulated data into a \code{\link{phyloseq}}
#' object.
#'
#' @param obj
