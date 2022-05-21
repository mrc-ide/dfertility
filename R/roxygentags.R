#' @title Export TMB models
#'
#' This is a dummy function who's purpose is to hold the useDynLib roxygen tag.
#' This tag will populate the namespace with compiled c++ functions upon
#' package install.
#'
#' @useDynLib model7_eth
#' @useDynLib model7_mwi
#' @useDynLib model7_rwa
#' @useDynLib model7_zwe
#' @useDynLib model7_sd_caf
#'
#' @noRd
#' @keywords internal
.dummy <- function() {
  return(NULL)
}
