#' Fit TMB model
#'
#' @export

fit_sample_tmb <- function(data, par, random) {
  obj <-  MakeADFun(data = data,
                    parameters = par,
                    DLL = "dfertility",
                    random = random,
                    hessian = FALSE)

  f <- nlminb(obj$par, obj$fn, obj$gr)
  f$par.fixed <- f$par
  f$par.full <- obj$env$last.par

  fit <- c(f, obj = list(obj))
  fit$sdreport <- sdreport(fit$obj, fit$par)

  class(fit) <- "naomi_fit"  # this is hacky...
  fit <- sample_tmb(fit, random_only=FALSE)
  return(fit)
}
