#' Fit TMB model
#'
#' @export

fit_sample_tmb <- function(data, par, random) {
  obj <-  TMB::MakeADFun(data = data,
                    parameters = par,
                    DLL = "dfertility",
                    random = random,
                    hessian = FALSE)

  f <- stats::nlminb(obj$par, obj$fn, obj$gr)
  f$par.fixed <- f$par
  f$par.full <- obj$env$last.par

  fit <- c(f, obj = list(obj))
  # fit$sdreport <- sdreport(fit$obj, fit$par)

  class(fit) <- "naomi_fit"  # this is hacky...
  fit <-  naomi::sample_tmb(fit, random_only = TRUE)
  return(fit)
}

#' Sample TMB
sample_tmb <- function (fit, nsample = 1000, rng_seed = NULL, random_only = TRUE,
          verbose = FALSE)
{
  set.seed(rng_seed)
  stopifnot(methods::is(fit, "naomi_fit"))
  stopifnot(nsample > 1)
  if (!random_only) {
    if (verbose)
      print("Calculating joint precision")
    hess <- sdreport_joint_precision(fit$obj, fit$par.fixed)
    if (verbose)
      print("Inverting precision for joint covariance")
    cov <- solve(hess)
    if (verbose)
      print("Drawing sample")
    smp <- mvtnorm::rmvnorm(nsample, fit$par.full, cov)
  }
  else {
    r <- fit$obj$env$random
    par_f <- fit$par.full[-r]
    par_r <- fit$par.full[r]
    hess_r <- fit$obj$env$spHess(fit$par.full, random = TRUE)
    smp_r <-  rmvnorm_sparseprec(nsample, par_r, hess_r)
    smp <- matrix(0, nsample, length(fit$par.full))
    smp[, r] <- smp_r
    smp[, -r] <- matrix(par_f, nsample, length(par_f), byrow = TRUE)
    colnames(smp)[r] <- colnames(smp_r)
    colnames(smp)[-r] <- names(par_f)
  }
  if (verbose)
    print("Simulating outputs")
  sim <- apply(smp, 1, fit$obj$report)
  r <- fit$obj$report()
  if (verbose)
    print("Returning sample")
  fit$sample <- Map(vapply, list(sim), "[[", lapply(lengths(r),
                                                    numeric), names(r))
  is_vector <- vapply(fit$sample, inherits, logical(1), "numeric")
  fit$sample[is_vector] <- lapply(fit$sample[is_vector], matrix,
                                  nrow = 1)
  names(fit$sample) <- names(r)
  fit
}
