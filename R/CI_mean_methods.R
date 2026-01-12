

#' Print and coerce a NAVAE_CI_Mean object
#'
#' @param x the object
#'
#' @param verbose if zero, only basic printing is done. Higher values corresponds
#' to more detailed output.
#'
#' @param ... other arguments, currently ignored.
#'
#'
#' @returns
#' \code{print.Navae_ci_ols} prints information about \code{x} and returns it
#' invisibly.
#'
#' \code{as.data.frame} returns a \code{data.frame} with 2 rows.
#'
#'
#' @references
#' Derumigny, A., Girard, L., & Guyonvarch, Y. (2025).
#' Can we have it all?
#' Non-asymptotically valid and asymptotically exact confidence intervals
#' for expectations and linear regressions.
#' ArXiv preprint, \doi{10.48550/arXiv.2507.16776}
#'
#' @seealso
#' The function to generate such objects \code{\link{Navae_ci_mean}}.
#'
#' The corresponding methods for the regression (OLS):
#' \code{\link{print.NAVAE_CI_OLS}} and
#' \code{\link{as.data.frame.NAVAE_CI_OLS}}.
#'
#' @examples
#' n = 10000
#' x = rexp(n, 1)
#' myCI = Navae_ci_mean(x, bound_K = 9, alpha = 0.2)
#'
#' print(myCI)
#' as.data.frame(myCI)
#'
#'
#' @export
print.NAVAE_CI_Mean <- function(x, verbose = 0, ...){

  cat("Call: ")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"))
  cat("\n\n")

  cat("CLT-based confidence interval:\n")

  print(x$ci_asymp)

  cat("\n")

  cat("NAVAE confidence interval:\n")

  print(x$ci_navae)

  if (verbose >= 1){
    cat("\n")
    cat("n =", x$n, "\n")
    cat("a = ", x$a, " (", x$properties_a$method, ")\n", sep = "")
    cat("K = ", x$bound_K_value, " (", x$bound_K_method, ")\n", sep = "")
  }

  return (invisible(x))
}


#' @rdname print.NAVAE_CI_Mean
#' @export
as.data.frame.NAVAE_CI_Mean <- function(x, ...){
  # ci_asymp = as.data.frame(x$ci_asymp)
  # ci_navae = as.data.frame(x$ci_navae)

  result = rbind(x$ci_asymp, x$ci_navae)

  result = as.data.frame(result)

  result$method = c("Asymptotic (CLT-based)",
                    paste0("NAVAE (", x$delta_n_from, "-based)"))

  return (result)
}
