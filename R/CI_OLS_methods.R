

#' Print and coerce a NAVAE_CI_Regression object
#'
#' This also displays CLT-based confidence intervals. The results are different
#' from the confidence intervals that can be obtained via \code{confint(lm( ))}
#' since they are robust to heteroscedasticity.
#'
#' @param x the object
#'
#' @param verbose if zero, only basic printing is done. Higher values corresponds
#' to more detailed output.
#'
#' @param ... additional arguments, currently ignored.
#'
#' @returns
#' \code{print.Navae_ci_ols} prints information about \code{x} and returns it
#' invisibly.
#'
#' \code{as.data.frame.NAVAE_CI_Regression} returns a data frame consisting
#' of two observations for each vector u given as a line of \code{matrix_u},
#' with the following columns:
#' \itemize{
#'    \item \code{name}: name of the estimateed coefficient in the linear model
#'    \item \code{lower}: lower bound of the confidence interval
#'    \item \code{upper}: upper bound of the confidence interval
#'    \item \code{estimate}: the estimated value of the coefficient
#'    \item \code{length}: the length of the interval
#'
#'    \item \code{method}: the method used for the computation of the confidence
#'    intervals. This is either "Asymptotic (CLT-based), or "NAVAE (BE-based)",
#'    or "NAVAE (EE-based)".
#'
#'    \item \code{regime}: the regime used for the computation of the CI
#'    (only applicable for NAVAE confidence intervals).
#'    Four regimes are possible: \itemize{
#'        \item the degenerate regimes \code{R1} and \code{R2} in which
#'        the confidence interval is \code{(-Inf, Inf)}.
#'        \item the exponential regime \code{Exp}
#'        \item the Edgeworth regime \code{Edg}.
#'    }
#' }
#'
#' @examples
#' n = 4000
#' X1 = rnorm(n, sd = 1)
#' true_eps = rnorm(n)
#' Y = 8 * X1 + true_eps
#' X = cbind(X1)
#'
#' myCI <- Navae_ci_ols(Y, X, K_xi = 3, intercept = TRUE, a = 1.1)
#'
#' print(myCI)
#' as.data.frame(myCI)
#'
#'
#' @export
print.NAVAE_CI_Regression <- function(x, verbose = 0, ...){
  cat("Call: ")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"))
  cat("\n\n")

  cat("CLT-based confidence intervals:\n")

  print(x$ci_asymp)

  cat("\n")

  cat("NAVAE confidence intervals:\n")

  print(x$ci_navae)

  if (verbose >= 1){
    cat("\n")
    print(x$allTuningParameters)
  }

  if (verbose >= 1){
    cat("Bounds: \n")
    print(x$allBounds, row.names = FALSE)
  }

  return (invisible(x))
}


#' @rdname print.NAVAE_CI_Regression
#' @export
as.data.frame.NAVAE_CI_Regression <- function(x, ...){
  ci_asymp = x$ci_asymp
  ci_navae = x$ci_navae

  ci_asymp$method = "Asymptotic (CLT-based)"
  ci_navae$method = paste0("NAVAE (", x$about_delta_n$delta_n_from_u, "-based)")

  # Needs to define this so that ci_asymp and ci_navae have the same columns.
  ci_asymp$regime = NA

  result = rbind(ci_asymp, ci_navae)
  rownames(result) <- NULL

  result = cbind(name = rownames(ci_asymp), result)

  return (result)
}
