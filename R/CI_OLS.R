#' Compute valid CI for an ordinary least squares regression
#'
#' @param Y vector of observations of the explained variables
#' @param X matrix of explanatory variables.
#' \code{X} must not contain a constant column.
#' @param alpha 1 - level of confidence of the CI
#' @param omega the tuning parameter \eqn{\omega} of the interval
#' @param a the tuning parameter \eqn{a} of the interval
#' @param bounds list of bounds for the DGP.
#' It can contain the following items: \itemize{
#'    \item \code{lambda_m}
#'    \item \code{K_eps}
#'    \item \code{K_xi}
#'    \item \code{K3_xi}
#'    \item \code{lambda3_xi}
#'    \item \code{K3tilde_xi}
#'    \item \code{B}, \code{C} Bounds for the concentration of || Xi tilde %*% Xi tilde'||
#'    \item \code{K_reg} Bound on
#'    \eqn{ E[ || vec( \widetilde{X}\widetilde{X}'- \mathbb{I}_p ) ||^2 ] }
#'    Defined in Assumption 3.2 (ii).
#' }
#' The bounds that are not given are replaced by plug-ins.
#' For K3_xi, lambda3_xi and K3tilde_xi, the bounds are obtained
#' from K_xi (= K4_xi)
#'
#' @param matrix_u each row of this matrix is understood as a new vector u
#' for which a confidence interval should be computed.
#' By default \code{matrix_u} is the identity matrix, corresponding
#' to the canonical basis of \eqn{R^p}.
#'
#' @return This function returns a data frame consisting of an observation
#' for each vector u given as a line of \code{matrix_u},
#' with the following columns: \itemize{
#'    \item \code{lower}: lower bound of the confidence interval
#'    \item \code{upper}: upper bound of the confidence interval
#'    \item \code{regime}: the regime used for the computation of the CI.
#'    Four regimes are possible: \itemize{
#'        \item the degenerate regimes \code{R1} and \code{R2} in which
#'        the confidence interval is \code{(-Inf, Inf)}.
#'        \item the exponential regime \code{Exp}
#'        \item the Edgeworth regime \code{Edg}.
#'    }
#' }
#'
#' @examples
#' n = 20000
#' X1 = rnorm(n, sd = 4)
#' true_eps = rnorm(n)
#' Y = 3 + 8 * X1 + true_eps
#' X = cbind(X1)
#'
#' myCI <- CI.OLS(Y, X, alpha = 0.05, omega = 0.2, a = 2,
#'   bounds = list(lambda_m = 1, K_reg = 5, K_eps = 5, K_xi = 10),
#'   setup = list(continuity = FALSE, no_skewness = FALSE) )
#'
#' print(myCI)
#'
#' myCI <- CI.OLS(Y, X, alpha = 0.05, omega = 0.2, a = 2,
#'   bounds = list(lambda_m = 1, K_reg = 5, K_eps = 5, K_xi = 20),
#'   setup = list(continuity = TRUE, no_skewness = TRUE) )
#'
#' print(myCI)
#'
#' myCI <- CI.OLS(Y, X, alpha = 0.01, omega = 0.2, a = 2,
#'   bounds = list(lambda_m = 1, K_reg = 5, K_eps = 5, K_xi = 10),
#'   setup = list(continuity = TRUE, no_skewness = TRUE) )
#'
#' print(myCI)
#'
#' @export
#'
CI.OLS <- function(
    Y, X,
    alpha = 0.05,
    omega = NULL, a = NULL,
    bounds = list(lambda_m = NULL,
                  K_reg = NULL,
                  K_eps = NULL,
                  K_xi = NULL,
                  C = NULL,
                  B = NULL),
    setup = list(continuity = FALSE, no_skewness = FALSE),
    regularity = list(C0 = 1,
                      p = 2),
    eps = 0.1,
    options = list(center = TRUE,
                   bounded_case = FALSE),
    matrix_u = diag(NCOL(X)+1) )
{
  # 1- Checking the validity of inputs ==================================

  # Force X to be a matrix, even if there is only one variable
  # Same for matrix_u
  if(NCOL(X) == 1) {X <- matrix(X, ncol = 1)}
  if(NCOL(matrix_u) == 1) {matrix_u <- matrix(matrix_u, ncol = 1)}

  if(!is.vector(Y)) {
    stop("Y must be a vector.")
  }
  if( length(Y) != nrow(X) ){
    stop("Y must have the same number of observations as X.")
  }
  if (ncol(matrix_u) != ncol(X) + 1){
    stop("matrix_u must have exactly one more column than X. ",
         "The first column of matrix_u is then interpreted ",
         "as the coefficient of the intercept. The other columns ",
         "of matrix_u correspond respectively to the columns of X.")
  }
  if (any(unlist(lapply(1:ncol(X), FUN = function(i) {length(unique(X[,i])) == 1})))){
    stop("X must not contain any constant columns.")
  }


  # 2- Computing fundamental quantities that will be useful later ========================

  number_u <- nrow(matrix_u)
  n <- nrow(X)

  # Add a column of ones and take the empirically recentered X
  X <- cbind(matrix(1, nrow = n, ncol = 1),
             scale(X, center = options$center, scale = FALSE) )

  # Estimation of crossproducts and other useful matrices
  XXt <- crossprod(X)
  XXtbar <- (1 / n) * XXt
  minvpXXtbar <- min(eigen(XXtbar, only.values = TRUE)$values)
  inverse_XXt <- solve(XXt)
  inverse_XXtbar <- n * inverse_XXt
  inverse_sqrt_XXtbar <- expm::sqrtm(inverse_XXtbar) # Estimate of E(XX')^{-1/2}

  # X_i tilde
  # this is a list of n elements which are all vectors of size p
  list_Xtilde_i <- lapply(
    X = 1:n,
    FUN = function(i) { inverse_sqrt_XXtbar %*% matrix(X[i, ], ncol = 1) } )

  norms_row_X = apply(X = X, MARGIN = 1, FUN = function(u){sqrt(sum(u^2))})
  norms_row_X_tilde = unlist(lapply(X = list_Xtilde_i,
                                    FUN = function(x) {sqrt(sum(x^2))} ) )

  # Regression
  reg = stats::lm(Y ~ X - 1)
  betahat = reg$coefficients
  OLSestimate_u = matrix_u %*% matrix(betahat, ncol = 1)  # u' OLS estimate

  # Computation of Vhat
  Vhat = n * inverse_XXt %*% (crossprod(X * reg$residuals)) %*% inverse_XXt
  Vhat_u = apply(X = matrix_u, MARGIN = 1,
                 FUN = function(u){t(u) %*% Vhat %*% u})


  # 3- Completing the environment by computing plug-ins for missing bounds ===================

  # Choice of omega, a, and bounds (if they were not given)
  env <- environment()
  OLS.updateTuningParameters(env = env)
  OLS.updateBounds(env = env)

  # Computation of delta
  delta = alpha * omega / 2


  # 4- Computing concentration, Rlin, Rnvar ==================================================

  # Concentration of XXt
  concentr_XXtranspose = concentrationXXt(
    bounded_case = options$bounded_case, bounds = bounds,
    n = n, d = ncol(X), delta = delta)

  # Rn_lin
  Rnlin_u <- RnLin(concentr_XXtranspose = concentr_XXtranspose,
                   bounds = bounds, delta = delta, matrix_u = matrix_u)

  RnVar_u <- RnVar(
    delta = delta, n = n, norms_row_X = norms_row_X,
    residuals = reg$residuals,
    bounds = bounds,
    concentr_XXtranspose = concentr_XXtranspose,
    X = X, inverse_XXtbar = inverse_XXtbar, matrix_u = matrix_u)


  # 5- Computing asymptotic CIs ==============================================================

  # Additional part (at least for the simulations) to return a list with also
  # the plug-in bounds (to have an idea)
  # the usual asymptotic CIs
  result_asymp = matrix(ncol = 2, nrow = number_u)
  colnames(result_asymp) = c("lower", "upper")
  # rownames(result_asymp) = c("intercept", colnames(X))
  CIs.Asymp.extend = (stats::qnorm(1 - alpha/2) / sqrt(n)) * sqrt(Vhat_u)
  result_asymp[,1] = OLSestimate_u - CIs.Asymp.extend
  result_asymp[,2] = OLSestimate_u + CIs.Asymp.extend

  # 6- Preparing the final matrix ============================================================
  result = matrix(nrow = number_u, ncol = 4)
  colnames(result) <- c("lower", "upper", "regime", "estimate")
  # rownames(result) <- c("intercept", colnames(X))
  result = as.data.frame(result)

  # u' OLS estimate (temporary, for check)
  result[,4] = OLSestimate_u

  # 1st regime: R1 =======================================================================

  # First condition (independent of u) for returning \Rb
  if (concentr_XXtranspose >= 1){
    result[,1] = -Inf
    result[,2] = Inf
    result[,3] = "R1"

    return( list(finite = result,
                 asymp = result_asymp,
                 bounds = bounds) )
  }


  # 2nd regime: R2 =======================================================================

  # Second condition for returning \Rb
  nu_nExp_u = OLS.Nu_nExp(alpha = alpha, omega = omega,
                          a = a, K_xi = bounds$K_xi_u, n = n)

  which_regime_R <- which(nu_nExp_u >= alpha / 2)

  if (length(which_regime_R) > 0){

    result[which_regime_R, 1] = -Inf
    result[which_regime_R, 2] = Inf
    result[which_regime_R, 3] = "R2"

    if (length(which_regime_R) == number_u){
      return( list(finite = result,
                   asymp = result_asymp,
                   bounds = bounds) )
    }
  }

  # Preparing for 3rd and 4th regime ==========================================

  setup$iid = TRUE # we are always in the iid framework

  if (!is.null(bounds$K_xi)){
    # K_xi is provided => K_xi and delta_nE are uniform across vectors u
    delta_nE_u = rep(
      BoundEdgeworth::Bound_EE1(
        setup = setup, regularity = regularity, n = n,
        K4 = bounds$K_xi, K3 = NULL, lambda3 = NULL, K3tilde = NULL,
        eps = eps),
      times = number_u)

  } else {
    # Plug-in for K_xi => u-specific K_xi hence delta_nE
    delta_nE_u = apply(
      X = array(1:number_u, dim = number_u),
      MARGIN = 1,
      FUN = function(index_u){
        return( BoundEdgeworth::Bound_EE1(
          setup = setup, regularity = regularity, n = n,
          K4 = bounds$K_xi_u[[index_u]],
          K3 = bounds$K3_xi_u[[index_u]],
          lambda3 = bounds$lambda3_xi_u[[index_u]],
          K3tilde = bounds$K3tilde_xi_u[[index_u]],
          eps = eps) ) } )
  }

  nu_nEdg_u = nu_nExp_u + delta_nE_u

  # Control u' oracle variance u

  bound_Voracle <- Vhat_u + RnVar_u

  nuApprox_u = Rnlin_u / sqrt(bound_Voracle)


  # 3rd regime: Exp =================================================================

  which_regime_Exp = which((nu_nExp_u < alpha / 2) & (nu_nEdg_u >= alpha / 2))

  if (length(which_regime_Exp) > 0){
    CIs.Exp.extend = OLS.CIs.Exp.extend(
      n = n, alpha = alpha, a = a, K_reg = bounds$K_reg, delta = delta,
      nu_nExp = nu_nExp_u[which_regime_Exp],
      nuApprox_u = nuApprox_u[which_regime_Exp],
      bound_Voracle = bound_Voracle[which_regime_Exp])

    result[which_regime_Exp, 1] = OLSestimate_u[which_regime_Exp] - CIs.Exp.extend
    result[which_regime_Exp, 2] = OLSestimate_u[which_regime_Exp] + CIs.Exp.extend
    result[which_regime_Exp, 3] = "Exp"
  }


  # 4th regime: Edg ==============================================================

  which_regime_Edg = which(nu_nEdg_u < alpha / 2)
  # Implies automatically nu_nExp_u < alpha / 2 as well.

  if (length(which_regime_Edg) > 0){
    CIs.Edg.extend = OLS.CIs.Edg.extend(
      n = n, alpha = alpha, a = a, K_reg = bounds$K_reg, delta = delta,
      nu_nEdg = nu_nEdg_u[which_regime_Edg],
      nuApprox_u = nuApprox_u[which_regime_Edg],
      bound_Voracle = bound_Voracle[which_regime_Edg])

    result[which_regime_Edg, 1] = OLSestimate_u[which_regime_Edg] - CIs.Edg.extend
    result[which_regime_Edg, 2] = OLSestimate_u[which_regime_Edg] + CIs.Edg.extend
    result[which_regime_Edg, 3] = "Edg"
  }

  return( list(finite = result,
               asymp = result_asymp,
               bounds = bounds) )
}

OLS.CIs.Exp.extend <- function(
    n, alpha, a, K_reg, delta,
    nu_nExp, nuApprox_u, bound_Voracle)
{

  part1_Qn = sqrt( 2 * (1 + a) * ( 1 - log(alpha/2 - nu_nExp) ) )

  QnExp_u = part1_Qn + nuApprox_u

  result = ( QnExp_u / sqrt(n) ) * sqrt(bound_Voracle)

  return(result)
}

OLS.CIs.Edg.extend <- function(
    n, alpha, a, K_reg, delta,
    nu_nEdg, nuApprox_u, bound_Voracle)
{

  part1_Qn = sqrt(a) * stats::qnorm(1 - alpha/2 + nu_nEdg)

  QnEdg_u = part1_Qn + nuApprox_u

  result = ( QnEdg_u / sqrt(n) ) * sqrt(bound_Voracle)

  return(result)
}

#' Computation of nu_nExp
#'
#' @noRd
#'
OLS.Nu_nExp <- function(alpha, omega, a, K_xi, n)
{
  result = (omega * alpha + exp(-n * (1 - 1/a)^2 / (2 * K_xi) ) ) / 2
  return (result)
}


RnLin <- function(concentr_XXtranspose, bounds, delta, matrix_u){

  cc <- concentr_XXtranspose

  RnLin_without_norm_u <-
    sqrt(2) * cc /
    ((1 - cc) * sqrt(bounds$lambda_m)) *
    (bounds$K_eps / delta)^(1/4)

  Rnlin_u <- apply(X = matrix_u, MARGIN = 1,
                   FUN = function(u){sqrt(sum(u^2))}) * RnLin_without_norm_u

  return(Rnlin_u)
}


RnVar <- function(
    delta, n, norms_row_X, residuals,
    bounds,
    concentr_XXtranspose,
    X, inverse_XXtbar, matrix_u)
{

  cc <- concentr_XXtranspose
  lambda_m <- bounds$lambda_m
  K_eps <- bounds$K_eps

  part1 = (2 / (n * lambda_m^3)) *
    (cc / (1 - cc) + 1)^2 *
    sqrt(K_eps / delta) *
    mean(norms_row_X^4)

  part2 = (2 * sqrt(2) / (lambda_m^(5/2) * sqrt(n))) *
    (cc / (1 - cc) + 1) *
    (K_eps / delta)^(1/4) *
    mean(norms_row_X^3 * abs(residuals))

  part3 = (1 / lambda_m^2) * (cc / (1 - cc))^2 *
    mean(norms_row_X^2 * residuals^2)

  part4 = 2 / lambda_m * (cc / (1 - cc)) *
    base::norm(x = n^(-1) * crossprod(X * residuals) %*% inverse_XXtbar,
               type = "2")

  result = apply(X = matrix_u, MARGIN = 1,
                 FUN = function(u){sum(u^2)}) *
    (part1 + part2 + part3 + part4)

  return(result)
}



concentrationXXt <- function(bounded_case, bounds, n, d, delta)
{
  if (bounded_case){
    # Concentration of XX transpose, assuming bounded regressors

    concentr_XXtranspose = concentration_Bernstein(
      B = bounds$B, C = bounds$C, delta = delta, n = n, d = d)

  } else {
    # Concentration of XX transpose without assuming bounded regressors

    concentr_XXtranspose <- sqrt(bounds$K_reg / (n * delta))
  }

  return (concentr_XXtranspose)
}


#' Compute the concentration bound from Lemma C.2
#' (Bernstein-type lemma for concentration of random, bounded matrices)
#'
#' @param B \eqn{:= || \E [ (A-\E[A]) (A-\E[A]) ] ||}.
#' @param C almost sure upper bound on the operator 2-norm of A.
#' @param delta confidence level
#' @param n sample size
#' @param d dimension of A
#'
concentration_Bernstein <- function(B, C, delta, n, d)
{
  return( sqrt(2 * B * log(2 * d / delta) / n) +
            4 * C * log(2 * d / delta) / (3*n) )
}
