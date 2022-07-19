
#' Compute valid CI for an ordinary least squares regression
#'
#' @param Y vector of observations of the explained variables
#' @param X matrix of explanatory variables
#' @param alpha 1 - level of confidence of the CI
#' @param omega the tuning parameter \eqn{\omega} of the interval
#' @param a the tuning parameter \eqn{a} of the interval
#' @param bounds list of bounds for the DGP.
#' It can contain the following items: \itemize{
#'    \item lambda_m
#'    \item K_X
#'    \item K_eps
#'    \item K_xi, K3_xi, lambda3_xi, K3tilde_xi
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
#' n = 20
#' X1 = rnorm(n, sd = 4)
#' true_eps = rnorm(n)
#' Y = 3 + 8 * X1 + true_eps
#' X = cbind(intercept = 1, X1)
#'
#' myCI <- CI.OLS(Y, X, alpha = 0.05, omega = 0.2, a = 2,
#'   bounds = list(lambda_m = 1, K_X = 5, K_eps = 5, K_xi = 10),
#'   setup = list(continuity = FALSE, no_skewness = FALSE) )
#'
#' print(myCI)
#'
#' myCI <- CI.OLS(Y, X, alpha = 0.05, omega = 0.2, a = 2,
#'   bounds = list(lambda_m = 1, K_X = 5, K_eps = 5, K_xi = 20),
#'   setup = list(continuity = TRUE, no_skewness = TRUE) )
#'
#' print(myCI)
#'
#' myCI <- CI.OLS(Y, X, alpha = 0.01, omega = 0.2, a = 2,
#'   bounds = list(lambda_m = 1, K_X = 5, K_eps = 5, K_xi = 10),
#'   setup = list(continuity = TRUE, no_skewness = TRUE) )
#'
#' print(myCI)
#'
#' @export
#'
CI.OLS <- function(Y_data, X_data, alpha = 0.05, omega, a,
                   regression_made,
                   bounds = list(lambda_m = NULL,
                                 K_X = NULL,
                                 K_eps = NULL,
                                 K_xi = NULL),
                   setup = list(continuity = FALSE, no_skewness = FALSE),
                   regularity = list(C0 = 1,
                                     p = 2,
                                     kappa = 0.99),
                   eps = 0.1,
                   use_uniform_bounds = TRUE,
                   matrix_u = diag(NCOL(X_data)) )
{
  if( NROW(Y_data) != NROW(X_data) ){
    stop("Y_data should have the same number of observations as X_data.")
  }
  n = NROW(Y_data)
  if (NCOL(matrix_u) != NCOL(X_data)){
    stop("matrix_u should have the same number of columns as X_data")
  }
  number_u = NROW(matrix_u)

  # Regression en ecart a la moyenne pour de meilleures bornes plug-in
  # Differentes possibilites
  if (regression_made == "Y_X"){
    # Regresssion of given Y on given X with a constant in X implicitly
    X = X_data
    Y = Y_data
  } else if (regression_made == "Y_X_centered"){
    # Regression of Y recentered over X recentered (with a constant/intercept)
    X_centered_wo_cst = scale(X_data[,-1], center = TRUE, scale = FALSE)
    X = cbind(matrix(1, nrow = n, ncol = 1), X_centered_wo_cst)
    Y = scale(Y_data, center = TRUE, scale = FALSE)
  } else if (regression_made == "Y_X_centered_wo_cst"){
    # Regression of Y recentered over X recentered without a constant/intercept
    X = scale(X_data[,-1], center = TRUE, scale = FALSE)
    Y = scale(Y_data, center = TRUE, scale = FALSE)
    # Need to delete first row and column of u corresponding to the constant in this case
    matrix_u <- matrix_u[-1,]
    matrix_u <- matrix_u[,-1]
  } else if (regression_made == "X_centered"){
    # Regression of Y (as initial data) over X recentered (with a constant/intercept)
    X_centered_wo_cst = scale(X_data[,-1], center = TRUE, scale = FALSE)
    X = cbind(matrix(1, nrow = n, ncol = 1), X_centered_wo_cst)
    Y = Y_data
  } else {
    stop("'regression_made' should be one of 'Y_X', 'Y_X_centered', 'Y_X_centered_wo_cst', 'X_centered'.")
  }

  # - 1 in the formula of lm() because if wanted X should contain an intercept
  reg = stats::lm(Y ~ X - 1)

  betahat = reg$coefficients

  XXt <- crossprod(X)
  XXtbar = (1 / n) * XXt
  minvpXXtbar = min(eigen(XXtbar, only.values = TRUE)$values)
  inverse_XXt <- solve(XXt)
  inverse_XXtbar <- n * inverse_XXt

  norms_row_X = apply(X = X, MARGIN = 1, FUN = function(u){sqrt(sum(u^2))})

  Vhat = n * inverse_XXt %*% (crossprod(X * reg$residuals)) %*% inverse_XXt

  Vhat_u = apply(X = matrix_u, MARGIN = 1,
                 FUN = function(u){t(u) %*% Vhat %*% u})

  # Replacing NULL values by plug-ins whenever necessary
  env <- environment()
  OLS.updateBounds(env)
  # After this, we can safely assume that all bounds have the right length \code{number_u}.

  # Replacing NULL tuning parameters by default choices
  OLS.updateTuningParameters(env)


  # Preparing the final matrix
  result = matrix(nrow = number_u, ncol = 4)
  colnames(result) <- c("lower", "upper", "regime", "estimate")
  rownames(result) <- colnames(X)
  result = as.data.frame(result)

  # u' OLS estimate (temporary, for check)
  result[,4] = matrix_u %*% matrix(betahat, ncol = 1)

  # Additional part (at least for the simulations) to return a list with also
  # the plug-in bounds (to have an idea)
  # the usual asymptotic CIs
  result_asymp = matrix(ncol = 2, nrow = number_u)
  rownames(result_asymp) = rownames(result)
  colnames(result_asymp) = colnames(result)[1:2]
  CIs.Asymp.extend = (stats::qnorm(1 - alpha/2) / sqrt(n)) * sqrt(Vhat_u)
  result_asymp[,1] = result$estimate - CIs.Asymp.extend
  result_asymp[,2] = result$estimate + CIs.Asymp.extend

  # First condition (independent of u) for returning \Rb
  if ( n <= ( 2 * bounds$K_X / ( omega * alpha * bounds$lambda_m^2 ) ) ){
    result[,1] = -Inf
    result[,2] = Inf
    result[,3] = "R1"

    return( list(finite = result,
                 asymp = result_asymp,
                 bounds = bounds) )
  }

  # Second condition for returning \Rb
  nu_nExp_u = OLS.Nu_nExp(alpha = alpha, omega = omega,
                          a = a, K_xi = bounds$K_xi_u, n = n)

  which_regime_R <- which((nu_nExp_u >= alpha / 2))

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

  # V1
  # Proposition below to avoid ... and rely on default behavior
  # of Bound_Edgeworth_Expansion if NULL for K3, kambda3, K3tilde

  # delta_nE_u = apply(
  #   FUN = function(i_u, ...) {
  #     return(Bound_Edgeworth_Expansion(
  #       K4 = bounds$K_xi_u[[i_u]], K3 = bounds$K3_xi_u[i_u],
  #       lambda3 = bounds$lambda3_xi_u[i_u], K3tilde = bounds$K3tilde_xi_u[i_u], ...))},
  #   MARGIN = 1,
  #   X = array(1:number_u,dim = number_u),
  #   continuity = args_Edg$continuity, iid = TRUE, noskewness = args_Edg$noskewness,
  #   n = n,
  #   bound_C0_fSn_cont_inid = args_Edg$bound_C0_fSn_cont_inid,
  #   bound_p_fSn_cont_inid = args_Edg$bound_p_fSn_cont_inid,
  #   bound_kappa_fXoversigma_cont_iid = args_Edg$bound_kappa_fXoversigma_cont_iid,
  #   eps = args_Edg$eps,
  #   use_uniform_bounds_in_Omega1_Omega2 = args_Edg$uniform_bounds)

  setup$iid = TRUE # Always in the iid framework

  if (!is.null(bounds$K_xi)){
  # K_xi is provided => K_xi and delta_nE are uniform across vectors u
    delta_nE_u = rep(
      Bound_Edgeworth_Expansion(
        setup = setup, regularity = regularity, n = n,
        K4 = bounds$K_xi, K3 = NULL, lambda3 = NULL, K3tilde = NULL,
        use_uniform_bounds = use_uniform_bounds, eps = eps),
      times = number_u)

  } else {
  # Plug-in for K_xi => u-specific K_xi hence delta_nE
    delta_nE_u = apply(
      X = array(1:number_u, dim = number_u),
      MARGIN = 1,
      FUN = function(index_u){
        return( Bound_Edgeworth_Expansion(
          setup = setup, regularity = regularity, n = n,
          K4 = bounds$K_xi_u[[index_u]],
          K3 = bounds$K3_xi_u[[index_u]],
          lambda3 = bounds$lambda3_xi_u[[index_u]],
          K3tilde = bounds$K3tilde_xi_u[[index_u]],
          use_uniform_bounds = use_uniform_bounds,
          eps = eps) ) } )
  }

  nu_nEdg_u = nu_nExp_u + delta_nE_u

  basic_Rnlin = OLS.basic_RnLin(
    delta = alpha * omega / 2, n = n, minvpXXtbar = minvpXXtbar,
    lambda_m = bounds$lambda_m, K_X = bounds$K_X , K_eps = bounds$K_eps)

  basic_Rnvar = OLS.basic_RnVar(
    delta = alpha * omega / 2, n = n,
    norms_row_X = norms_row_X, minvpXXtbar = minvpXXtbar,
    residuals = reg$residuals,
    lambda_m = bounds$lambda_m, K_X = bounds$K_X, K_eps = bounds$K_eps)

  # Computation of the Rn and final quantities depending on u

  Rnlin_u = apply(X = matrix_u, MARGIN = 1,
                  FUN = function(u){sqrt(sum(u^2))}) * basic_Rnlin

  RnVar_u = apply(X = matrix_u, MARGIN = 1,
                  FUN = function(u){sum(u^2)}) * basic_Rnvar

  nuApprox_u = Rnlin_u / sqrt(Vhat_u + RnVar_u)

  # Final part: computation of the confidence intervals
  which_regime_Exp = which((nu_nExp_u < alpha / 2) & (nu_nEdg_u >= alpha / 2))

  if (length(which_regime_Exp) > 0){
    CIs.Exp.extend = OLS.CIs.Exp.extend(
      n = n, alpha = alpha, a = a,
      nu_nExp = nu_nExp_u[which_regime_Exp],
      nuApprox_u = nuApprox_u[which_regime_Exp],
      Vhat_u = Vhat_u[which_regime_Exp],
      RnVar_u = RnVar_u[which_regime_Exp])

    result[which_regime_Exp, 1] =
      matrix_u[which_regime_Exp, ] %*% matrix(betahat, ncol = 1) - CIs.Exp.extend

    result[which_regime_Exp, 2] =
      matrix_u[which_regime_Exp, ] %*% matrix(betahat, ncol = 1) + CIs.Exp.extend

    result[which_regime_Exp, 3] = "Exp"
  }

  which_regime_Edg = which(nu_nEdg_u < alpha / 2)
  # Implies automatically nu_nExp_u < alpha / 2 as well.

  if (length(which_regime_Edg) > 0){
    CIs.Edg.extend = OLS.CIs.Edg.extend(
      n = n, alpha = alpha, a = a,
      nu_nEdg = nu_nEdg_u[which_regime_Edg],
      nuApprox_u = nuApprox_u[which_regime_Edg],
      Vhat_u = Vhat_u[which_regime_Edg],
      RnVar_u = RnVar_u[which_regime_Edg])

    result[which_regime_Edg, 1] =
      matrix_u[which_regime_Edg, ] %*% matrix(betahat, ncol = 1) - CIs.Edg.extend

    result[which_regime_Edg, 2] =
      matrix_u[which_regime_Edg, ] %*% matrix(betahat, ncol = 1) + CIs.Edg.extend

    result[which_regime_Edg, 3] = "Edg"
  }

  return( list(finite = result,
               asymp = result_asymp,
               bounds = bounds) )
}


#==============================================================================#
#### Miscellaneous functions ####
#==============================================================================#

#' Computation of nu_nExp
#'
#' @noRd
#'
OLS.Nu_nExp <- function(alpha, omega, a, K_xi, n)
{
  result = (omega * alpha + exp(-n * (1 - 1/a)^2 / (2 * K_xi) ) ) / 2
  return (result)
}

#' Half-length of the CI for the OLS with the Exp regime
#'
#' @noRd
OLS.CIs.Exp.extend <- function(
  n, alpha, a,
  nu_nExp, nuApprox_u, Vhat_u, RnVar_u)
{
  QnExp_u = sqrt( 2 * (1 + a) * ( 1 - log(alpha/2 - nu_nExp) ) ) + nuApprox_u

  result = ( QnExp_u / sqrt(n) ) * sqrt(Vhat_u + RnVar_u)

  return(result)
}

#' Half-length of the CI for the OLS with the Edg regime
#'
#' @noRd
OLS.CIs.Edg.extend <- function(
  n, alpha, a,
  nu_nEdg, nuApprox_u, Vhat_u, RnVar_u)
{
  QnEdg_u = sqrt(a) * stats::qnorm(1 - alpha/2 + nu_nEdg) + nuApprox_u

  result = ( QnEdg_u / sqrt(n) ) * sqrt(Vhat_u + RnVar_u)

  return(result)
}

# "basic": without the factor ||u||
OLS.basic_RnVar <- function(
  delta, n, norms_row_X, minvpXXtbar, residuals,
  lambda_m, K_X, K_eps)
{
  part1 = (2 / (n * lambda_m^2 * minvpXXtbar^2)) *
    sqrt(K_eps / delta) * mean(norms_row_X^4)

  part2 = (2 * sqrt(2) / (lambda_m^2 * minvpXXtbar * sqrt(n))) *
    (K_eps / delta)^(1/4) *
    mean(norms_row_X^3 * abs(residuals))

  mean_norms_square_X_residuals <- mean(norms_row_X^2 * residuals^2)

  part3 = (K_X / (n * delta * lambda_m^2 * minvpXXtbar^2)) *
    mean_norms_square_X_residuals

  part4 = 2 / (lambda_m * minvpXXtbar^2) *
    sqrt(K_X / (n * delta)) * mean_norms_square_X_residuals

  return(part1 + part2 + part3 + part4)
}

# "basic": without the factor ||u||
OLS.basic_RnLin <- function(
  delta, n, minvpXXtbar, lambda_m, K_X, K_eps)
{
  basic_Rnlin = 1 / (lambda_m * minvpXXtbar) *
    sqrt(2 * K_X / (n * delta)) *
    (K_eps / delta)^{1/4}

  return(basic_Rnlin)
}

# EOF
