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
#' X = cbind(X1)
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
CI.OLS <- function(
    Y, X,
    alpha = 0.05,
    omega = NULL, a = NULL,
    C = NULL, # Borne sur les || Xi tilde %*% Xi tilde'||
    bounds = list(lambda_m = NULL,
                  K_X = NULL,
                  K_eps = NULL,
                  K_xi = NULL),
    setup = list(continuity = FALSE, no_skewness = FALSE),
    regularity = list(C0 = 1,
                      p = 2,
                      kappa = 0.99),
    eps = 0.1,
    options = list(center = TRUE,
                   scale = FALSE,
                   bounded_case = FALSE),
    matrix_u = diag(NCOL(X)+1) )
{
  # Force X to be a matrix, even if there is only one variable
  # Same for matrix_u
  if(NCOL(X) == 1) {X <- matrix(X, ncol = 1)}
  if(NCOL(matrix_u) == 1) {matrix_u <- matrix(matrix_u, ncol = 1)}

  if(!is.vector(Y)) {
    stop("Y should be a vector.")
  }
  if( nrow(Y) != nrow(X) ){
    stop("Y should have the same number of observations as X.")
  }
  if (ncol(matrix_u) != ncol(X) + 1){
    stop("matrix_u should have exactly one more column than X. ",
         "The first column of matrix_u is then interpreted ",
         "as the coefficient of the intercept.")
  }
  if (any(unlist(lapply(1:ncol(X), FUN = function(i) {length(unique(X[,i])) == 1})))){
    stop("X should not contain any constant columns.")
  }

  number_u <- nrow(matrix_u)
  n <- nrow(X)

  # Choice of omega and a, if they are not provided yet.
  env <- environment()
  OLS.updateTuningParameters(env = env, bounded_case = bounded_case)

  # Add a column of ones and take the empirically recentered X
  X <- cbind(matrix(1, nrow = n, ncol = 1),
            scale(X, center = options$center, scale = options$scale))

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
                                    FUN = function(x) {sqrt(sum((x)^2))} ) )

  # Regression
  reg = stats::lm(Y ~ X - 1)
  betahat = reg$coefficients

  # Updating bounds (if they were not given)
  OLS.updateBounds(env)

  # Delta, omega, a
  delta = alpha * omega / 2

  # Control linearization term Rn lin dans le cas borne

  if (isTRUE(bounded_case)){

    # Borne C
    if (is.null(C)){
      C = max(norms_row_X_tilde^2)
    }

    # Plug-in de la borne B (A = X_i tilde X_i tilde transpose)
    list_A_i = purrr::map(1:n,
                          ~ matrix(list_Xtilde_i[[.x]]) %*%
                            t(matrix(list_Xtilde_i[[.x]])))

    # connu theoriquement par definition des X_i tilde
    d = ncol(X)
    expectation_A = diag(x = 1, nrow = d, ncol = d)

    # Liste des (A - E(A))(A - E(A))
    list_A_mEA_sq = purrr::map(1:n,
                               ~ (list_A_i[[.x]] - expectation_A) %*%
                                 (list_A_i[[.x]] - expectation_A))

    B_before_norm = purrr::reduce(list_A_mEA_sq, `+`, .init = matrix(0, ncol = d, nrow = d)) / n
    B = base::norm(x = B_before_norm, type = "2")

    concentr_XXtranspose = sqrt(2 * B * log(2 * d / delta) / n) +
      4 * C * log(2 * d / delta) / (3*n)

  } else {

    concentr_XXtranspose <- sqrt(bounds$K_X / (n * delta))
  }

  # without the product by norm of u
  base_Rnlin = sqrt(2 / n) * concentr_XXtranspose /
    ((1 - concentr_XXtranspose) * sqrt(bounds$lambda_m)) *
    (bounds$K_eps / delta)^(1/4)

  Rnlin_u = apply(X = matrix_u, MARGIN = 1,
                  FUN = function(u){sqrt(sum(u^2))}) * base_Rnlin

  # Vhat

  Vhat = n * inverse_XXt %*% (crossprod(X * reg$residuals)) %*% inverse_XXt

  Vhat_u = apply(X = matrix_u, MARGIN = 1,
                 FUN = function(u){t(u) %*% Vhat %*% u})

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
  if (concentr_XXtranspose >= 1){
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

  # Control u' oracle variance u

  bound_Voracle <-
    Vhat_u +
    apply(X = matrix_u, MARGIN = 1,
          FUN = function(u){sum(u^2)}) *
    OLS.RnVar_bounded(
      delta = delta, n = n, norms_row_X = norms_row_X,
      residuals = reg$residuals,
      lambda_m = bounds$lambda_m, K_X = bounds$K_X, K_eps = bounds$K_eps,
      concentr_XXtranspose = concentr_XXtranspose,
      X = X, inverse_XXtbar = inverse_XXtbar)

  nuApprox_u = Rnlin_u / sqrt(bound_Voracle)

  # Final part: computation of the confidence intervals
  which_regime_Exp = which((nu_nExp_u < alpha / 2) & (nu_nEdg_u >= alpha / 2))

  if (length(which_regime_Exp) > 0){
    CIs.Exp.extend = OLS.CIs.Exp.extend_new(
      n = n, alpha = alpha, a = a, K_X = bounds$K_X, delta = delta,
      nu_nExp = nu_nExp_u[which_regime_Exp],
      nuApprox_u = nuApprox_u[which_regime_Exp],
      bound_Voracle = bound_Voracle[which_regime_Exp])

    result[which_regime_Exp, 1] =
      matrix_u[which_regime_Exp, ] %*% matrix(betahat, ncol = 1) - CIs.Exp.extend

    result[which_regime_Exp, 2] =
      matrix_u[which_regime_Exp, ] %*% matrix(betahat, ncol = 1) + CIs.Exp.extend

    result[which_regime_Exp, 3] = "Exp"
  }

  which_regime_Edg = which(nu_nEdg_u < alpha / 2)
  # Implies automatically nu_nExp_u < alpha / 2 as well.

  if (length(which_regime_Edg) > 0){
    CIs.Edg.extend = OLS.CIs.Edg.extend_new(
      n = n, alpha = alpha, a = a, K_X = bounds$K_X, delta = delta,
      nu_nEdg = nu_nEdg_u[which_regime_Edg],
      nuApprox_u = nuApprox_u[which_regime_Edg],
      bound_Voracle = bound_Voracle[which_regime_Edg])

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

OLS.CIs.Exp.extend_new <- function(
  n, alpha, a, K_X, delta,
  nu_nExp, nuApprox_u, bound_Voracle,
  choice_new_bound_RnLin)
{

  part1_Qn = sqrt( 2 * (1 + a) * ( 1 - log(alpha/2 - nu_nExp) ) )

  QnExp_u = part1_Qn + nuApprox_u

  result = ( QnExp_u / sqrt(n) ) * sqrt(bound_Voracle)

  return(result)
}

OLS.CIs.Edg.extend_new <- function(
  n, alpha, a, K_X, delta,
  nu_nEdg, nuApprox_u, bound_Voracle,
  choice_new_bound_RnLin)
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

OLS.RnVar_bounded <- function(
  delta, n, norms_row_X, residuals,
  lambda_m, K_X, K_eps,
  concentr_XXtranspose,
  X, inverse_XXtbar)
{

  cc <- concentr_XXtranspose

  part1 = (2 / (n * lambda_m^3)) *
    (cc / (1 - cc) + 1)^2 *
    sqrt(K_eps / delta) *
    mean(norms_row_X^4)

  part2 = (2 * sqrt(2) / (lambda_m^(5/2) * sqrt(n))) *
    (cc / (1 - cc) + 1) *
    (K_eps / delta)^(1/4) *
    mean(norms_row_X^3 * abs(residuals))

  mean_norms_square_X_residuals <- mean(norms_row_X^2 * residuals^2)
  part3 = (1 / lambda_m^2) * (cc / (1 - cc))^2 *
    mean_norms_square_X_residuals

  part4 = 2 / lambda_m *
    (1 - cc)^(-1) * cc *
    base::norm(x = n^(-1) * crossprod(X * residuals) %*% inverse_XXtbar,
               type = "2")
  return(part1 + part2 + part3 + part4)
}
