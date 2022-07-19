

CI.OLS_new <- function(Y_data, X_data, alpha = 0.05, omega, a,
                       choice_new_bound_RnLin,
                       choice_new_bound_RnVar,
                       M = NULL,
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

  stopifnot(choice_new_bound_RnLin %in% c("1", "2"))
  # "1" : Equation (1) temp
  # "2" : Equation (2) temp (sans le terme multiplicatif)

  stopifnot(choice_new_bound_RnVar %in% c("1", "2"))
  # "1" : 4 termes, haut page 2
  # "2" : avec astuce majoration |ab| par demi-somme des carres

  number_u = NROW(matrix_u)

  n = nrow(X_data)

  # Take the empirically recentred X
  X_centered_wo_cst = scale(X_data[,-1], center = TRUE, scale = FALSE)
  X = cbind(matrix(1, nrow = n, ncol = 1), X_centered_wo_cst)
  Y = Y_data

  XXt <- crossprod(X)
  XXtbar = (1 / n) * XXt
  minvpXXtbar = min(eigen(XXtbar, only.values = TRUE)$values)
  inverse_XXt <- solve(XXt)
  inverse_XXtbar <- n * inverse_XXt
  # Estimate of E(XX')^{-1/2}
  inverse_sqrt_XXtbar <- expm::sqrtm(inverse_XXtbar)

  # X_i tilde

  list_Xtilde_i <- purrr::map(1:n,
                              ~ inverse_sqrt_XXtbar %*% matrix(X[.x,], ncol = 1))


  # Bound lambda_m
  bounds$lambda_m = minvpXXtbar

  # Bound K_X (K_reg)

  # To mean across observations i
  To_mean_over_obs_i_for_KX_new <- function(index_obs_i, list_Xtilde_i, p){
    return(
      sum( c(
        list_Xtilde_i[[index_obs_i]] %*% t(list_Xtilde_i[[index_obs_i]]) -
          diag(x = 1, nrow = p, ncol = p) )^2 )
    )
  }
  veca_To_mean_over_obs_i_for_KX_new <- purrr::map_dbl(
    1:n, To_mean_over_obs_i_for_KX_new,
    list_Xtilde_i = list_Xtilde_i, p = ncol(X))

  bound_KX_estimated_plug_in <- mean(veca_To_mean_over_obs_i_for_KX_new)
  bounds$K_X = bound_KX_estimated_plug_in

  # Regression
  reg = stats::lm(Y ~ X - 1)
  betahat = reg$coefficients

  norms_row_X = apply(X = X, MARGIN = 1, FUN = function(u){sqrt(sum(u^2))})

  norms_row_X_tilde = purrr::map_dbl(list_Xtilde_i, ~ sqrt(sum(c(.x)^2)))

  # Bound K_eps
  bounds$K_eps <- mean( norms_row_X_tilde^4 * reg$residuals^4 )

  # Bound Xi
  bounds$K_xi_u = rep(bounds$K_xi, length.out = number_u)

  # Vhat

  Vhat = n * inverse_XXt %*% (crossprod(X * reg$residuals)) %*% inverse_XXt

  Vhat_u = apply(X = matrix_u, MARGIN = 1,
                 FUN = function(u){t(u) %*% Vhat %*% u})

  # Delta, omega, a
  delta = alpha * omega / 2

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
  if (sqrt(bounds$K_X / (n * delta)) >= 1){
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

  # Control of linearization error

  if (choice_new_bound_RnLin == "1"){

    Rnlin_u = apply(X = matrix_u, MARGIN = 1,
                    FUN = function(u){sum(u^2)}) *
      OLS.RnLin_new1_part1(
        delta = delta, n = n,
        lambda_m = bounds$lambda_m, K_X = bounds$K_X , K_eps = bounds$K_eps) +
      OLS.RnLin_new1_part2(
        delta = delta, n = n,
        lambda_m = bounds$lambda_m, K_X = bounds$K_X , K_eps = bounds$K_eps)

  } else {

    Rnlin_u = apply(X = matrix_u, MARGIN = 1,
                    FUN = function(u){sqrt(sum(u^2))}) *
      OLS.RnLin_new2(
        delta = delta, n = n,
        lambda_m = bounds$lambda_m, K_X = bounds$K_X , K_eps = bounds$K_eps)

  }

  # Control u' oracle variance u

  if (choice_new_bound_RnVar == "1"){

    bound_Voracle <-
      Vhat_u +
      apply(X = matrix_u, MARGIN = 1,
            FUN = function(u){sum(u^2)}) *
      OLS.RnVar_new1(
        delta = delta, n = n, norms_row_X = norms_row_X,
        residuals = reg$residuals,
        lambda_m = bounds$lambda_m, K_X = bounds$K_X, K_eps = bounds$K_eps)

  } else {

    bound_Voracle <-
      (1 + M^2)^2 * Vhat_u +
      apply(X = matrix_u, MARGIN = 1,
            FUN = function(u){sum(u^2)}) *
      OLS.RnVar_new2(
        M = M,
        delta = delta, n = n, norms_row_X = norms_row_X,
        residuals = reg$residuals,
        lambda_m = bounds$lambda_m, K_X = bounds$K_X, K_eps = bounds$K_eps)

  }

  nuApprox_u = Rnlin_u / sqrt(bound_Voracle)

  # Final part: computation of the confidence intervals
  which_regime_Exp = which((nu_nExp_u < alpha / 2) & (nu_nEdg_u >= alpha / 2))

  if (length(which_regime_Exp) > 0){
    CIs.Exp.extend = OLS.CIs.Exp.extend_new(
      n = n, alpha = alpha, a = a, K_X = bounds$K_X, delta = delta,
      nu_nExp = nu_nExp_u[which_regime_Exp],
      nuApprox_u = nuApprox_u[which_regime_Exp],
      bound_Voracle = bound_Voracle[which_regime_Exp],
      choice_new_bound_RnLin = choice_new_bound_RnLin)

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
      bound_Voracle = bound_Voracle[which_regime_Edg],
      choice_new_bound_RnLin = choice_new_bound_RnLin)

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

  if (choice_new_bound_RnLin == "1"){

    QnExp_u = (1 +  sqrt(K_X / (n * delta))) * part1_Qn + nuApprox_u

  } else {

    QnExp_u = part1_Qn + nuApprox_u

  }

  result = ( QnExp_u / sqrt(n) ) * sqrt(bound_Voracle)

  return(result)
}

OLS.CIs.Edg.extend_new <- function(
  n, alpha, a, K_X, delta,
  nu_nEdg, nuApprox_u, bound_Voracle,
  choice_new_bound_RnLin)
{

  part1_Qn = sqrt(a) * stats::qnorm(1 - alpha/2 + nu_nEdg)

  if (choice_new_bound_RnLin == "1"){

    QnEdg_u = (1 +  sqrt(K_X / (n * delta))) * part1_Qn + nuApprox_u

  } else {

    QnEdg_u = part1_Qn + nuApprox_u

  }

  result = ( QnEdg_u / sqrt(n) ) * sqrt(bound_Voracle)

  return(result)
}


OLS.Rnlin_bounded_case <- function(

)
{




}



OLS.RnLin_new1_part1 <- function(
  delta, n, lambda_m, K_X, K_eps)
{
  sqrt_K_X_div_n_delta <- sqrt(K_X / (n * delta))
  sq <- sqrt_K_X_div_n_delta

  return(sq / (2 * lambda_m))
}

OLS.RnLin_new1_part2 <- function(
  delta, n, lambda_m, K_X, K_eps)
{
  sqrt_K_X_div_n_delta <- sqrt(K_X / (n * delta))
  sq <- sqrt_K_X_div_n_delta

  return((sq / n) * sqrt(K_eps / delta))
}

OLS.RnLin_new2 <- function(
  delta, n, lambda_m, K_X, K_eps)
{
  sqrt_K_X_div_n_delta <- sqrt(K_X / (n * delta))
  sq <- sqrt_K_X_div_n_delta

  return(sqrt(2 * K_X / delta) / (n * (1 - sq) * sqrt(lambda_m)) *
           (K_eps / delta)^(1/4))
}


OLS.RnVar_new1 <- function(
  delta, n, norms_row_X, residuals,
  lambda_m, K_X, K_eps)
{

  sqrt_K_X_div_n_delta <- sqrt(K_X / (n * delta))
  sq <- sqrt_K_X_div_n_delta

  part1 = (2 / (n * lambda_m^2)) *
    (sq / (1 - sq) + 1)^2 *
    sqrt(K_eps / delta) *
    mean(norms_row_X^4)

  part2 = (2 * sqrt(2) / (lambda_m^2 * sqrt(n))) *
    (sq / (1 - sq) + 1) *
    (K_eps / delta)^(1/4) *
    mean(norms_row_X^3 * abs(residuals))

  mean_norms_square_X_residuals <- mean(norms_row_X^2 * residuals^2)

  part3 = (K_X / (n * delta * lambda_m^2)) *
    (1 - sq)^(-2) * sq^2 *
    mean_norms_square_X_residuals

  part4 = 2 / lambda_m *
    (1 - sq)^(-2) * sq *
    mean_norms_square_X_residuals

  return(part1 + part2 + part3 + part4)
}

OLS.RnVar_new2 <- function(
  M, delta, n, norms_row_X, residuals,
  lambda_m, K_X, K_eps)
{

  sqrt_K_X_div_n_delta <- sqrt(K_X / (n * delta))
  sq <- sqrt_K_X_div_n_delta

  part1 = (1 + M^2) * (1 + M^(-2)) *
    sq^2 / (lambda_m^2 * (1 - sq)^2) *
    mean(norms_row_X^2 * residuals^2)

  part2 = (2 * (1 + M^(-2)) / (n * lambda_m^2)) *
    (sq / (1 - sq) + 1)^2 *
    sqrt(K_eps / delta) *
    mean(norms_row_X^4)

  return(part1 + part2)
}

