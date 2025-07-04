#' Compute NAVAE CI for coefficients of a linear regression based on
#' the OLS estimator and Berry-Esseen (BE) or Edgeworth Expansions (EE) bounds

#' @param Y vector of observations of the explained variables
#' @param X matrix of explanatory variables.
#' \code{X} must not contain a constant column.
#' @param alpha 1 - level of confidence of the CI
#' @param omega the tuning parameter \eqn{\omega} of the interval
#' @param a the tuning parameter \eqn{a} of the interval
#' @param bounds list of bounds for the DGP.
#' It can contain the following items: \itemize{
#'    \item \code{lambda_reg}
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
#' n = 4000
#' X1 = rnorm(n, sd = 1)
#' true_eps = rnorm(n)
#' Y = 8 * X1 + true_eps
#' X = cbind(X1)
#'
#' myCI <- Navae_ci_ols(Y, X, K_xi = 3, intercept = TRUE, verbose = 2, a = 1.1)
#'
#' print(myCI)
#'
#'
#' X1 = rnorm(n, sd = 4)
#' X2 = X1 + rnorm(n, sd = 0.1)
#' true_eps = rnorm(n)
#' Y = 3 + 8 * X1 + 4 * X2 + true_eps
#' X = cbind(X1, X2)
#'
#' myCI <- Navae_ci_ols(Y, X)
#'
#' print(myCI)
#'
#' @export
#'
Navae_ci_ols <- function(
    Y, X,
    alpha = 0.05,
    a = NULL, power_of_n_for_b = NULL,
    omega = NULL, power_of_n_for_omega = NULL,
    bounds = list(lambda_reg = NULL,
                  K_reg = NULL,
                  K_eps = NULL,
                  K_xi = NULL,
                  C = NULL,
                  B = NULL), K_xi = NULL,
    param_BE_EE = list(
      choice = "best",
      setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
      regularity = list(C0 = 1, p = 2),
      eps = 0.1),
    intercept = TRUE,
    options = list(center = FALSE, bounded_case = FALSE, with_Exp_regime = FALSE),
    matrix_u = NULL,
    verbose = 2)
{

  # 1- Checking the validity of inputs ==================================

  # Force X to be a matrix, even if there is only one variable
  # Same for matrix_u
  if (is.vector(X)) {X <- matrix(X, ncol = 1)}
  if (NCOL(matrix_u) == 1) {matrix_u <- matrix(matrix_u, ncol = 1)}

  if (!is.vector(Y)) {
    stop("Y must be a vector.")
  }
  if (length(Y) != nrow(X)) {
    stop("Y must have the same number of observations as X.")
  }
  if (is.null(matrix_u)){
    if (intercept){
      matrix_u = diag(NCOL(X) + 1)
    } else {
      matrix_u = diag(NCOL(X))
    }
  }
  if ( (ncol(matrix_u) != ncol(X) + 1) && intercept) {
    stop("matrix_u must have exactly one more column than X. ",
         "The first column of matrix_u is then interpreted ",
         "as the coefficient of the intercept. The other columns ",
         "of matrix_u correspond respectively to the columns of X.")
  }
  if (any(unlist(lapply(1:ncol(X), FUN = function(i) {length(unique(X[,i])) == 1})))) {
    stop("X must not contain any constant columns.")
  }
  if ((!is.logical(options$center)) || (length(options$center) != 1)) {
    stop("`options$center' must be either TRUE or FALSE.")
  }

  # 2- Computing fundamental quantities that will be useful later ========================

  number_u <- nrow(matrix_u)
  n <- nrow(X)

  # Add a column of ones and take the empirically recentered X if option center set to TRUE.
  if (intercept){
    X <- cbind(matrix(1, nrow = n, ncol = 1),
               scale(X, center = options$center, scale = FALSE))
    if (!is.null(colnames(X))) {
      colnames(X)[1] <- "intercept"
    } else {
      colnames(X) <- c("intercept", paste0("X", 1:(ncol(X) - 1)))
    }
  }
  p <- ncol(X)

  if (verbose >= 2){
    cat("Dimension of the problem: \n")
    cat("*  n = ", n, "\n")
    cat("*  p = ", p, "   ")
    if (intercept){
      cat("(including the intercept)\n")
    } else {
      cat("(model without intercept)\n")
    }
    cat("*  number of u = ", number_u, "\n")
    cat("\n")
  }

  # Estimation of crossproducts and other useful matrices
  XXt <- base::crossprod(X)
  XXtbar <- (1/n) * XXt # shortcut notation "S" in the article.

  # plug-in \hat{\lambda_{reg}} in the article.
  minvpXXtbar <- min(eigen(XXtbar, only.values = TRUE)$values)
  inverse_XXt <- solve(XXt)
  inverse_XXtbar <- n * inverse_XXt
  inverse_sqrt_XXtbar <- expm::sqrtm(inverse_XXtbar) # Estimate of E(XX')^{-1/2}

  # X_i tilde
  # this is a list of n elements which are all vectors of size p
  list_Xtilde_i <- lapply(
    X = 1:n,
    FUN = function(i) {inverse_sqrt_XXtbar %*% matrix(X[i, ], ncol = 1)})

  norms_row_X = apply(X = X, MARGIN = 1, FUN = function(u){sqrt(sum(u^2))})
  norms_row_X_tilde = unlist(lapply(X = list_Xtilde_i, FUN = function(u) {sqrt(sum(u^2))}))

  # Regression
  reg = stats::lm(Y ~ X - 1)
  betahat = reg$coefficients
  OLSestimate_u = matrix_u %*% matrix(betahat, ncol = 1)  # u' OLS estimate

  # Computation of Vhat
  Vhat = n * inverse_XXt %*% (crossprod(X * reg$residuals)) %*% inverse_XXt
  Vhat_u = apply(X = matrix_u, MARGIN = 1, FUN = function(u){t(u) %*% Vhat %*% u})

  # 3- Setting the parameters omega and a if not provided ======================

  # Parameter a = a_n
  if (is.null(power_of_n_for_b)) {
    power_of_n_for_b <- -2/5
  } else {
    stopifnot((is.numeric(power_of_n_for_b)) && length(power_of_n_for_b) == 1)
    if ((power_of_n_for_b > 0) || (power_of_n_for_b <= -1/2)) {
      warning("The choice of 'power_of_n_for_b' does not satisfy the ",
              "requirements for asymptotic properties of the resulting CI.")
    }
  }
  if (is.null(a)) {
    b_n <- n^power_of_n_for_b
    a <- 1 + b_n
  }

  # Parameter omega = omega_n
  if (is.null(power_of_n_for_omega)) {
    power_of_n_for_omega <- -1/5
  } else {
    stopifnot((is.numeric(power_of_n_for_omega)) && length(power_of_n_for_omega) == 1)
    if ((power_of_n_for_omega > 0) || (power_of_n_for_omega <= -2/3)) {
      warning("The choice of 'power_of_n_for_omega' does not satisfy the ",
              "requirements for asymptotic properties of the resulting CI.")
    }
  }
  if (is.null(omega)) {
    omega <- n^power_of_n_for_omega
    if (verbose >= 2){
      cat("omega = n^", power_of_n_for_omega, "\n")
    }
  }

  if (verbose >= 2){
    cat("Tuning parameters: \n")
    cat("*  a: "    , a    , "\n")
    cat("*  omega: ", omega, "\n\n")
  }


  # 4- Computing plug-ins for missing bounds ===================================

  env <- environment()
  allBounds = computeBounds(env = env, verbose = verbose)


  # 5- Computing concentration, Rlin, Rnvar ====================================

  # Computation of gamma
  gamma = alpha * omega / 2
  # just a name for the value at which we will compute Rnlin and Rnvar
  # in the article: omega_n alpha / 2.

  # Concentration of XXt
  # quantity called gamma_tilde as of now in the paper,
  # in the baseline case (when we do not assume bounded X).
  concentr_XXtranspose = Compute_concentrationXXt(
    bounded_case = options$bounded_case, bounds = bounds,
    delta = gamma, # only used for the bounded case, not present in the paper.
                   # to be checked later.
    n = n, d = ncol(X), gamma = gamma)

  Rnlin_u <- Compute_RnLin(
    gamma = gamma, gammatilde = concentr_XXtranspose,
    bounds = bounds, matrix_u = matrix_u)

  Rnvar_u <- Compute_Rnvar(
    gamma = gamma, n = n, norms_row_X = norms_row_X,
    residuals = reg$residuals, bounds = bounds,
    gammatilde = concentr_XXtranspose,
    X = X, inverse_XXtbar = inverse_XXtbar, matrix_u = matrix_u)

  Rnvar_u_times_norm_u_squared = Rnvar_u *
    apply(X = matrix_u, MARGIN = 1, FUN = function(x){sum(x^2)})


  if (verbose >= 2){
    cat("Concentration of XXt: \n")
    cat("*  gamma:"  ,  gamma  , "\n" )
    cat("*  Rnlin_u:",  Rnlin_u, "\n" )
    cat("*  Rnvar_u:",  Rnvar_u, "\n" )

    cat("\n")
  }


  # 6- Computing the standard asymptotic CIs ===================================
  # In order to compare (and used in the simulation)

  result_asymp = matrix(nrow = number_u, ncol = 3)
  colnames(result_asymp) = c("lower", "upper", "length")
  rownames(result_asymp) = colnames(X)
  result_asymp = as.data.frame(result_asymp)
  CIs.Asymp.extend = (stats::qnorm(1 - alpha/2) / sqrt(n)) * sqrt(Vhat_u)
  result_asymp[, 1] = OLSestimate_u - CIs.Asymp.extend
  result_asymp[, 2] = OLSestimate_u + CIs.Asymp.extend
  result_asymp[, 3] = 2 * CIs.Asymp.extend

  # 7- Preparing the final matrix with our CI ==================================
  # and first case of R regime (R1)

  result = matrix(nrow = number_u, ncol = 5)
  colnames(result) <- c("lower", "upper", "regime", "estimate", "length")
  rownames(result) <- colnames(X)
  result = as.data.frame(result)

  # u' OLS estimate (temporary, for check)
  result[, 4] = OLSestimate_u

  # First condition (independent of u) for returning \Rb
  # correspond to n > 2 K_reg / (omega_n alpha) in the paper to exit regime \Rb.
  if (concentr_XXtranspose >= 1) {

    result[, 1] = -Inf
    result[, 2] = Inf
    result[, 3] = "R1"
    result[, 5] = Inf
  }

  # 8- Computation of the bound delta_n ========================================

  param_BE_EE$setup$iid = TRUE # we always consider the iid framework

  if (allBounds["K_xi", "method"] == "provided by user") {
    # K_xi is provided => K_xi and delta_nE are uniform across vectors u

    delta_n_BE <- BE_bound_Shevtsova(bounds$K_xi, n)

    delta_n_EE <- BoundEdgeworth::Bound_EE1(
      setup = param_BE_EE$setup,
      regularity = param_BE_EE$regularity,
      eps = param_BE_EE$eps, n = n,
      K4 = bounds$K_xi, K3 = NULL, lambda3 = NULL, K3tilde = NULL)

    if (param_BE_EE$choice == "best") {
      if (delta_n_BE < delta_n_EE) {
        delta_n <- delta_n_BE; delta_n_from <- "BE"
      } else {
        delta_n <- delta_n_EE; delta_n_from <- "EE"
      }
    } else if (param_BE_EE$choice == "EE") {
      delta_n <- delta_n_EE; delta_n_from <- "EE"
    } else if (param_BE_EE$choice == "BE") {
      delta_n <- delta_n_BE; delta_n_from <- "BE"
    } else {
      stop("Invalid specification of the argument 'param_BE_EE$choice'.")
    }

    # Replicate the bound for each u (to have the same object with or wo plug-in)
    delta_n_u <- rep(delta_n, times = number_u)
    delta_n_from_u <- rep(delta_n_from, times = number_u)

  } else {

    delta_n_BE_u <- apply(
      X = array(1:number_u, dim = number_u), MARGIN = 1,
      FUN = function(index_u){BE_bound_Shevtsova(
        bound_K = bounds$K_xi_u[[index_u]], n = n)})

    delta_n_EE_u <- apply(
      X = array(1:number_u, dim = number_u), MARGIN = 1,
      FUN = function(index_u){BoundEdgeworth::Bound_EE1(
        setup = param_BE_EE$setup,
        regularity = param_BE_EE$regularity,
        eps = param_BE_EE$eps, n = n,
        K4 = bounds$K_xi_u[[index_u]], K3 = NULL, lambda3 = NULL, K3tilde = NULL)})

    # a priori, the ranking BE versus EE can be u-specific; hence the comparison
    # component by component when the choice is "best".
    if (param_BE_EE$choice == "best") {
      delta_n <- delta_n_from <- NULL # only for output with verbose = TRUE.
      delta_n_u <- vector(mode = "numeric", length = number_u)
      delta_n_from_u <- vector(mode = "character", length = number_u)
      for (index_u in seq_along(delta_n_u)) {
        if (delta_n_BE_u[[index_u]] < delta_n_EE_u[[index_u]]) {
          delta_n_u[[index_u]] <- delta_n_BE_u[[index_u]]
          delta_n_from_u[[index_u]] <- "BE"
        } else {
          delta_n_u[[index_u]] <- delta_n_EE_u[[index_u]]
          delta_n_from_u[[index_u]] <- "EE"
        }
      }
    } else if (param_BE_EE$choice == "EE") {
      delta_n_u <- delta_n_EE_u; delta_n_from <- "EE"
    } else if (param_BE_EE$choice == "BE") {
      delta_n_u <- delta_n_BE_u; delta_n_from <- "BE"
    } else {
      stop("Invalid specification of the argument 'param_BE_EE$choice'.")
    }

  }


  if (verbose >= 2){
    cat("Bound on Berry-Esseen / Edgeworth expansions: \n")
    cat("*  delta_n:"  ,  delta_n  , "\n" )
    cat("*  delta_n_u:",  delta_n_u  , "\n" )

    cat("\n")
  }

  # 9- Computing parts of our CI ===============================================

  nu_n_Exp_u = OLS.Nu_nExp(alpha = alpha, omega = omega, a = a,
                           K_xi = bounds$K_xi_u, n = n, verbose = verbose)

  nu_n_Edg_u = nu_n_Exp_u + delta_n_u

  bound_Voracle <- Vhat_u + Rnvar_u_times_norm_u_squared

  nu_n_Approx_u = Rnlin_u / sqrt(bound_Voracle)


  if (verbose >= 2){
    cat("Computing components of the CI: \n")
    cat("*  nu_n_Exp_u: "   ,  nu_n_Exp_u  , "\n" )
    cat("*  nu_n_Edg_u: "   ,  nu_n_Edg_u  , "\n" )
    cat("*  bound_Voracle: ",  bound_Voracle  , "\n" )
    cat("*  nu_n_Approx_u: ",  nu_n_Approx_u  , "\n" )

    cat("\n")
  }

  # 10- Determining the regime, with or without regime Exp, for each u =========

  if (options$with_Exp_regime) {

    which_regime_R <- which(nu_n_Exp_u >= alpha / 2)

    which_regime_Exp = which((nu_n_Exp_u < alpha / 2) & (nu_n_Edg_u >= alpha / 2))

    which_regime_Edg = which(nu_n_Edg_u < alpha / 2)
    # automatically implies nu_n_Exp_u < alpha / 2 too, since nu_n_Exp_u < nu_n_Edg_u

  } else {

    which_regime_Edg = which(nu_n_Edg_u < alpha / 2)

    which_regime_R <- setdiff((1:number_u), which_regime_Edg)

  }

  # 11- Constructing our CI ====================================================

  if (concentr_XXtranspose < 1) {
    # Otherwise, already in the case R1, R regime for all u.

    if (length(which_regime_R) > 0) {

      result[which_regime_R, 1] = -Inf
      result[which_regime_R, 2] = Inf
      result[which_regime_R, 3] = "R2"
      result[which_regime_R, 5] = Inf
    }

    if (length(which_regime_Edg) > 0) {

      CIs.Edg.extend = OLS.CIs.Edg.extend(
        n = n, alpha = alpha, a = a, K_reg = bounds$K_reg,
        nu_n_Edg = nu_n_Edg_u[which_regime_Edg],
        nu_n_Approx = nu_n_Approx_u[which_regime_Edg],
        bound_Voracle = bound_Voracle[which_regime_Edg])

      result[which_regime_Edg, 1] = OLSestimate_u[which_regime_Edg] - CIs.Edg.extend
      result[which_regime_Edg, 2] = OLSestimate_u[which_regime_Edg] + CIs.Edg.extend
      result[which_regime_Edg, 3] = "Edg"
      result[which_regime_Edg, 5] = 2 * CIs.Edg.extend
    }

    if (options$with_Exp_regime) {

      if (length(which_regime_Exp) > 0) {

        CIs.Exp.extend = OLS.CIs.Exp.extend(
          n = n, alpha = alpha, a = a, K_reg = bounds$K_reg,
          nu_n_Exp = nu_n_Exp_u[which_regime_Exp],
          nu_n_Approx = nu_n_Approx_u[which_regime_Exp],
          bound_Voracle = bound_Voracle[which_regime_Exp])

        result[which_regime_Exp, 1] = OLSestimate_u[which_regime_Exp] - CIs.Exp.extend
        result[which_regime_Exp, 2] = OLSestimate_u[which_regime_Exp] + CIs.Exp.extend
        result[which_regime_Exp, 3] = "Exp"
        result[which_regime_Exp, 5] = 2 * CIs.Exp.extend
      }
    }

  }

  # 12- Returning output ===============================================================

  about_delta_n = list(
    delta_n = delta_n,
    delta_n_u = delta_n_u,
    delta_n_from = delta_n_from,
    delta_n_from_u = delta_n_from_u)

  ratio_length_wrt_ci_asymp <- result[, 5] / result_asymp[, 3]

  minimal_alpha_to_enter_Edg_regime <- 2 * nu_n_Edg_u

  if (verbose >= 2){
    cat("Minimal alpha to enter Edg regime: ", minimal_alpha_to_enter_Edg_regime)
    cat("\n")
    cat("Ratio length wrt asymptotic CI: ", ratio_length_wrt_ci_asymp)
    cat("\n\n")
  }

  if (verbose >= 1){
    cat("Classical CI obtained from the CLT: \n")
    print(result_asymp)

    cat("\n NAVAE CI: \n")
    print(result)
    cat("\n")
  }

  return(list(ci_navae = result,
              ci_asymp = result_asymp,
              allBounds = allBounds,
              about_delta_n = about_delta_n,
              ratio_length_wrt_ci_asymp = ratio_length_wrt_ci_asymp,
              nu_n_Edg_u = nu_n_Edg_u,
              nu_n_Approx_u = nu_n_Approx_u,
              bound_Voracle = bound_Voracle,
              Rnlin_u = Rnlin_u,
              Rnvar_u = Rnvar_u,
              Rnvar_u_times_norm_u_squared = Rnvar_u_times_norm_u_squared,
              minimal_alpha_to_enter_Edg_regime = minimal_alpha_to_enter_Edg_regime,
              S = XXtbar))
}

# Auxiliary functions ----------------------------------------------------------

#' Auxiliary function used for the plug-in estimate of K_xi.
#' Computation of xi for a given observation i
#' and for each u, return a vector for the different u that are
#' the rows of matrix_u
#'
#' (Argument X named dataX to avoid confusion with the X argument of lapply.)
#'
#' @noRd
#'
Compute_xi_u_for_one_obs_i <- function(
    index_obs_i, dataX, inverse_XXtbar, matrix_u, residuals)
{
  return( matrix_u %*%
            (inverse_XXtbar %*%
               matrix(dataX[index_obs_i,], ncol = 1) * residuals[[index_obs_i]]) )
}


OLS.CIs.Exp.extend <- function(
    n, alpha, a, K_reg,
    nu_n_Exp, nu_n_Approx, bound_Voracle)
{

  part1_Qn = sqrt( 2 * (1 + a) * ( 1 - log(alpha/2 - nu_n_Exp) ) )

  QnExp_u = part1_Qn + nu_n_Approx

  result = ( QnExp_u / sqrt(n) ) * sqrt(bound_Voracle)

  return(result)
}

OLS.CIs.Edg.extend <- function(
    n, alpha, a, K_reg,
    nu_n_Edg, nu_n_Approx, bound_Voracle)
{

  part1_Qn = sqrt(a) * stats::qnorm(1 - alpha/2 + nu_n_Edg)

  QnEdg_u = part1_Qn + nu_n_Approx

  result = ( QnEdg_u / sqrt(n) ) * sqrt(bound_Voracle)

  return(result)
}

#' Computation of nu_nExp
#'
#' @noRd
OLS.Nu_nExp <- function(alpha, omega, a, K_xi, n, verbose)
{
  term_omega_alpha = omega * alpha / 2
  term_exp = exp( - n * (1 - 1/a)^2 / (2 * K_xi) ) / 2

  if (verbose >= 2){
    cat("Computing components of nu_n_Exp_u: \n")
    cat("*  term_omega_alpha: ", term_omega_alpha, "\n")
    cat("*  term_exp: ", term_exp, "\n\n")
  }

  result = term_omega_alpha + term_exp

  return (result)
}


#' Function to compute R_{n, lin}(gamma) for a free variable gamma
#' gammatilde corresponds to the sqrt(K_reg / (n * gamma))
#' in the paper.
#' It is possible to sharpen the concentration with additional assumptions
#' (bounded support, etc.), see function Compute_concentrationXXt.
#'
#' @noRd
Compute_RnLin <- function(
    gamma, gammatilde, bounds, matrix_u)
{
  RnLin_without_norm_u <-
    sqrt(2) * (1 / sqrt(bounds$lambda_reg)) * gammatilde / (1 - gammatilde)  *
    (bounds$K_eps / gamma)^(1/4)

  Rnlin_u <- apply(X = matrix_u, MARGIN = 1,
                   FUN = function(x){sqrt(sum(x^2))}) * RnLin_without_norm_u

  return(Rnlin_u)
}

#' Function to compute R_{n, var}(gamma) for a free variable gamma
#' (The different part follows the order of the paper.)
#' NB: in fact, it returns ||u||^2 * Rnvar(gamma), since Rnvar is always
#' used multiplied by ||u||^2 (equivalently, the part ||u||^2 could be put in
#' the definition of Rnvar as for Rnlin).
#'
#' @noRd
Compute_Rnvar <- function(
    gamma, n, norms_row_X, residuals,
    bounds,
    gammatilde,
    X, inverse_XXtbar, matrix_u)
{
  lambda_reg <- bounds$lambda_reg
  K_eps <- bounds$K_eps


  part1 = 2 / (n * lambda_reg^3) *
    (gammatilde / (1 - gammatilde) + 1)^2 *
    sqrt(K_eps / gamma) *
    mean(norms_row_X^4)

  part2 = 2 * sqrt(2) / (lambda_reg^(5/2) * sqrt(n)) *
    (gammatilde / (1 - gammatilde) + 1) *
    (K_eps / gamma)^(1/4) *
    mean(norms_row_X^3 * abs(residuals))

  part3 = (1 / lambda_reg^2) * (gammatilde / (1 - gammatilde))^2 *
    mean(norms_row_X^2 * residuals^2)

  part4 = 2 / lambda_reg * (gammatilde / (1 - gammatilde)) *
    base::norm(x = n^(-1) * crossprod(X * residuals) %*% inverse_XXtbar,
               type = "2")

  Rnvar = part1 + part2 + part3 + part4

  return(Rnvar)
}


Compute_concentrationXXt <- function(bounded_case, bounds, delta, n, d, gamma)
{
  if (bounded_case){
    # Concentration of XX transpose, assuming bounded regressors

    concentr_XXtranspose = Compute_concentration_Bernstein(
      B = bounds$B, C = bounds$C, delta = delta, n = n, d = d)

  } else {
    # Concentration of XX transpose without assuming bounded regressors

    concentr_XXtranspose <- sqrt(bounds$K_reg / (n * gamma))
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
Compute_concentration_Bernstein <- function(B, C, delta, n, d)
{
  return( sqrt(2 * B * log(2 * d / delta) / n) +
            4 * C * log(2 * d / delta) / (3*n) )
}



computeBounds <- function(env, verbose = 2)
{
  # Bound lambda_reg

  if (is.null(env$bounds$lambda_reg)) {

    env$bounds$lambda_reg <- env$minvpXXtbar
    allBounds = data.frame(name   = "lambda_reg",
                           value  = env$bounds$lambda_reg,
                           method = "plug-in")

  } else {
    allBounds = data.frame(name   = "lambda_reg",
                           value  = env$bounds$lambda_reg,
                           method = "provided by user")
  }


  # Bound K_reg

  if (is.null(env$bounds$K_reg)) {

    veca_to_mean_over_obs_i_for_K_reg <- unlist(lapply(
      X = 1:env$n, FUN = function(i) {
        sum((env$list_Xtilde_i[[i]] %*% t(env$list_Xtilde_i[[i]]) -
               diag(x = 1, nrow = env$p, ncol = env$p))^2)
      }))

    env$bounds$K_reg <- mean(veca_to_mean_over_obs_i_for_K_reg)

    allBounds[2, ] = data.frame(name   = "K_reg",
                                value  = env$bounds$K_reg,
                                method = "plug-in")
  } else {

    allBounds[2, ] = data.frame(name   = "K_reg",
                                value  = env$bounds$K_reg,
                                method = "provided by user")
  }

  # Bound K_epsilon

  if (is.null(env$bounds$K_eps)) {

    env$bounds$K_eps <- mean(env$norms_row_X_tilde^4 * env$reg$residuals^4)
    allBounds[3, ] = data.frame(name   = "K_eps",
                                value  = env$bounds$K_eps,
                                method = "plug-in")

  } else {
    allBounds[3, ] = data.frame(name   = "K_eps",
                                value  = env$bounds$K_eps,
                                method = "provided by user")
  }

  # Bound K_xi

  if (!is.null(env$bounds$K_xi) && !is.null(env$K_xi)) {
    stop("The bound K_xi was provided twice; use either K_xi or bounds argument")
  }

  if (is.null(env$bounds$K_xi) && is.null(env$K_xi)) {

    xi_u_i <- do.call(
      what = cbind,
      args = lapply(X = 1:env$n,
                    FUN = Compute_xi_u_for_one_obs_i,
                    dataX = env$X, inverse_XXtbar = env$inverse_XXtbar,
                    matrix_u = env$matrix_u, residuals = env$reg$residuals))

    mean_xi4_u <- rowMeans(xi_u_i^4)
    mean_xi2_u <- rowMeans(xi_u_i^2)
    empirical_kurtosis_xi_u <- mean_xi4_u / mean_xi2_u^2
    env$bounds$K_xi_u <- empirical_kurtosis_xi_u

    allBounds[4, ] = list(name   = "K_xi",
                          value  = list(list(env$bounds$K_xi_u)),
                          method = "plug-in")

    # use of name 'K_xi_u' to stress it is a vector (u-specific element)
    # and keep a scalar bounds$K_xi if given (as memory)

  } else {
    if (is.null(env$bounds$K_xi)) {
      env$bounds$K_xi <- env$K_xi
    }

    # If a bound on K_xi is provided, simply replicate it number_u times
    # to have a vector in bounds$K_xi_u as in the plug-in case.
    env$bounds$K_xi_u = rep(env$bounds$K_xi, length.out = env$number_u)

    allBounds[4, ] = list(name   = "K_xi",
                          value  = list(list(env$bounds$K_xi_u)),
                          method = "provided by user")
  }


  # Bound B and C (used when options$bounded_case is TRUE)
  # (used for Bernstein concentration of square matrices
  # applied to A = X_i tilde X_i tilde')

  if (is.null(env$bounds$C)) {

    env$bounds$C <- max(env$norms_row_X_tilde^2)
    allBounds[5, ] = data.frame(name   = "C",
                                value  = env$bounds$C,
                                method = "plug-in")
  } else {

    allBounds[5, ] = data.frame(name   = "C",
                                value  = env$bounds$C,
                                method = "provided by user")
  }

  if (is.null(env$bounds$B)) {

    # By definition of X_i tilde,
    # E[A] = E[X_i tilde X_i tilde'] = identity matrix of size p
    expectation_A = diag(env$p)
    B_before_norm = 0
    for (i in 1:env$n) {
      A_i = matrix(env$list_Xtilde_i[[i]]) %*% t(matrix(env$list_Xtilde_i[[i]]))
      # Computation of (A - E(A))(A - E(A))
      A_mEA_sq = (A_i - expectation_A) %*% (A_i - expectation_A)
      B_before_norm = B_before_norm + A_mEA_sq
    }

    env$bounds$B <- base::norm(x = B_before_norm, type = "2")
    allBounds[6, ] = data.frame(name   = "B",
                                value  = env$bounds$B,
                                method = "plug-in")
  } else {

    allBounds[6, ] = data.frame(name   = "B",
                                value  = env$bounds$B,
                                method = "provided by user")
  }

  rownames(allBounds) <- c("lambda_reg", "K_reg", "K_eps", "K_xi", "C", "B")

  if (verbose >= 2){
    cat("Bounds: \n")
    print(allBounds, row.names = FALSE)
    cat("\n")
  }

  return(allBounds)
}


