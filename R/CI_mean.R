#' Compute NAVAE CI for the expectation based on empirical mean estimator
#' and Berry-Esseen (BE) or Edgeworth Expansions (EE) bounds
#' 
#' @param data vector of univariate observations
#' 
#' @param alpha 1 - level of confidence of the CI.
#' By default, alpha is set to 0.05, yielding a 95% CI.
#' 
#' @param a the free parameter \eqn{a} of the interval.
#' If the argument is not provided, the following default choice is made
#' a = 1 + \eqn{b_n} with \eq{b_n = n^{-1/5}}.
#' Such a choice satisfies the conditions on the sequence of free parameters
#' for the CI to be asymptotically exact pointwise and uniformly over the set
#' of distributions with a bounded kurtosis.
#' 
#' @param bound_K bound on the kurtosis K_4(theta) of the distribution of the
#' observations that are assumed to be i.i.d.
#' The choice of 9 covers most "usual" distribution.
#' Otherwise, if the argument is not provided, the value used is the plug-in
#' counterpart \widehat{K}, that is, the empirical kurtosis of the observations.
#' 
#' @param param_BE_EE parameters to compute the BE or EE bound \delta_n used
#' to construct the confidence interval.
#' If \code{param_BE_EE} is exactly equal to \code{"BE"}, then the bound used is 
#' the best up-to-date BE bound from Shevtsova (2013) combined with a convexity 
#' inequality.
#' Otherwise, \code{param_BE_EE} is a list of four objects:
#' 1- \code{choice}: 
#' If equal to "EE", the bound used is Derumigny et al. (2023)'s bound
#' computed using the parameters specified by the rest of \code{param_BE_EE},
#' namely:
#' 2- \code{setup}: itself a logical vector of size 3,
#' 3- \code{regularity}: itself a list of length up to 3,
#' 4- \code{eps}: value between 0 and 1/3.
#' as described in the arguments of the function \code{BoundEdegworth::Bound_EE1}.
#' Together, they specify the bounds and assumptions used to compute the 
#' bound \eqn{\delta_n} from Derumigny et al. (2023).
#' Finally, if choice is equal to "best", the bound used is the minimum between
#' the previous one (with choice = "EE") and the bound "BE".
#' By default, following Remark 3.3 of the article, "best" is used and 
#' Derumigny et al. (2023)'s bounds is computed assuming i.i.d data and no other 
#' regularity assumptions (continuous or unskewed distribution) and the bound on
#' kurtosis used is the one specified in the previous the argument \code{bound_K}.
#' 
#' @return this function returns a numerical vector of size 2: the first element
#' is the lower end of the CI, the second the upper end.
#' In the R regime, the CI is \code{(-Inf, Inf)}.
#' 
#' @example 
#' n = 1000
#' data = rexp(n, 1)
#' Navae_ci_mean(data, bound_K = 9)
#' Navae_ci_mean(data) # plug-in for K
#' 
#' n = 50 * 10^3
#' data = rexp(n, 1)
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.2, a = 1 + n^(-0/5))
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.2)
#' Navae_ci_mean(data, alpha = 0.2) # plug-in for K
#' Navae_ci_mean(data, alpha = 0.1, param_BE_EE = list(choice = "best", setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE), regularity = list(C0 = 1, p = 2), eps = 0.1))
#' Navae_ci_mean(data, alpha = 0.1, param_BE_EE = list(choice = "best", setup = list(continuity = TRUE, iid = TRUE, no_skewness = FALSE), regularity = list(kappa = 0.99), eps = 0.1))
#' Navae_ci_mean(data, alpha = 0.05, param_BE_EE = list(choice = "best", setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE), regularity = list(C0 = 1, p = 2), eps = 0.1))
#' Navae_ci_mean(data, alpha = 0.05, param_BE_EE = list(choice = "best", setup = list(continuity = TRUE, iid = TRUE, no_skewness = FALSE), regularity = list(kappa = 0.99), eps = 0.1))
#' 
#' n = 100 * 10^3
#' data = rexp(n, 1)
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.10)
#' Navae_ci_mean(data, alpha = 0.10)
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.05)
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.05, a = 1 + n^(-1/3))
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.05, a = 1 + n^(-1/5))
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.05, a = 1 + n^(-0.5/5))
#' 
#' TODO: possibly change the choice of a when not provided to select, if exist
#' the best (for the resulting length of the IC) i.e. lowest a such that 
#' we are not in the R regime?
#' But in all cases, it requires quite large sample sizes. 
#' 
Navae_ci_mean <- function(
    data,
    alpha = 0.05, 
    a = NULL, 
    bound_K = NULL,
    param_BE_EE = list(
      choice = "best",
      setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
      regularity = list(C0 = 1, p = 2),
      eps = 0.1))
{
  
  stopifnot(is.vector(data, mode = "numeric"))
  xi <- data; # shortcut and to follow the notation of the paper
  
  n <- length(xi)
  xi_bar <- mean(xi); sigma_hat <- sqrt(var(xi))
  if (is.null(a)) {b_n <- n^(-1/5); a <- 1 + b_n}
  
  if (is.null(bound_K)) {
    bound_K <- Empirical_kurtosis(xi, xi_bar, sigma_hat)
  }
  
  if (is.vector(param_BE_EE, mode = "character") && (param_BE_EE == "BE")) {
    
    delta_n <- BE_bound_Shevtsova(bound_K, n)
    
  } else {
    
    delta_n_EE <- BoundEdgeworth::Bound_EE1(
      setup = param_BE_EE$setup, 
      regularity = param_BE_EE$regularity, 
      eps = param_BE_EE$eps, n = n, 
      K4 = bound_K, K3 = NULL, lambda3 = NULL, K3tilde = NULL)
    
    if (param_BE_EE$choice == "best") {
      delta_n <- min(delta_n_EE, BE_bound_Shevtsova(bound_K, n))
    } else if (param_BE_EE$choice == "EE") {
      delta_n <- delta_n_EE
    } else {
      stop("Invalid specification of the argument 'param_BE_EE$choice'.")
    }
  }
  
  nu_n_var <- exp(-n*(1 - 1/a)^2 / (2*bound_K))
  
  arg_modif_quant <- 1 - alpha/2 + delta_n + nu_n_var/2
  
  I_R_regime <- arg_modif_quant >= pnorm(sqrt(n / a))
  
  if (I_R_regime) {
    
    return(c(-Inf, Inf))
    
  } else {
    
    q <- qnorm(arg_modif_quant)
    C_n <- 1 / sqrt(1/a - (1/n) * q^2)
    half_length <- sigma_hat/sqrt(n) * C_n * q
    
    return(xi_bar + c(-1, 1) * half_length)
  }
}

# Auxiliary functions -----------------------------------------------------

Empirical_kurtosis <- function(xi, xi_bar, sigma_hat) {
  return(mean((xi - xi_bar)^4 / sigma_hat^4))
}

BE_bound_Shevtsova <- function(bound_K, n) {
  constant_BE <- 0.4690
  return(constant_BE * bound_K^(3/4) / sqrt(n))
}
