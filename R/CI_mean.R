#' Compute NAVAE CI for the expectation based on empirical mean estimator
#' and Berry-Esseen (BE) or Edgeworth Expansions (EE) bounds
#'
#' @param data vector of univariate observations.
#'
#' @param alpha 1 - level of confidence of the CI; nominal level = 1 - alpha.
#' By default, alpha is set to 0.05, yielding a 95\% CI.
#'
#' @param a the free parameter \eqn{a} (or \eqn{a_n}) of the interval.
#' @param power_of_n_for_b
#' If the argument \code{a} is not provided, the following default choice is
#' made for the free parameter:
#' a_n = a = 1 + \eqn{b_n} with \eqn{b_n = n^{\code{power_of_n_for_b}}}.
#' If \code{power_of_n_for_b} is not provided (default argument NULL),
#' the default choice is \code{power_of_n_for_b} = -2/5.
#' Such a choice satisfies the conditions on the sequence of free parameters
#' for the CI to be asymptotically exact pointwise and uniformly over the set
#' of distributions with a bounded kurtosis. If the specified
#' \code{power_of_n_for_b} violates those condition, a warning is output.
#'
#' @param bound_K bound on the kurtosis K_4(theta) of the distribution of the
#' observations that are assumed to be i.i.d.
#' The choice of 9 covers most "usual" distributions.
#' If the argument is not provided (default argument NULL), the value used is
#' the plug-in counterpart \eqn{\widehat{K}}, that is, the empirical kurtosis
#' of the observations.
#'
#' @param known_variance by default NULL, in this case, the function computes
#' the CI in the general case with an unknown variance (which is estimated).
#' Otherwise, a scalar numeric vector equal to the (assumed/known) variance.
#' (NB: if the option is used, one must provide the variance and not the standard
#' deviation.)
#'
#' @param param_BE_EE parameters to compute the BE or EE bound \eqn{\delta_n} used
#' to construct the confidence interval.
#' If \code{param_BE_EE} is exactly equal to \code{"BE"}, then the bound used is
#' the best up-to-date BE bound from Shevtsova (2013) combined with a convexity
#' inequality.
#' Otherwise, \code{param_BE_EE} is a list of four objects: \itemize{
#'   \item \code{choice}:
#'   If equal to \code{"EE"}, the bound used is Derumigny et al. (2023)'s bound
#'   computed using the parameters specified by the rest of \code{param_BE_EE},
#'   namely
#'   \item \code{setup}: itself a logical vector of size 3,
#'   \item \code{regularity}: itself a list of length up to 3,
#'   \item \code{eps}: value between 0 and 1/3,
#' }
#' as described in the arguments of the function \code{BoundEdegworth::Bound_EE1}.
#' Together, they specify the bounds and assumptions used to compute the
#' bound \eqn{\delta_n} from Derumigny et al. (2023).
#' Finally, if \code{choice} is equal to \code{"best"}, the bound used is the minimum
#' between the previous one (with \code{choice = "EE"}) and the bound \code{"BE"}.
#'
#' By default, following Remark 3.3 of the article, \code{"best"} is used and
#' Derumigny et al. (2023)'s bounds is computed assuming i.i.d data and no other
#' regularity assumptions (continuous or unskewed distribution) and the bound on
#' kurtosis used is the one specified in the previous the argument \code{bound_K}.
#'
#' @param na.rm logical, should missing values in \code{data} be removed?
#'
#' @param verbose FALSE by default, if TRUE, several additional elements are reported
#'
#' @return with \code{verbose = FALSE}, it returns a numerical vector of size 2:
#' the first element is the lower end of the CI, the second the upper end.
#' In the \eqn{\mathbb{R}} regime, the CI is \code{(-Inf, Inf)}; otherwise
#' it is a numeric 2-length vector.
#' with \code{verbose = TRUE}, a list is output with various elements:
#' - the CI;
#' - the classical "asymptotic" CI based on CLT (as a comparison);
#' - the regime of our CI (indicator of being in \eqn{\mathbb{R}} regime)
#' - the bound delta_n used: its numerical value and whether it comes from BE or EE;
#' - the value of the argument of the modified quantile;
#' - the minimal alpha to exit the \eqn{\mathbb{R}} regime.
#' - the value K used and the method
#'
#' @examples
#' n = 1000
#' data = rexp(n, 1)
#' Navae_ci_mean(data, bound_K = 9)
#' Navae_ci_mean(data, bound_K = 9, a = 6)
#' Navae_ci_mean(data) # plug-in for K
#'
#' n = 1000
#' data = rexp(n, 1)
#' Navae_ci_mean(data, bound_K = 9)
#' Navae_ci_mean(data, bound_K = 9, a = 1 + 5)
#' Navae_ci_mean(data) # plug-in for K
#'
#' n = 5 * 10^3
#' data = rexp(n, 1)
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.2, a = 1 + n^(-0/5))
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.2)
#' Navae_ci_mean(data, alpha = 0.2) # plug-in for K
#'
#' listParams1 = list(
#'   choice = "best",
#'   setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
#'   regularity = list(C0 = 1, p = 2),
#'   eps = 0.1)
#'
#' listParams2 = list(
#'   choice = "best",
#'   setup = list(continuity = TRUE, iid = TRUE, no_skewness = FALSE),
#'   regularity = list(kappa = 0.99), eps = 0.1)
#'
#' Navae_ci_mean(data, alpha = 0.1, param_BE_EE = listParams1)
#' Navae_ci_mean(data, alpha = 0.1, param_BE_EE = listParams2)
#' Navae_ci_mean(data, alpha = 0.05, param_BE_EE = listParams1)
#' Navae_ci_mean(data, alpha = 0.05, param_BE_EE = listParams2)
#'
#' n = 1000 * 10^3
#' data = rexp(n, 1)
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.10)
#' Navae_ci_mean(data, alpha = 0.10)
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.05)
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.05, a = 1 + n^(-1/3))
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.05, a = 1 + 4)
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.05, a = 1 + n^(-1/5))
#' Navae_ci_mean(data, bound_K = 9, alpha = 0.05, a = 1 + n^(-0.5/5))
#'
#' @export
#
# TODO: possibly change the choice of a when not provided to select, if exist
# the best (for the resulting length of the IC) i.e. lowest a such that
# we are not in the R regime?
# But in all cases, it requires quite large sample sizes.
#
Navae_ci_mean <- function(
  data, alpha = 0.05,
  a = NULL, power_of_n_for_b = NULL, 
  optimize_in_a = FALSE, max_a_tested_for_optim = 50,
  bound_K = NULL, known_variance = NULL,
  param_BE_EE = list(
    choice = "best",
    setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
    regularity = list(C0 = 1, p = 2),
    eps = 0.1),
  na.rm = FALSE, verbose = FALSE)
{

  # 1- Checks on data and computation of empirical mean and variance -----------

  if (!is.vector(data, mode = "numeric")) {
    stop("'data' must be a numeric vector.")
  }

  if (na.rm) {
    xi <- stats::na.exclude(data)
  } else {
    xi <- data
  }
  n <- length(xi)

  if (n < 2) {
    stop("'data' must be at least of length 2.")
  }

  xi_bar <- mean(xi)
  sigma_hat <- sqrt(stats::var(xi))

  # 2- Determination of the bound K and the bound delta_n ----------------------
  
  if (is.null(bound_K)) {
    bound_K_method <- "plug-in"
    bound_K <- Empirical_kurtosis(xi = xi, xi_bar = xi_bar, sigma_hat = sigma_hat)
  } else {
    bound_K_method <- "bound"
  }
  
  delta_n_BE <- BE_bound_Shevtsova(bound_K = bound_K, n = n)
  
  if (is.vector(param_BE_EE, mode = "character") && (param_BE_EE == "BE")) {
    
    delta_n <- delta_n_BE; delta_n_from <- "BE"
    
  } else {
    
    delta_n_EE <- BoundEdgeworth::Bound_EE1(
      setup = param_BE_EE$setup,
      regularity = param_BE_EE$regularity,
      eps = param_BE_EE$eps, n = n,
      K4 = bound_K, K3 = NULL, lambda3 = NULL, K3tilde = NULL)
    
    if (param_BE_EE$choice == "best") {
      if (delta_n_BE < delta_n_EE) {
        delta_n <- delta_n_BE; delta_n_from <- "BE"
      } else {
        delta_n <- delta_n_EE; delta_n_from <- "EE"
      }
    } else if (param_BE_EE$choice == "EE") {
      delta_n <- delta_n_EE; delta_n_from <- "EE"
    } else {
      stop("Invalid specification of the argument 'param_BE_EE$choice'.")
    }
  }
  
  # 3- Determination of the free parameter a_n ---------------------------------

  if (optimize_in_a) {
    
    if (is.null(bound_K_method == "plug-in")) {
      warning("Optimization in a is performed while using plug-in.")
    }
    
    res_optimal_a <- Get_optimal_a(
      n = n, bound_K = bound_K, alpha = alpha, delta_n = delta_n,
      max_a_tested_for_optim = max_a_tested_for_optim)
    
    cannot_find_a_to_exit_R_regime <- res_optimal_a$cannot_find_a_to_exit_R_regime
    minimal_a_to_exit_R_regime <- res_optimal_a$minimal_a_to_exit_R_regime
    optimal_a_to_minimize_width <- res_optimal_a$optimal_a_to_minimize_width
    optimal_width_up_to_sigmahat_sqrtn <- res_optimal_a$optimal_width_up_to_sigmahat_sqrtn
    
    if (is.na(optimal_a_to_minimize_width)) {
      warning("Optimization in a failed; default choice for a is used.")
      optimized_a_is_used <- FALSE
    } else {
      optimized_a_is_used <- TRUE
      a <- optimal_a_to_minimize_width
      b_n <- a - 1
    }
    
  }
  
  if (isFALSE(optimize_in_a) || isFALSE(optimized_a_is_used)) {
    
    optimized_a_is_used <- FALSE
    cannot_find_a_to_exit_R_regime <- NA
    minimal_a_to_exit_R_regime <- NA
    optimal_a_to_minimize_width <- NA
    optimal_width_up_to_sigmahat_sqrtn <- NA
    
    if (is.null(power_of_n_for_b)) {
      power_of_n_for_b <- -2/5
    } else {
      stopifnot((is.numeric(power_of_n_for_b)) && length(power_of_n_for_b) == 1)
      if ((power_of_n_for_b > 0) || (power_of_n_for_b <= -1/2)) {
        warning("The choice of 'power_of_n_for_b' does not satisfy the requirements for asymptotic properties of the resulting CI.")
      }
    }
    
    if (is.null(a)) {
      b_n <- n^power_of_n_for_b
      a <- 1 + b_n
    }

  }
  
  # 4- Computation of our CI ---------------------------------------------------

  if (is.null(known_variance)) {

    # 4a- CI in the general case with unknown variance -------------------------

    nu_n_var <- Compute_nu_n_var(n = n, a = a, bound_K = bound_K)
    
    arg_modif_quant <- 1 - alpha/2 + delta_n + nu_n_var/2

    indicator_R_regime <- arg_modif_quant >= stats::pnorm(sqrt(n / a))

    minimal_alpha_to_exit_R_regime <-
      2*(1 + delta_n + nu_n_var/2 - stats::pnorm(sqrt(n / a)))

    if (indicator_R_regime) {

      ci <- c(-Inf, Inf)
      ratio_length_wrt_CLT <- Inf

    } else {

      q <- stats::qnorm(arg_modif_quant)
      C_n <- 1 / sqrt(1/a - (1/n) * q^2)
      half_length <- sigma_hat * C_n * q / sqrt(n)
      ci <- xi_bar + c(-1, 1) * half_length
      ratio_length_wrt_CLT <- C_n * q / stats::qnorm(1 - alpha/2)
    }

  } else {

    # 4b- CI in the particular case with known variance ------------------------

    if ((!is.numeric(known_variance)) || (length(known_variance) > 1)) {
      stop("Argument 'known_variance' must be a scalar numeric vector.")
    }

    arg_modif_quant <- 1 - alpha/2 + delta_n

    indicator_R_regime <- arg_modif_quant >= 1

    minimal_alpha_to_exit_R_regime <- 2 * delta_n

    if (indicator_R_regime) {

      ci <- c(-Inf, Inf)
      ratio_length_wrt_CLT <- Inf

    } else {

      q <- stats::qnorm(arg_modif_quant)
      half_length <- sqrt(known_variance) * q / sqrt(n)
      ci <- xi_bar + c(-1, 1) * half_length
      ratio_length_wrt_CLT <- q / stats::qnorm(1 - alpha/2)
    }

  }

  # 5- Output ------------------------------------------------------------------

  if (verbose) {

    # Comparison with usual "asymptotic" CI (based on CLT + Slutsky)
    if (is.null(known_variance)) {
      sigma_used_for_CLT <- sigma_hat
    } else {
      sigma_used_for_CLT <- sqrt(known_variance)
    }
    half_length_CLT <- sigma_used_for_CLT * stats::qnorm(1 - alpha/2) / sqrt(n)
    ci_CLT <- xi_bar + c(-1, 1) * half_length_CLT

    return(list(ci = ci,
                indicator_R_regime = indicator_R_regime,
                ci_CLT = ci_CLT,
                ratio_length_wrt_CLT = ratio_length_wrt_CLT,
                delta_n = delta_n,
                delta_n_from = delta_n_from,
                arg_modif_quant = arg_modif_quant,
                minimal_alpha_to_exit_R_regime = minimal_alpha_to_exit_R_regime,
                bound_K_method = bound_K_method,
                bound_K_value = bound_K,
                a = a,
                b_n = b_n,
                optimized_a_is_used = optimized_a_is_used,
                cannot_find_a_to_exit_R_regime = cannot_find_a_to_exit_R_regime,
                minimal_a_to_exit_R_regime = minimal_a_to_exit_R_regime,
                optimal_a_to_minimize_width = optimal_a_to_minimize_width,
                optimal_width_up_to_sigmahat_sqrtn = optimal_width_up_to_sigmahat_sqrtn))
    
  } else {
    return(ci)
  }
}

# Auxiliary functions ----------------------------------------------------------

Compute_nu_n_var <- function(n, a, bound_K)
{
  return(exp(-n*(1 - 1/a)^2 / (2*bound_K)))  
}

Empirical_kurtosis <- function(xi, xi_bar, sigma_hat)
{
  return(mean((xi - xi_bar)^4 / sigma_hat^4))
}

BE_bound_Shevtsova <- function(bound_K, n)
{
  constant_BE_iid <- 0.4690
  return(constant_BE_iid * bound_K^(3/4) / sqrt(n))
}

Get_optimal_a <- function(
  n, bound_K, alpha, delta_n,
  max_a_tested_for_optim = 50)
{
  
  # Step 1: find the region to exit R regime
  
  f_R_regime_if_non_negative <- function(a) {
    nu_n_var <- Compute_nu_n_var(n = n, a = a, bound_K = bound_K)
    return(1 - alpha/2 + delta_n + nu_n_var/2 - stats::pnorm(sqrt(n / a)))
  }
  
  # Check if the set of a exiting R regime is empty or not
  res_optim_condition_R_regime <- optimise(
    f = f_R_regime_if_non_negative, 
    lower = 1, upper = max_a_tested_for_optim)
  
  if (res_optim_condition_R_regime$objective >= 0) {
    warning("Optimization in a: cannot find a to exit R regime.")
    return(list(cannot_find_a_to_exit_R_regime = TRUE,
                minimal_a_to_exit_R_regime = NA,
                optimal_a_to_minimize_width = NA))
  }
  
  res_uniroot_condition_R_regime <- uniroot(
    f = f_R_regime_if_non_negative, 
    lower = 1, upper = res_optim_condition_R_regime$minimum)
  
  minimal_a_to_exit_R_regime <- res_uniroot_condition_R_regime$root
  
  # Step 2: minimize the width within the relevant set of a exiting R regime
  
  f_CI_width_up_to_sigmahat_sqrtn <- function(a) {
    nu_n_var <- Compute_nu_n_var(n = n, a = a, bound_K = bound_K)
    arg_modif_quant <- 1 - alpha/2 + delta_n + nu_n_var/2
    q <- stats::qnorm(arg_modif_quant)
    C_n <- 1 / sqrt(1/a - (1/n) * q^2)
    return(C_n * q)
  }
  
  res_optim_CI_width <- optimise(
    f = f_CI_width_up_to_sigmahat_sqrtn, 
    lower = minimal_a_to_exit_R_regime, upper = max_a_tested_for_optim)
  optimal_a_to_minimize_width <- res_optim_CI_width$minimum
  optimal_width_up_to_sigmahat_sqrtn <- res_optim_CI_width$objective
  
  # Check we find a relevant a, i.e., exiting R regime 
  if (f_R_regime_if_non_negative(optimal_a_to_minimize_width) >= 0) {
    warning("Optimization in a: failure, the a found does not exit the R regime.")
    return(list(cannot_find_a_to_exit_R_regime = FALSE,
                minimal_a_to_exit_R_regime = NA,
                optimal_a_to_minimize_width = NA))
  }
  
  return(list(cannot_find_a_to_exit_R_regime = FALSE,
              minimal_a_to_exit_R_regime = minimal_a_to_exit_R_regime,
              optimal_a_to_minimize_width = optimal_a_to_minimize_width,
              optimal_width_up_to_sigmahat_sqrtn = optimal_width_up_to_sigmahat_sqrtn))
}

