#' Compute NAVAE CI for the expectation based on empirical mean estimator
#' and Berry-Esseen (BE) or Edgeworth Expansions (EE) bounds
#'
#' @param data vector of univariate observations.
#'
#' @param alpha this is 1 minus the confidence level of the CI; in other words,
#' the nominal level is 1 - alpha.
#' By default, \code{alpha} is set to \code{0.05}, yielding a 95\% CI.
#'
#' @param a the free parameter \eqn{a} (or \eqn{a_n}) of the interval.
#' @param power_of_n_for_b
#' If the argument \code{a} is not provided, the following default choice is
#' made for the free parameter:
#' \code{a_n = a = 1 + b_n} with \code{b_n = n^{power_of_n_for_b}}.
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
#' If the argument is not provided (default argument \code{NULL}), the value used
#' is the plug-in counterpart \eqn{\widehat{K}}, that is, the empirical kurtosis
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
#' as described in the arguments of the function
#' \code{BoundEdgeworth::\link[BoundEdgeworth]{Bound_EE1}}.
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
#' @return A list is output with various elements:
#' - the CI;
#' - the classical "asymptotic" CI based on CLT (as a comparison);
#' - the regime of our CI (indicator of being in \eqn{\mathbb{R}} regime)
#' - the bound delta_n used: its numerical value and whether it comes from BE or EE;
#' - the value of the argument of the modified quantile;
#' - the minimal alpha to exit the \eqn{\mathbb{R}} regime.
#' - the value K used and the method
#'
#' @references
#' For the Edgeworth expansion bounds:
#'
#' Derumigny A., Girard L., and Guyonvarch Y. (2023).
#' Explicit non-asymptotic bounds for the distance to the first-order Edgeworth expansion,
#' Sankhya A. \doi{10.1007/s13171-023-00320-y}
#' \href{https://arxiv.org/abs/2101.05780}{arxiv:2101.05780}.
#'
#'
#' @seealso \code{\link{Navae_ci_ols}} the corresponding function for the linear
#' regression case.
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
  a = NULL, power_of_n_for_b = NULL, optimize_in_a = FALSE,
  bound_K = NULL, known_variance = NULL,
  param_BE_EE = list(
    choice = "best",
    setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
    regularity = list(C0 = 1, p = 2),
    eps = 0.1),
  na.rm = FALSE)
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
    stop("'data' have at least 2 non-missing observations.")
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

    properties_optimal_a <- Get_optimal_a(
      n = n, bound_K = bound_K, alpha = alpha, delta_n = delta_n)

    if (is.na(properties_optimal_a$optimal_a_to_minimize_width)) {
      warning("Optimization in a failed; default choice for a is used.")

      # In this case, we revert to a default choice of a since the optimization
      # does not work.
      optimized_a_is_used <- FALSE
      properties_optimal_a <- list(use_optimal_a = FALSE)
    } else {
      optimized_a_is_used <- TRUE
      a <- properties_optimal_a$optimal_a_to_minimize_width
      b_n <- a - 1
    }
  } else {
    optimized_a_is_used <- FALSE
    properties_optimal_a <- list(use_optimal_a = FALSE)
  }

  if (isFALSE(optimized_a_is_used)) {
    # Then we do the default choice of a

    if (is.null(power_of_n_for_b)) {
      power_of_n_for_b <- -2/5
    } else {
      if (!is.numeric(power_of_n_for_b) || (length(power_of_n_for_b) != 1)){
        stop("`power_of_n_for_b` must be a numeric vector of length 1. ",
             "Here it is: ", power_of_n_for_b)
      }
      if ((power_of_n_for_b > 0) || (power_of_n_for_b <= -1/2)) {
        warning("The choice of 'power_of_n_for_b' does not satisfy the ",
                "requirements for asymptotic properties of the resulting CI.")
      }
    }

    if (is.null(a)) {
      b_n <- n^power_of_n_for_b
      a <- 1 + b_n
    } else {
      # Define b_n also in this case for consistency since it is also reported
      # in the verbose output.
      b_n <- a - 1
    }
  }

  # 4- Computation of our CI ---------------------------------------------------

  if (is.null(known_variance)) {

    # 4a- CI in the general case with unknown variance -------------------------

    known_variance_used <- FALSE

    nu_n_var <- Compute_nu_n_var(n = n, a = a, bound_K = bound_K)

    arg_modif_quant <- 1 - alpha/2 + delta_n + nu_n_var/2

    indicator_R_regime <- arg_modif_quant >= stats::pnorm(sqrt(n / a))

    minimal_alpha_to_exit_R_regime <-
      2*(1 + delta_n + nu_n_var/2 - stats::pnorm(sqrt(n / a)))

    if (indicator_R_regime) {

      ci <- c(-Inf, Inf)
      ratio_length_wrt_CLT <- Inf
      C_n <- NA # put as NA for output in verbose, to be defined in all cases.

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

    known_variance_used <- TRUE

    arg_modif_quant <- 1 - alpha/2 + delta_n

    C_n <- NA # not defined in the case of known variance; put as NA for output in verbose.

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
  length_ci = ifelse(test = indicator_R_regime,
                     yes = Inf, no = ci[[2]] - ci[[1]])

  # 4c- Comparison with usual "asymptotic" CI (based on CLT + Slutsky) ---------

  if (is.null(known_variance)) {
    sigma_used_for_CLT <- sigma_hat
  } else {
    sigma_used_for_CLT <- sqrt(known_variance)
  }
  half_length_CLT <- sigma_used_for_CLT * stats::qnorm(1 - alpha/2) / sqrt(n)
  ci_CLT <- xi_bar + c(-1, 1) * half_length_CLT


  # 5- Output ------------------------------------------------------------------

  result = list(ci = ci,
                indicator_R_regime = indicator_R_regime,
                known_variance_used = known_variance_used,
                length_ci = length_ci,
                ci_CLT = ci_CLT,
                ratio_length_wrt_CLT = ratio_length_wrt_CLT,
                delta_n = delta_n,
                delta_n_from = delta_n_from,
                arg_modif_quant = arg_modif_quant,
                C_n = C_n,
                minimal_alpha_to_exit_R_regime = minimal_alpha_to_exit_R_regime,
                bound_K_method = bound_K_method,
                bound_K_value = bound_K,
                a = a,
                b_n = b_n,
                properties_optimal_a = properties_optimal_a,
                n = n,
                call = match.call())

  class(result) <- "NAVAE_CI_Mean"

  return (result)
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


#' Find the optimal a
#'
#'
#' @param n sample size
#' @param bound_K bound on the kurtosis
#' @param alpha confidence level
#' @param delta_n value of the bound $delta_n$ from Berry-Esseen or
#' first-order Edgeworth expansion.
#'
#' @examples
#' # Does not find any value of a
#' Get_optimal_a(n = 10000, bound_K = 9, alpha = 0.05,
#'               delta_n = BoundEdgeworth::Bound_EE1(n = 10000) )
#'
#' # Find something
#' Get_optimal_a(n = 100000, bound_K = 9, alpha = 0.05,
#'               delta_n = BoundEdgeworth::Bound_EE1(n = 100000) )
#'
#'
#' @noRd
Get_optimal_a <- function(n, bound_K, alpha, delta_n)
{

  # Step 1: find the region to exit R regime

  f_R_regime_if_non_negative <- function(a) {
    nu_n_var <- Compute_nu_n_var(n = n, a = a, bound_K = bound_K)
    return(1 - alpha/2 + delta_n + nu_n_var/2 - stats::pnorm(sqrt(n / a)))
  }

  # Check if the set of a exiting R regime is empty or not
  res_optim_condition_R_regime <- stats::optim(
    f = f_R_regime_if_non_negative,
    lower = 1, upper = Inf,
    method = "L-BFGS-B",
    par = 1 + n^(-1/5) # Initial value of the optimizer
  )

  if (res_optim_condition_R_regime$convergence > 0) {
    warning("Optimization in a: cannot find a to exit R regime. ",
            res_optim_condition_R_regime$message)

    return(list(cannot_find_a_to_exit_R_regime = TRUE,
                minimal_a_to_exit_R_regime = NA,
                optimal_a_to_minimize_width = NA,
                res_optim_condition_R_regime = res_optim_condition_R_regime))
  }

  if (res_optim_condition_R_regime$value > 0) {
    warning("Optimization in a: cannot find a to exit R regime.")

    return(list(cannot_find_a_to_exit_R_regime = TRUE,
                minimal_a_to_exit_R_regime = NA,
                optimal_a_to_minimize_width = NA,
                res_optim_condition_R_regime = res_optim_condition_R_regime))
  }

  # We try to find a_1, which is the smallest value of a such that we enter the
  # R regime. This is the zero of the function f_R_regime_if_non_negative that
  # is at the left of the current minimum.
  res_uniroot_a_1 <- stats::uniroot(
    f = f_R_regime_if_non_negative,
    lower = 1, upper = res_optim_condition_R_regime$par)

  a_1 <- res_uniroot_a_1$root

  try_find_a_2 <- TRUE
  max_a_tested_for_optim = 100 * a_1

  while(try_find_a_2){
    res_uniroot_a_2  <- tryCatch(
      stats::uniroot(f = f_R_regime_if_non_negative,
              # We do not start at 1 if not we may find the other root a_1 again (!)
              lower = res_optim_condition_R_regime$par,
              upper = max_a_tested_for_optim),
      error = function(e) e
    )
    if (inherits(res_uniroot_a_2, "error")){
      max_a_tested_for_optim = 1000 * max_a_tested_for_optim
      if (max_a_tested_for_optim == Inf){
        try_find_a_2 = FALSE
        a_2 = Inf
      }
    } else {
      a_2 = res_uniroot_a_2$root
      try_find_a_2 = FALSE
    }
  }

  # Step 2: minimize the width within the relevant set of a exiting R regime

  f_CI_width_up_to_sigmahat_sqrtn <- function(a) {
    nu_n_var <- Compute_nu_n_var(n = n, a = a, bound_K = bound_K)
    arg_modif_quant <- 1 - alpha/2 + delta_n + nu_n_var/2
    # cat("a = ", a, "\n")
    # cat("arg_modif_quant = ", arg_modif_quant, "\n\n")
    q <- stats::qnorm(arg_modif_quant)
    C_n <- 1 / sqrt(1/a - (1/n) * q^2)
    return(C_n * q)
  }

  res_optim_CI_width <- stats::optimize(
    f = f_CI_width_up_to_sigmahat_sqrtn,
    lower = a_1, upper = a_2)

  optimal_a_to_minimize_width <- res_optim_CI_width$minimum
  optimal_width_up_to_sigmahat_sqrtn <- res_optim_CI_width$objective

  return(list(cannot_find_a_to_exit_R_regime = FALSE,
              minimal_a_to_exit_R_regime = a_1,
              maximal_a_to_exit_R_regime = a_2,
              optimal_a_to_minimize_width = optimal_a_to_minimize_width,
              optimal_width_up_to_sigmahat_sqrtn = optimal_width_up_to_sigmahat_sqrtn))
}

