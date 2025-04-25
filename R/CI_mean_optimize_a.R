



n = 100
bound_K = 9

f <- function(x){
  a <- x
  nu_n_var <- exp(-n*(1 - 1/a)^2 / (2*bound_K))
  return(nu_n_var)
}

curve(f(x), from = 1, to = 2)

curve(qnorm(f(x)), from = 1, to = 2)

curve(qnorm(x), from = 0, to = 1)


f_to_zero_cond_R_reg <- function(a){
  nu_n_var <- exp(-n*(1 - 1/a)^2 / (2*bound_K))
  return(1 - alpha/2 + delta_n + nu_n_var/2 - pnorm(sqrt(n / a)))
}

f_nu_n_var <- function(a){
  return(exp(-n*(1 - 1/a)^2 / (2*bound_K)))
}

# Choice K, n, alpha
bound_K = 9
n = 700
alpha = 0.2

delta_n <- BE_bound_Shevtsova(bound_K = bound_K, n = n)

curve(f_to_zero_cond_R_reg(a), from = 1, to = 200, xname = "a")

res_root <- uniroot(f = f_to_zero_cond_R_reg, lower = 1, upper = 10)


f_prop_width <- function(a){
  nu_n_var <- f_nu_n_var(a)
  arg_modif_quant <- 1 - alpha/2 + delta_n + nu_n_var/2
  q <- stats::qnorm(arg_modif_quant)
  C_n <- 1 / sqrt(1/a - (1/n) * q^2)
  return(C_n * q)
}

res_opt <- optimise(f = f_prop_width, lower = res_root$root, upper = 10)

2.6 / qnorm(1 - alpha/2)
res_opt$objective / qnorm(1 - alpha/2)

curve(f_prop_width(a), from = 1, to = 3, xname = "a")


Get_delta_n_fixed_bound_K <- function(
    n, bound_K,
    param_BE_EE = list(
      choice = "best",
      setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
      regularity = list(C0 = 1, p = 2),
      eps = 0.1))
{
  
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
  
  return(delta_n)
}


Get_optimal_non_stochastic_a <- function(
    n, bound_K, alpha,
    param_BE_EE = list(
      choice = "best",
      setup = list(continuity = FALSE, iid = TRUE, no_skewness = FALSE),
      regularity = list(C0 = 1, p = 2),
      eps = 0.1), max_a_tested = 10)
{
  
  delta_n <- Get_delta_n_fixed_bound_K(n, bound_K, param_BE_EE)
  
  # Step 1: find the region to exit R regime
  
  f_to_zero_cond_R_reg <- function(a) {
    nu_n_var <- exp(-n*(1 - 1/a)^2 / (2*bound_K))
    return(1 - alpha/2 + delta_n + nu_n_var/2 - pnorm(sqrt(n / a)))
  }
  
  if (f_to_zero_cond_R_reg(max_a_tested) > 0) {
    warning("Cannot find a to exit R regime.")
    return(list(cannot_find_a_to_exit_R_regime = TRUE,
                minimal_a_to_exit_R_regime = NA,
                optimal_a_to_minimize_width = NA))
  }
  
  res_uniroot <- uniroot(f = f_to_zero_cond_R_reg, lower = 1, upper = max_a_tested)
  minimal_a_to_exit_R_regime <- res_uniroot$root
  
  # Step 2: minimize the width within the relevant region
  
  f_prop_width <- function(a) {
    nu_n_var <- exp(-n*(1 - 1/a)^2 / (2*bound_K))
    arg_modif_quant <- 1 - alpha/2 + delta_n + nu_n_var/2
    q <- stats::qnorm(arg_modif_quant)
    C_n <- 1 / sqrt(1/a - (1/n) * q^2)
    return(C_n * q)
  }
  
  res_opt <- optimise(f = f_prop_width, lower = minimal_a_to_exit_R_regime, upper = max_a_tested)
  optimal_a_to_minimize_width <- res_opt$minimum
  optimal_width_up_to_sigmahat_sqrtn <- res_opt$objective
  
  return(list(cannot_find_a_to_exit_R_regime = FALSE,
              minimal_a_to_exit_R_regime = minimal_a_to_exit_R_regime,
              optimal_a_to_minimize_width = optimal_a_to_minimize_width,
              optimal_width_up_to_sigmahat_sqrtn = optimal_width_up_to_sigmahat_sqrtn))
}

# Find minimal a to exit R regime and optimal a ----------------------------------------------------------

bound_K = 9
alpha = 0.10
# vec_n <- c((1:50), 100, 200, 300, 400, 500, 1000) * 10^3
vec_n <- c(seq(from = 50, to = 1000, by = 10)) * 10^3


res_a <- matrix(data = NA, nrow = length(vec_n), ncol = 4)
colnames(res_a) <- c(
  "n", "minimal_a_to_exit_R_regime", "optimal_a_to_minimize_width",
  "optimal_width_up_to_sigmahat_sqrtn")

for (i in seq_along(vec_n)) {
  
  res_a[i, 1] <- vec_n[[i]]
  
  res_one_n <- Get_optimal_non_stochastic_a(
    n = vec_n[[i]], bound_K = bound_K, alpha = alpha)
  
  if (isFALSE(res_one_n$cannot_find_a_to_exit_R_regime)) {
    res_a[i, 2] <- res_one_n$minimal_a_to_exit_R_regime
    res_a[i, 3] <- res_one_n$optimal_a_to_minimize_width
    res_a[i, 4] <- res_one_n$optimal_width_up_to_sigmahat_sqrtn
  }
    
  rm(res_one_n)
}

b_n <- res_a[,3] - 1

res_lm <- lm(I(log10(b_n)) ~ I(log10(vec_n)))
summary(res_lm)

plot(log10(vec_n), log10(b_n))
abline(a = res_lm$coefficients[[1]], b = res_lm$coefficients[[2]], col = "red")

plot(vec_n, res_a[,4] / qnorm(1 - alpha/2))


# Function g, to study the function ----------------------------------------------------------

n <- 500
alpha <- 0.5
bound_K <- 10

delta_n <- Get_delta_n_fixed_bound_K(n = n, bound_K = bound_K)

g <- function(a) 
{
  nu_n_var <- exp(-n*(1 - 1/a)^2 / (2*bound_K))
  return(1 - alpha/2 + delta_n + nu_n_var/2 - pnorm(sqrt(n / a)))
}

g_prime <- function(a)
{
  nu_n_var <- exp(-n*(1 - 1/a)^2 / (2*bound_K))
  part1 <- -n / (2*bound_K) * (1 - (1/a)) / a^2 * nu_n_var
  part2 <- sqrt(n) / 2 * (1 / (a * sqrt(a))) * dnorm(sqrt(n / a))
  return(part1 + part2)
}

g_prime_Alexis <- function(a)
{
  nu_n_var <- exp(-n*(1 - 1/a)^2 / (2*bound_K))
  return(n / (2*bound_K) * (a - 1) * nu_n_var - dnorm(sqrt(n / a)))
}

g_prime_numeric <- function(a, h = 0.000001)
{
  return((g(a + h) - g(a)) / h)
}

lim_at_plusinf_g <- 1/2 - alpha/2 + delta_n + exp(-n / (2*bound_K))/2
lim_at_one_g <- 1/2 - alpha/2 + delta_n

curve(g(a), from = 1, to = 1000, xname = "a")

g(1)
g(1 + 10^(-3))
g(1.4)
g(100)
lim_at_plusinf_g

curve(g_prime(a), from = 1, to = 100, xname = "a")

ylim_choice <- 0.01
curve(g_prime(a), from = 1, to = 200, xname = "a",
      ylim = c(-ylim_choice, ylim_choice))

eval_at <- 10.3
g_prime_numeric(eval_at)
g_prime(eval_at)
g_prime_Alexis(eval_at)

curve(g(a), from = 1, to = 100*n, xname = "a")
