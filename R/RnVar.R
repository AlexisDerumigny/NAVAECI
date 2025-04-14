#'
#'
#'
#'
Compute_Rnvar <- function(
  gamma, n, norms_row_X, residuals,
  bounds,
  concentr_XXtranspose,
  inverse_XXtbar,
  minvpXXtbar, X,
  version = "v7")
{
  cc <- concentr_XXtranspose
  lambda_reg <- bounds$lambda_reg
  K_eps <- bounds$K_eps
  K_reg <- bounds$K_reg
  
  # v3 of overleaf for the first three terms and fourth one from v7
  # a priori, the smallest/most optimized one for practical use.
  
  part1_v3 <- 2 / (n * lambda_reg^2 * minvpXXtbar^2) *
    sqrt(K_eps / gamma) * mean(norms_row_X^4)
  
  part2_v3 <- 2 * sqrt(2) / (lambda_reg^2 * minvpXXtbar * sqrt(n)) *
    (K_eps / gamma)^(1/4) * mean(norms_row_X^3 * abs(residuals))
  
  part3_v3 <- K_reg / (n * gamma * lambda_reg^2 * minvpXXtbar^2) *
    mean(norms_row_X^2 * residuals^2)
  
  part4_v3 <- 2 / (lambda_reg * minvpXXtbar^2) * sqrt(K_eps / gamma) *
    (1 / sqrt(n)) * mean(norms_row_X^2 * residuals^2)
  
  # current version v7 (less optimized, but cleaner conceptually in the sense
  # that we use the bound lambda_reg instead of lambda_min(XX') empirique)
  
  part1_v7 = (2 / (n * lambda_reg^3 * (1 - cc)^2)) *
    sqrt(K_eps / gamma) * mean(norms_row_X^4)
  
  part2_v7 = (2 * sqrt(2) / (lambda_reg^(5/2) * sqrt(n) * (1 - cc))) *
    (K_eps / gamma)^(1/4) * mean(norms_row_X^3 * abs(residuals))
  
  part3_v7 = (cc / (lambda_reg * (1 - cc)))^2 *
    mean(norms_row_X^2 * residuals^2)
  
  part4_v7 = 2 * cc / (lambda_reg * (1 - cc)) *
    base::norm(x = n^(-1) * crossprod(X * residuals) %*% inverse_XXtbar,
               type = "2")
  
  if (version == "v7") {
    
    return(part1_v7 + part2_v7 + part3_v7 + part4_v7)
    
  } else if (version == "v3_optimized") {
    
    return(part1_v3 + part2_v3 + part3_v3 + part4_v7)
    
  } else if (version == "v3") {
  
    return(part1_v3 + part2_v3 + part3_v3 + part4_v3)
    
  } else {
    
    stop("Invalid argument 'version'.")
  }
}
