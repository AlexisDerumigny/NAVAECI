



# This is in fact Rnvar multiplied by the square of the norm of u.
Rnvar_u <- apply(X = matrix_u, MARGIN = 1, 
                 FUN = function(x){sum(x^2)}) * Compute_Rnvar(
                   gamma = delta, n = n, norms_row_X = norms_row_X,
                   residuals = reg$residuals, bounds = bounds, X = X,
                   concentr_XXtranspose = concentr_XXtranspose, minvpXXtbar = minvpXXtbar,
                   inverse_XXtbar = inverse_XXtbar)

### TEMPORARY, for different checks and comparisons

Rnvar_v3 <- Compute_Rnvar(
  gamma = delta, n = n, norms_row_X = norms_row_X,
  residuals = reg$residuals, bounds = bounds, X = X,
  concentr_XXtranspose = concentr_XXtranspose, minvpXXtbar = minvpXXtbar,
  inverse_XXtbar = inverse_XXtbar, version = "v3")

Rnvar_v7 <- Compute_Rnvar(
  gamma = delta, n = n, norms_row_X = norms_row_X,
  residuals = reg$residuals, bounds = bounds, X = X,
  concentr_XXtranspose = concentr_XXtranspose,  minvpXXtbar = minvpXXtbar,
  inverse_XXtbar = inverse_XXtbar, version = "v7")

Rnvar_v3_optimized <- Compute_Rnvar(
  gamma = delta, n = n, norms_row_X = norms_row_X,
  residuals = reg$residuals, bounds = bounds, X = X,
  concentr_XXtranspose = concentr_XXtranspose,  minvpXXtbar = minvpXXtbar,
  inverse_XXtbar = inverse_XXtbar, version = "v3_optimized")


### TEMPORARY, for check different constructions
cc <- concentr_XXtranspose
lambda_reg <- bounds$lambda_reg

modif_part1_Rnvar_v4 <- (cc / (1 - cc) + 1)^2

modif_part1_Rnvar_v7 <- (lambda_reg * (1 - cc)^2)^(-1)

### TEMPORARY test du 4e et dernier terme de Rnvar
part4_Rnvar_v4 <- 2 / (lambda_reg * (1 - cc)^2) * cc *
  mean(norms_row_X^2 * (reg$residuals)^2)

part4_Rnvar_v7 = 2 * cc / (lambda_reg * (1 - cc)) *
  base::norm(x = n^(-1) * crossprod(X * reg$residuals) %*% inverse_XXtbar,
             type = "2")







