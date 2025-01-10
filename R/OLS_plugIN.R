
#' Decides whether or not to use plug-ins
#' and update the bounds if necessary
#'
#' @param env the environment of CI.OLS
#' in which the computations should be done
#'
#'
#' @noRd
#'
OLS.updateBounds <- function(env)
{

  if (is.null(env$bounds$lambda_m)){
    env$bounds$lambda_m <- env$minvpXXtbar
  }

  if (is.null(env$bounds$K_reg)){

    p <- ncol(env$X)

    veca_To_mean_over_obs_i_for_KX <- lapply(
      X = 1:env$n,
      FUN = function(i){
        sum( ( env$list_Xtilde_i[[i]] %*% t( env$list_Xtilde_i[[i]] ) -
                 diag(x = 1, nrow = p, ncol = p) )^2 )
      }
    )

    env$bounds$K_reg <- mean(veca_To_mean_over_obs_i_for_KX)
  }

  if (is.null(env$bounds$K_eps)){
    env$bounds$K_eps <- mean( env$norms_row_X^4 * env$reg$residuals^4 )
  }

  if (is.null(env$bounds$K_xi)){
    xi_u_i <- do.call(
      what = cbind,
      args = lapply(X = 1:env$n,
                    FUN = Compute_xi_u_for_one_obs_i,
                    dataX = env$X,
                    inverse_XXtbar = env$inverse_XXtbar,
                    matrix_u = env$matrix_u,
                    residuals = env$reg$residuals) )

    mean_xi4_u <- rowMeans(xi_u_i^4)
    mean_xi2_u <- rowMeans(xi_u_i^2)
    empirical_kurtosis_xi_u <- mean_xi4_u / mean_xi2_u^2
    env$bounds$K_xi_u <- empirical_kurtosis_xi_u
    # use of name 'K_xi_u' to stress it is a vector (u-specific element)
    # and keep a scalar bounds$K_xi if given (as memory)
  } else {
    # If a bound on K_xi is provided, simply replicate it number_u times
    # to have again a vector in bounds$K_xi_u
    env$bounds$K_xi_u = rep(env$bounds$K_xi,
                            length.out = env$number_u)
  }

  if (!is.null(env$bounds$K3_xi)){
    env$bounds$K3_xi_u = rep(env$bounds$K3_xi,
                             length.out = env$number_u)
  }

  if (!is.null(env$bounds$lambda3_xi)){
    env$bounds$lambda3_xi_u = rep(env$bounds$lambda3_xi,
                                  length.out = env$number_u)
  }

  if (!is.null(env$bounds$K3tilde_xi)){
    env$bounds$K3tilde_xi_u = rep(env$bounds$K3tilde_xi,
                                  length.out = env$number_u)
  }

  # Bound on C
  # (used for Bernstein concentration of square matrices, applied to A = X_i tilde X_i tilde')
  if (is.null(env$bounds$C)){
    env$bounds$C = max(env$norms_row_X_tilde^2)
  }

  # Bound on B
  # (used for Bernstein concentration of square matrices, applied to A = X_i tilde X_i tilde')
  if (is.null(env$bounds$B)){

    d = ncol(env$X)

    # By definition of X_i tilde,
    # E[A] = E[X_i tilde X_i tilde'] = identity matrix of size d
    expectation_A = diag(d)
    B_before_norm = 0

    for (i in 1:env$n){
      A_i = matrix(env$list_Xtilde_i[[i]]) %*%
        t(matrix(env$list_Xtilde_i[[i]]))

      # Computation of (A - E(A))(A - E(A))
      A_mEA_sq = (A_i - expectation_A) %*% (A_i - expectation_A)

      B_before_norm = B_before_norm + A_mEA_sq
    }
    env$bounds$B = base::norm(x = B_before_norm, type = "2")
  }
}


## Other functions used for plug-ins ===============================


#' Function used for the plug-in of K_xi
#' Computation of xi for a given observation i
#' and for each u, return a vector for the different u that are
#' the rows of matrix_u
#'
#' @noRd
#'
Compute_xi_u_for_one_obs_i <- function(index_obs_i,
                                       dataX, inverse_XXtbar, matrix_u, residuals){
  return( matrix_u %*%
            (inverse_XXtbar %*%
               matrix(dataX[index_obs_i,], ncol = 1) * residuals[[index_obs_i]]) )
}

