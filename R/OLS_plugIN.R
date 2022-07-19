
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

  if (is.null(env$bounds$K_X)){
    env$bounds$K_X <- mean( purrr::map_dbl(1:env$n,
                                           To_mean_over_observation_i_for_K_X,
                                           X = env$X, XXtbar = env$XXtbar) )
  }

  if (is.null(env$bounds$K_eps)){
    env$bounds$K_eps <- mean( env$norms_row_X^4 * env$reg$residuals^4 )
  }

  if (is.null(env$bounds$K_xi)){
    xi_u_i <- purrr::reduce( (purrr::map(1:env$n,
                                         Compute_xi_u_for_one_obs_i,
                                         X = env$X,
                                         inverse_XXtbar = env$inverse_XXtbar,
                                         matrix_u = env$matrix_u,
                                         residuals = env$reg$residuals)
                              ), cbind)
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

}


## Other functions used for plug-ins ===============================

#' Function used for the plug-in of K_X
#'
#' @noRd
#'
To_mean_over_observation_i_for_K_X <- function(index_obs_i,
                                               X, XXtbar){
  return( sum( c( ( matrix(X[index_obs_i,], ncol = 1) %*%
                      t(matrix(X[index_obs_i,], ncol = 1)) ) - XXtbar )^2 ) )
}

#' Function used for the plug-in of K_xi
#' Computation of xi for a given observation i
#' and for each u, return a vector for the different u that are
#' the rows of matrix_u
#'
#' @noRd
#'
Compute_xi_u_for_one_obs_i <- function(index_obs_i,
                                       X, inverse_XXtbar, matrix_u, residuals){
  return( matrix_u %*%
            (inverse_XXtbar %*%
               matrix(X[index_obs_i,], ncol = 1) * residuals[[index_obs_i]]) )
}

# Essai pour un KX normalise avec des bornes relatives - commente pour le moment
# To_mean_KX_normalized_on_i <- function(i, X, inverse_XXtbar){
#   sum( c(
#     inverse_XXtbar %*%
#       matrix(X[i,], ncol = 1) %*% t(matrix(X[i,], ncol = 1)) -
#       diag(x = 1,
#            nrow = NROW(inverse_XXtbar), ncol = NCOL(inverse_XXtbar)) )^2)
# }
