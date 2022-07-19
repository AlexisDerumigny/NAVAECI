
#' Updates the tuning parameters
#' if they are not given by the user of CI.OLS
#'
#' @noRd
#'
OLS.updateTuningParameters <- function(env, bounded_case = FALSE)
{

  # Choice for omega
  if (is.null(env$omega)){
    # To have asymptotic pointwise exactness, we need:
    # omega_n = 1 / n^(p) with 0 < p < 2/3
    # default choice: p = 1/2, omega_n = 1 / sqrt(n) ?
    # Trade-off
    # Want omega_n large (p small) to have
    # n omega > 2 * bounds$K_X / (alpha * bounds$lambda_m^2 )
    # and avoid R1 regime
    # Want omega_n small (p large) to have nu_nExp small < alpha/2
    # nu_nExp is increasing in omega
    # and avoid R2 regime

    if (bounded_case){

      default_power_omega = 2/3 # doit etre plus grand que 1/2
      env$omega = 1 / env$n^default_power_omega

    } else {
      default_power_omega = 1/4
      env$omega = 1 / env$n^default_power_omega
    }

  }

  # # Choice for omega to be large enough (provided < 1) to avoid R1 regime
  # if (env$omega == "highest_possible_to_avoid_R1"){
  #
  #   add_to_be_strictly_above = 0.0001
  #
  #   threshold = 2 * env$bounds$K_X / (alpha * env$bounds$lambda_m^2 * n)
  #
  #   if (threshold < 1 - add_to_be_strictly_above){
  #     env$omega = threshold + add_to_be_strictly_above
  #   } else {
  #     env$omega = 1 - add_to_be_strictly_above
  #   }
  #
  # }

  # Choice for a
  if (is.null(env$a)){
    # To have asymptotic pointwise exactness, we need:
    # a_n = 1 + b_n = 1 + 1 / n^(p) with 0 < p < 1/2
    # Trade-off:
    # a_n large i.e. p small: more availability of CI (condition nu_n < alpha/2)
    # since nu_nExp (hence nu_nEdg) is decreasing in a.
    # a_n small i.e. p large: better precision
    # since Q_nExp and Q_nEdg are increasing in a.

    # idem dans le cas borne

    default_power_b_n = 1/4
    env$a = 1 + 1 / env$n^default_power_b_n

    }

}


