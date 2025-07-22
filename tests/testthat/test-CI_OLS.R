test_that("CI_OLS accept different choices for 'a'", {
  n = 4000
  X1 = rnorm(n, sd = 1)
  true_eps = rnorm(n)
  Y = 8 * X1 + true_eps
  X = cbind(X1)

  myCI <- Navae_ci_ols(Y, X, K_xi = 3, intercept = TRUE)
  myCI <- Navae_ci_ols(Y, X, K_xi = 3, intercept = TRUE, a = NULL)

  expect_identical(myCI$allTuningParameters$a$power_of_n_for_b$method,
                   expected = "default value")

  myCI_1 <- Navae_ci_ols(Y, X, K_xi = 3, intercept = TRUE,
                         a = 1 + 100 * n^(-1/5))

  expect_identical(myCI_1$allTuningParameters$a$method,
                   expected = "provided by user")

  myCI_2 <- Navae_ci_ols(Y, X, K_xi = 3, intercept = TRUE,
                         a = list(power_of_n_for_b = -1/5))

  expect_identical(myCI_2$allTuningParameters$a$method,
                   expected = "computed from 'power_of_n_for_b'")

  expect_identical(as.data.frame(myCI_1), as.data.frame(myCI_2))

})


test_that("CI_OLS accept different choices for 'omega'", {
  n = 4000
  X1 = rnorm(n, sd = 1)
  true_eps = rnorm(n)
  Y = 8 * X1 + true_eps
  X = cbind(X1)

  myCI <- Navae_ci_ols(Y, X, K_xi = 3, intercept = TRUE)
  myCI <- Navae_ci_ols(Y, X, K_xi = 3, intercept = TRUE, omega = NULL)

  expect_identical(myCI$allTuningParameters$omega$power_of_n_for_omega$method,
                   expected = "default value")

  myCI_1 <- Navae_ci_ols(Y, X, K_xi = 3, intercept = TRUE,
                         omega = n^(-1/5))

  expect_identical(myCI_1$allTuningParameters$omega$method,
                   expected = "provided by user")

  myCI_2 <- Navae_ci_ols(Y, X, K_xi = 3, intercept = TRUE,
                         omega = list(power_of_n_for_omega = -1/5))

  expect_identical(myCI_2$allTuningParameters$omega$power_of_n_for_omega$method,
                   expected = "provided by user")

  expect_identical(as.data.frame(myCI_1), as.data.frame(myCI_2))

  expect_warning({ myCI <- Navae_ci_ols(Y, X, K_xi = 3, intercept = TRUE,
                                        omega = list(power_of_n_for_omega = 2)) })
})

