
test_that("Navae_ci_mean works with different kinds of 'a' inputs", {
  n = 10000

  x = rexp(n, 1)

  result1 = Navae_ci_mean(x, bound_K = 9)

  result2 = Navae_ci_mean(x, bound_K = 9, a = 6)

  expect_warning(Navae_ci_mean(x, bound_K = 9, a = 0.2))

  expect_error(Navae_ci_mean(x, bound_K = 9, a = list("")))

  expect_error(Navae_ci_mean(x, bound_K = 9, a = list()))

  expect_warning(Navae_ci_mean(x, bound_K = 9, a = list(power_of_n_for_b = 1) ) )

  result3 <- suppressWarnings({Navae_ci_mean(x, bound_K = 9,
                                             a = list(power_of_n_for_b = 1) )})

  allMethods_a = c(result1$properties_a$method,
                   result2$properties_a$method,
                   result3$properties_a$method)

  # All these methods should be different
  expect_true(length(unique(allMethods_a)) == 3)


  # Two ways to say the same thing:
  result4 <- Navae_ci_mean(x, bound_K = 9, alpha = 0.2, a = 1 + n^(-2/5))

  result5 <- Navae_ci_mean(x, bound_K = 9, alpha = 0.2, a = list(power_of_n_for_b = -2/5))

  expect_identical(result4$a, result5$a)
})

test_that("Navae_ci_mean works in the known variance case", {
  n = 10000

  x = rexp(n, 1)

  result1 = Navae_ci_mean(x, bound_K = 9, known_variance = 1)

  result2 = Navae_ci_mean(x, bound_K = 9, known_variance = 1, a = 6)

  result3 = Navae_ci_mean(x, bound_K = 9, known_variance = 2)

  expect_identical(result1$ci_navae, result2$ci_navae)

  expect_true(! identical(result1$ci_navae , result3$ci_navae) )
})

