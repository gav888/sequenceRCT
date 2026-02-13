test_that("test_complexity_measures uses Kruskal-Wallis for more than two groups", {
  complexity_df <- data.frame(
    Entropy = c(0.20, 0.35, 0.45, 0.30, 0.55, 0.60, 0.42, 0.50, 0.62),
    Turbulence = c(1.2, 1.4, 1.5, 1.1, 1.6, 1.7, 1.3, 1.55, 1.8),
    Volatility = c(0.10, 0.25, 0.40, 0.15, 0.35, 0.45, 0.20, 0.30, 0.50),
    Group = factor(rep(c("Control", "TreatmentA", "TreatmentB"), each = 3))
  )

  results <- test_complexity_measures(complexity_df)

  expect_named(results, c("Entropy", "Turbulence", "Volatility"))
  expect_true(all(vapply(results, function(x) "kruskal" %in% names(x), logical(1))))
  expect_true(all(vapply(results, function(x) inherits(x$kruskal, "htest"), logical(1))))
})
