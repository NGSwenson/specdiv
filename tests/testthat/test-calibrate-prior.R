test_that("scores_to_prob rows are non-negative and sum to 1", {
  s <- matrix(c(0.2, -0.1, 0.9, 0.3, -5, -5), 3, 2, byrow = TRUE)
  P <- scores_to_prob(s)
  expect_true(all(P >= 0))
  expect_equal(unname(rowSums(P)), rep(1, 3))
  # all-nonpositive row becomes uniform
  expect_equal(unname(P[3, ]), c(0.5, 0.5))
})

test_that("temperature_scale sharpens and preserves the argmax", {
  P <- matrix(c(0.6, 0.4, 0.55, 0.45), 2, 2, byrow = TRUE)
  Pg <- temperature_scale(P, 4)
  expect_equal(unname(rowSums(Pg)), rep(1, 2))
  expect_true(all(apply(Pg, 1, max) >= apply(P, 1, max)))     # sharper
  expect_equal(max.col(Pg), max.col(P))                       # same winner
  expect_equal(temperature_scale(P, 1), P)                    # gamma=1 no-op
})

test_that("fit_temperature matches mean confidence toward accuracy", {
  set.seed(1)
  P <- matrix(runif(300), 100, 3); P <- P / rowSums(P)
  colnames(P) <- c("a", "b", "c")
  y <- colnames(P)[max.col(P)]                # perfectly accurate labels
  ft <- fit_temperature(P, y)
  expect_true(ft$gamma >= 1)
  expect_equal(ft$cv_accuracy, 1)
})

test_that("apply_abundance_prior tilts toward common species and renormalises", {
  P <- matrix(0.5, 2, 2); colnames(P) <- c("rare", "common")
  post <- apply_abundance_prior(P, c(common = 90, rare = 10))
  expect_equal(unname(rowSums(post)), rep(1, 2))
  expect_true(all(post[, "common"] > post[, "rare"]))
  # unknown prior names are ignored, not an error
  expect_message(apply_abundance_prior(P, c(common = 5, zzz = 1)), "Ignoring")
})

test_that("batch_align sets new-batch column means to the reference means", {
  set.seed(1)
  ref <- matrix(rnorm(200), 40, 5)
  new <- matrix(rnorm(30) + 3, 6, 5)          # offset batch
  aligned <- batch_align(new, ref)
  expect_equal(colMeans(aligned), colMeans(ref), tolerance = 1e-8)
})
