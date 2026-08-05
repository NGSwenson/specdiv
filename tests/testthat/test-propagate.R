make_toy <- function(seed = 1) {
  set.seed(seed)
  sp <- c("a", "b", "c", "d")
  tr <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  D <- trait_distance(rbind(a = c(0, 0), b = c(0, 1), c = c(3, 0), d = c(3, 1)))
  n <- 20; site <- rep(1:2, each = 10)
  P <- matrix(0.02, n, 4, dimnames = list(NULL, sp))
  P[1:10, "a"] <- 0.94                     # site 1: confident single species
  P[11:20, ] <- 0.25                       # site 2: uniform
  P <- P / rowSums(P)
  truth <- c(rep("a", 10), rep(c("a", "b", "c", "d"), length.out = 10))
  list(P = P, site = site, tree = tr, D = D, truth = truth, sp = sp)
}

test_that("propagate_diversity returns the expected shape and columns", {
  z <- make_toy()
  out <- propagate_diversity(z$P, z$site, tree = z$tree, trait_dist = z$D,
                             truth = z$truth, q = c(0, 1, 2),
                             n_draws = 100, seed = 1)
  expect_equal(nrow(out), 2)
  for (m in c("TD1", "PD1", "FD1"))
    for (suf in c("", "_lo", "_hi", ""))
      expect_true(paste0(m, suf) %in% names(out))
  expect_true(all(c("true_TD1", "true_PD1", "true_FD1") %in% names(out)))
})

test_that("a confident single-species site is less diverse than a uniform site", {
  z <- make_toy()
  out <- propagate_diversity(z$P, z$site, tree = z$tree, trait_dist = z$D,
                             n_draws = 200, seed = 1)
  expect_lt(out$TD1[1], out$TD1[2])
  expect_lt(out$PD1[1], out$PD1[2])
})

test_that("intervals bracket the mean and q-ordering holds", {
  z <- make_toy()
  out <- propagate_diversity(z$P, z$site, tree = z$tree, trait_dist = z$D,
                             n_draws = 200, seed = 1)
  expect_true(all(out$TD1_lo <= out$TD1 + 1e-9))
  expect_true(all(out$TD1 <= out$TD1_hi + 1e-9))
  expect_true(all(out$TD0 >= out$TD1 - 1e-6))
  expect_true(all(out$TD1 >= out$TD2 - 1e-6))
})

test_that("taxonomic-only works with no tree or trait_dist", {
  z <- make_toy()
  out <- propagate_diversity(z$P, z$site, n_draws = 50, seed = 1)
  expect_true(all(grepl("^TD|^site|^n_seedlings", names(out))))
  expect_false(any(grepl("^PD|^FD", names(out))))
})

test_that("mismatched species/tree/dist are caught", {
  z <- make_toy()
  bad_tree <- ape::read.tree(text = "((a:1,b:1):1,x:2);")
  expect_error(propagate_diversity(z$P, z$site, tree = bad_tree), "tip")
  expect_error(propagate_diversity(z$P, z$site[1:5]), "length")
})

test_that("diversity_from_labels agrees with a hand count", {
  lab <- c("a", "a", "b", "b", "c", "a")
  site <- c(1, 1, 1, 2, 2, 2)
  d <- diversity_from_labels(lab, site)
  expect_equal(d$TD0, c(2, 3))               # site1 {a,b}; site2 {b,c,a}
})
