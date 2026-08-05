test_that("taxonomic Hill numbers hit known values", {
  # even S-species assemblage -> effective number = S at every q
  for (S in c(2, 4, 7)) {
    p <- rep(1, S)
    expect_equal(hill_taxonomic(p, 0), S)
    expect_equal(hill_taxonomic(p, 1), S)
    expect_equal(hill_taxonomic(p, 2), S)
  }
  # q = 0 is richness regardless of evenness
  expect_equal(hill_taxonomic(c(0.9, 0.09, 0.01), 0), 3)
  # ordering q0 >= q1 >= q2 for an uneven assemblage
  p <- c(0.7, 0.2, 0.1)
  expect_gte(hill_taxonomic(p, 0), hill_taxonomic(p, 1))
  expect_gte(hill_taxonomic(p, 1), hill_taxonomic(p, 2))
})

test_that("functional Hill is 1 for a single species and >1 for distinct spp", {
  D <- rbind(a = c(0, 1, 2), b = c(1, 0, 1.5), c = c(2, 1.5, 0))
  colnames(D) <- rownames(D)
  expect_equal(hill_functional(c(a = 1), D, 1), 1)
  expect_gt(hill_functional(c(a = 1, b = 1, c = 1), D, 1), 1)
})

test_that("phylogenetic Hill equals richness for a star tree with even abundance", {
  # star tree, equal branch lengths, even abundances -> effective lineages = S
  tr <- ape::read.tree(text = "(a:1,b:1,c:1,d:1);")
  p <- c(a = 1, b = 1, c = 1, d = 1)
  expect_equal(hill_phylo(p, tr, 0), 4)
  expect_equal(hill_phylo(p, tr, 1), 4, tolerance = 1e-8)
})

test_that("phylogenetic Hill is bounded by taxonomic richness", {
  tr <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  p <- c(a = 0.4, b = 0.3, c = 0.2, d = 0.1)
  expect_lte(hill_phylo(p, tr, 0), hill_taxonomic(p, 0) + 1e-8)
})
