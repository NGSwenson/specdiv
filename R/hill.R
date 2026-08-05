# Hill numbers for the three diversity axes ----------------------------------
# All return an "effective number of species", so the three axes are directly
# comparable at a given order q. q = 0 richness, q = 1 Shannon, q = 2 Simpson.

#' Taxonomic Hill number
#'
#' Hill number of order `q` from a relative-abundance vector (abundances only).
#'
#' @param p Numeric vector of relative abundances (need not sum to 1; zeros are
#'   dropped).
#' @param q Order (`0` richness, `1` Shannon, `2` Simpson; any `q >= 0`).
#' @return Scalar effective number of species.
#' @examples
#' hill_taxonomic(c(0.5, 0.3, 0.2), q = 1)
#' @export
hill_taxonomic <- function(p, q) {
  p <- p[p > 0]; p <- p / sum(p)
  if (length(p) == 0) return(0)
  if (q == 1) exp(-sum(p * log(p))) else sum(p^q)^(1 / (1 - q))
}

# Internal: precompute per-edge descendant tips, scaled edge lengths, height.
.phylo_parts <- function(tree) {
  if (!inherits(tree, "phylo")) stop("`tree` must be an ape 'phylo' object.")
  Ntip <- length(tree$tip.label)
  edge <- tree$edge; EL <- tree$edge.length
  kids <- split(edge[, 2], edge[, 1])
  desc <- function(nd)
    if (nd <= Ntip) tree$tip.label[nd]
    else unlist(lapply(kids[[as.character(nd)]], desc))
  edge_tips <- lapply(edge[, 2], desc)
  H <- max(ape::node.depth.edgelength(tree)[seq_len(Ntip)])
  list(edge_tips = edge_tips, EL = EL, H = H, tips = tree$tip.label)
}

# Internal: phylogenetic Hill from precomputed parts and a tip-named p.
.hill_phylo_parts <- function(pn, q, parts) {
  a <- vapply(parts$edge_tips, function(tp) sum(pn[tp], na.rm = TRUE), numeric(1))
  k <- a > 0
  L <- parts$EL[k] / parts$H
  a <- a[k]
  if (q == 1) exp(-sum(L * a * log(a))) else sum(L * a^q)^(1 / (1 - q))
}

#' Phylogenetic Hill number (branch-based)
#'
#' Branch-based phylogenetic Hill number of Chao, Chiu & Jost (2010) on an
#' ultrametric tree. Returns an effective number of equally-abundant, maximally
#' distinct lineages.
#'
#' @param p Numeric vector of relative abundances, **named by tree tip label**.
#' @param tree An ape `phylo` object (ultrametric). Tip labels must cover
#'   `names(p)`.
#' @param q Order (`0`, `1`, `2`, ...).
#' @return Scalar effective number of lineages.
#' @references Chao A., Chiu C.-H. & Jost L. (2010) Phylogenetic diversity
#'   measures based on Hill numbers. Phil. Trans. R. Soc. B 365, 3599-3609.
#' @examples
#' tr <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
#' hill_phylo(c(a = 0.5, b = 0.3, c = 0.2), tr, q = 1)
#' @export
hill_phylo <- function(p, tree, q) {
  p <- p[p > 0]; p <- p / sum(p)
  if (length(p) == 0) return(0)
  .hill_phylo_parts(p, q, .phylo_parts(tree))
}

#' Functional Hill number (distance-based)
#'
#' Distance-based functional Hill number of Chiu & Chao (2014). Returns an
#' effective number of equally-abundant, functionally distinct species. For a
#' single species it is 1; note that at `q = 1` with exactly two species the
#' measure depends only on the abundances (the distance cancels), which is why
#' assemblage structure is better probed with [ses_structure()] (MPD/MNTD).
#'
#' @param p Numeric vector of relative abundances, named by species label.
#' @param dmat Symmetric functional distance matrix with dimnames covering
#'   `names(p)`.
#' @param q Order (`0`, `1`, `2`, ...).
#' @return Scalar effective number of functionally distinct species.
#' @references Chiu C.-H. & Chao A. (2014) Distance-based functional diversity
#'   measures and their decomposition. PLoS ONE 9, e100014.
#' @examples
#' D <- rbind(a = c(0, 1, 2), b = c(1, 0, 1.5), c = c(2, 1.5, 0))
#' colnames(D) <- rownames(D)
#' hill_functional(c(a = 0.5, b = 0.3, c = 0.2), D, q = 1)
#' @export
hill_functional <- function(p, dmat, q) {
  p <- p[p > 0]; p <- p / sum(p)
  if (length(p) < 2) return(1)
  dd <- dmat[names(p), names(p)]
  Q <- sum(outer(p, p) * dd)
  if (Q == 0) return(1)
  P2 <- outer(p, p)
  if (q == 1) exp(-0.5 * sum((dd / Q) * P2 * log(P2)))
  else sum((dd / Q) * P2^q)^(1 / (2 * (1 - q)))
}
