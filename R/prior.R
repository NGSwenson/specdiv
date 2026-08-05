# Abundance prior for the open-species-pool setting --------------------------

#' Apply a Bayesian abundance prior to a probability matrix
#'
#' Multiply each individual's likelihood (a row of `P`) by an independent
#' estimate of species frequencies and renormalise, giving a posterior. This
#' recovers most of the accuracy lost when the candidate library is broad rather
#' than closed to the species known present. The prior can come from a cheap,
#' non-destructive source such as counts from a second set of sites that were
#' never measured with the instrument.
#'
#' A small pseudocount `alpha` is added so no candidate is hard-zeroed. Because
#' any abundance prior down-weights uncommon species, it helps the realistic
#' open-pool case but taxes rare species when the pool is already correctly
#' closed.
#'
#' @param P Probability matrix (rows = individuals, columns = species).
#' @param prior_counts Named numeric vector of species counts (names = species
#'   labels). Species absent from the names get only the pseudocount. Extra
#'   names not in `colnames(P)` are ignored with a message.
#' @param alpha Pseudocount added to every candidate before normalising
#'   (default 0.5).
#' @return A posterior probability matrix the same shape as `P`, dimnames kept.
#' @examples
#' P <- matrix(c(0.5, 0.5, 0.5, 0.5), 2, 2, byrow = TRUE)
#' colnames(P) <- c("rare", "common")
#' apply_abundance_prior(P, c(common = 90, rare = 10))
#' @export
apply_abundance_prior <- function(P, prior_counts, alpha = 0.5) {
  P <- as.matrix(P)
  cls <- colnames(P)
  if (is.null(cls)) stop("`P` must have species column names.")
  extra <- setdiff(names(prior_counts), cls)
  if (length(extra))
    message(sprintf("Ignoring %d prior name(s) not among P's species: %s",
                    length(extra), paste(extra, collapse = ", ")))
  cnt <- prior_counts[cls]
  cnt[is.na(cnt)] <- 0
  pv <- as.numeric(cnt) + alpha
  pv <- pv / sum(pv)
  post <- sweep(P, 2, pv, "*")
  out <- post / rowSums(post)
  dimnames(out) <- dimnames(P)
  out
}

#' Functional distance matrix from species traits
#'
#' Standardise each trait to unit variance and return the Euclidean distance
#' matrix among species (or individuals). The row names of `traits` become the
#' dimnames, so they can be matched to a probability matrix's species columns.
#'
#' @param traits Numeric matrix or data frame, rows = species (row names =
#'   labels), columns = traits.
#' @param scale Logical; standardise each trait to unit variance first
#'   (default `TRUE`).
#' @return A symmetric numeric distance matrix with species dimnames.
#' @examples
#' tr <- rbind(a = c(1, 2), b = c(3, 1), c = c(2, 5))
#' trait_distance(tr)
#' @export
trait_distance <- function(traits, scale = TRUE) {
  M <- as.matrix(traits)
  if (is.null(rownames(M))) stop("`traits` must have species row names.")
  Z <- if (scale) base::scale(M) else M
  D <- as.matrix(stats::dist(Z))
  dimnames(D) <- list(rownames(M), rownames(M))
  D
}
