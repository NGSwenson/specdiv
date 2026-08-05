# Monte Carlo propagation of classification uncertainty into diversity --------

# Internal: build a function that maps a named count vector to all requested
# axis x order Hill numbers. Phylogeny parts are precomputed once.
.make_axis_fn <- function(species, q, tree, trait_dist) {
  parts <- if (!is.null(tree)) .phylo_parts(tree) else NULL
  function(cnt) {
    cnt <- cnt[cnt > 0]
    p <- cnt / sum(cnt)
    out <- numeric(0)
    for (qq in q) out[paste0("TD", qq)] <- hill_taxonomic(p, qq)
    if (!is.null(parts))
      for (qq in q) out[paste0("PD", qq)] <-
        .hill_phylo_parts(p, qq, parts)
    if (!is.null(trait_dist))
      for (qq in q) out[paste0("FD", qq)] <- hill_functional(p, trait_dist, qq)
    out
  }
}

.axis_names <- function(q, tree, trait_dist) {
  nm <- paste0("TD", q)
  if (!is.null(tree)) nm <- c(nm, paste0("PD", q))
  if (!is.null(trait_dist)) nm <- c(nm, paste0("FD", q))
  nm
}

#' Diversity of an assemblage from known identities
#'
#' Compute per-site Hill numbers on the taxonomic, and (if supplied)
#' phylogenetic and functional axes, from *known* species labels. This is the
#' validation target for [propagate_diversity()] and is also useful on its own.
#'
#' @param labels Character/factor vector of species labels, one per individual.
#' @param site Vector the same length as `labels` giving each individual's site.
#' @param tree Optional ape `phylo` (tip labels = species) for the phylogenetic
#'   axis.
#' @param trait_dist Optional functional distance matrix (dimnames = species)
#'   for the functional axis.
#' @param q Integer vector of orders (default `c(0, 1, 2)`).
#' @param species Character vector of the full species set (default the sorted
#'   unique `labels`).
#' @return A data frame, one row per site: `site`, `n` (individuals), and one
#'   column per axis-order metric (e.g. `TD1`, `PD1`, `FD1`).
#' @examples
#' lab <- c("a", "a", "b", "b", "c", "a")
#' site <- c(1, 1, 1, 2, 2, 2)
#' diversity_from_labels(lab, site)
#' @export
diversity_from_labels <- function(labels, site, tree = NULL, trait_dist = NULL,
                                  q = c(0, 1, 2),
                                  species = sort(unique(as.character(labels)))) {
  labels <- as.character(labels)
  axis_fn <- .make_axis_fn(species, q, tree, trait_dist)
  nm <- .axis_names(q, tree, trait_dist)
  sites <- sort(unique(site))
  rows <- lapply(sites, function(s) {
    idx <- which(site == s)
    cnt <- table(factor(labels[idx], levels = species))
    vals <- axis_fn(cnt)
    data.frame(site = s, n = length(idx), t(vals[nm]),
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Propagate classification uncertainty into per-site diversity
#'
#' The core of the framework. For each site, draw an identity for every
#' individual from its calibrated probability vector, reconstitute the
#' assemblage, and recompute every diversity axis; repeat `n_draws` times. The
#' mean across draws is the point estimate and the requested quantiles are the
#' confidence interval, so a site resting on confident identifications is
#' distinguished from one resting on coin-flips. The taxonomic axis is always
#' returned; the phylogenetic and functional axes are returned when `tree` and
#' `trait_dist` are supplied.
#'
#' The engine is agnostic to how `probs` was produced. For the leaf-spectra case
#' see [fit_spectral_classifier()] / [predict_probs()] and [fit_temperature()].
#'
#' @param probs Numeric probability matrix, rows = individuals, columns =
#'   species. Column names are the species labels and must match the tree tips
#'   and distance-matrix dimnames when those are supplied. Each row should sum
#'   to ~1 (see [scores_to_prob()], [temperature_scale()]).
#' @param site Vector of length `nrow(probs)` giving each individual's site.
#' @param tree Optional ape `phylo` (tip labels = species) for the phylogenetic
#'   axis.
#' @param trait_dist Optional functional distance matrix (dimnames = species)
#'   for the functional axis (see [trait_distance()]).
#' @param truth Optional vector of known labels (length `nrow(probs)`); when
#'   given, `true_*` columns are added for validation.
#' @param q Integer vector of orders (default `c(0, 1, 2)`).
#' @param n_draws Number of Monte Carlo draws (default 500).
#' @param ci Length-2 vector of lower/upper interval probabilities
#'   (default `c(0.025, 0.975)`).
#' @param min_n Drop sites with fewer than this many individuals (default 1).
#' @param seed Optional integer seed for reproducible draws.
#' @return A data frame, one row per retained site: `site`, `n_seedlings`, and
#'   for each metric a mean column plus `_lo`/`_hi` interval columns (and
#'   `true_*` columns if `truth` was given). Metric names are `TD`/`PD`/`FD`
#'   followed by the order, e.g. `PD1`, `PD1_lo`, `PD1_hi`, `true_PD1`.
#' @seealso [ses_structure()] for assemblage-structure effect sizes.
#' @examples
#' set.seed(1)
#' sp <- c("a", "b", "c")
#' P <- matrix(runif(30), 10, 3); P <- P / rowSums(P); colnames(P) <- sp
#' site <- rep(1:2, each = 5)
#' propagate_diversity(P, site, n_draws = 50, seed = 1)
#' @export
propagate_diversity <- function(probs, site, tree = NULL, trait_dist = NULL,
                                truth = NULL, q = c(0, 1, 2), n_draws = 500,
                                ci = c(0.025, 0.975), min_n = 1, seed = NULL) {
  probs <- as.matrix(probs)
  species <- colnames(probs)
  if (is.null(species)) stop("`probs` must have species column names.")
  if (length(site) != nrow(probs))
    stop("`site` must have length nrow(probs).")
  if (!is.null(truth) && length(truth) != nrow(probs))
    stop("`truth` must have length nrow(probs).")
  if (!is.null(tree) && !all(species %in% tree$tip.label))
    stop("Every species in colnames(probs) must be a tip of `tree`.")
  if (!is.null(trait_dist) && !all(species %in% rownames(trait_dist)))
    stop("Every species in colnames(probs) must be in `trait_dist` dimnames.")
  if (!is.null(seed)) set.seed(seed)

  axis_fn <- .make_axis_fn(species, q, tree, trait_dist)
  nm <- .axis_names(q, tree, trait_dist)
  qci <- function(x) c(mean(x), stats::quantile(x, ci[1]), stats::quantile(x, ci[2]))

  sites <- sort(unique(site))
  rows <- list()
  for (s in sites) {
    idx <- which(site == s); n <- length(idx)
    if (n < min_n) next
    Ppl <- probs[idx, , drop = FALSE]
    M <- matrix(NA_real_, n_draws, length(nm), dimnames = list(NULL, nm))
    for (i in seq_len(n_draws)) {
      drawn <- vapply(seq_len(n),
                      function(j) sample(species, 1L, prob = Ppl[j, ]), "")
      cnt <- table(factor(drawn, levels = species))
      M[i, ] <- axis_fn(cnt)[nm]
    }
    est <- apply(M, 2, qci)                    # 3 x metrics: mean, lo, hi
    r <- data.frame(site = s, n_seedlings = n,
                    stringsAsFactors = FALSE, check.names = FALSE)
    for (m in nm) r[[m]] <- est[1, m]
    for (m in nm) { r[[paste0(m, "_lo")]] <- est[2, m]
                    r[[paste0(m, "_hi")]] <- est[3, m] }
    if (!is.null(truth)) {
      tc <- table(factor(truth[idx], levels = species))
      tv <- axis_fn(tc)[nm]
      for (m in nm) r[[paste0("true_", m)]] <- tv[m]
    }
    rows[[length(rows) + 1]] <- r
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
